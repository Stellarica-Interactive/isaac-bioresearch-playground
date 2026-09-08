"""Schema invariants, tested on hand-built fixtures with no files involved."""

from __future__ import annotations

import pytest

from common.data.schemas import (
    Cell,
    CellCategory,
    Confidence,
    Connection,
    Connectome,
    Provenance,
    Scope,
    Sign,
    SIMRole,
    SynapseType,
)
from common.data.validate import validate

PROV = Provenance(
    source_id="test",
    kind="connectome",
    citation="Test fixture",
    url="https://example.invalid/test",
    license="CC0",
)


def make(cells: list[Cell], conns: list[Connection]) -> Connectome:
    return Connectome(
        id="test",
        organism="Test organism",
        sex="hermaphrodite",
        stage="adult",
        scope=Scope.WHOLE_ANIMAL,
        cells=tuple(cells),
        connections=tuple(conns),
        provenance=PROV,
    )


def neuron(name: str) -> Cell:
    return Cell(id=name, category=CellCategory.NEURON)


@pytest.fixture
def tiny() -> Connectome:
    cells = [neuron("A"), neuron("B"), neuron("C"), Cell(id="M1", category=CellCategory.MUSCLE)]
    conns = [
        Connection("A", "B", SynapseType.CHEMICAL, 3),
        Connection("B", "C", SynapseType.CHEMICAL, 5),
        Connection("C", "M1", SynapseType.CHEMICAL, 2),
        Connection("A", "C", SynapseType.ELECTRICAL, 4),  # canonical: A <= C
        Connection("B", "B", SynapseType.ELECTRICAL, 1),  # self junction
    ]
    return make(cells, conns)


class TestTotals:
    def test_cell_counts(self, tiny: Connectome) -> None:
        t = tiny.totals()
        assert (t.cells, t.neurons, t.muscles) == (4, 3, 1)

    def test_chemical(self, tiny: Connectome) -> None:
        t = tiny.totals()
        assert t.chemical_edges == 3
        assert t.chemical_weight == 10

    def test_electrical_stored_once(self, tiny: Connectome) -> None:
        t = tiny.totals()
        assert t.electrical_edges_undirected == 2
        assert t.electrical_weight_undirected == 5
        assert t.electrical_self_loops == 1

    def test_electrical_directed_view_matches_matrix_convention(self, tiny: Connectome) -> None:
        """A self junction is written once on the diagonal; every other one twice.

        This is exactly the arithmetic that reconciles our storage with published
        adjacency-matrix figures, and getting it wrong is silent.
        """
        t = tiny.totals()
        assert t.electrical_edges_directed == 2 * (2 - 1) + 1 == 3
        assert t.electrical_weight_directed == 2 * 5 - 1 == 9

    def test_other_cells_is_glia_plus_other(self) -> None:
        c = make(
            [
                Cell(id="G1", category=CellCategory.GLIA),
                Cell(id="O1", category=CellCategory.OTHER),
                Cell(id="O2", category=CellCategory.OTHER),
            ],
            [],
        )
        assert c.totals().other_cells == 3


class TestAdjacency:
    def test_electrical_matrix_is_symmetric(self, tiny: Connectome) -> None:
        nodes, m = tiny.adjacency(SynapseType.ELECTRICAL)
        assert (m == m.T).all()

    def test_chemical_matrix_is_not_symmetrized(self, tiny: Connectome) -> None:
        nodes, m = tiny.adjacency(SynapseType.CHEMICAL)
        i, j = nodes.index("A"), nodes.index("B")
        assert m[i, j] == 3
        assert m[j, i] == 0

    def test_node_order_is_sorted_not_insertion_order(self) -> None:
        c = make([neuron("Z"), neuron("A"), neuron("M")], [])
        nodes, _ = c.adjacency(SynapseType.CHEMICAL)
        assert nodes == ["A", "M", "Z"]


class TestSubgraph:
    def test_one_hop(self, tiny: Connectome) -> None:
        sub = tiny.subgraph(["B"], hops=1)
        assert set(sub.cell_ids()) == {"A", "B", "C"}

    def test_edges_restricted_to_kept_cells(self, tiny: Connectome) -> None:
        sub = tiny.subgraph(["A"], hops=1)
        assert all(e.pre in sub.cell_ids() and e.post in sub.cell_ids() for e in sub.connections)

    def test_unknown_seed_raises(self, tiny: Connectome) -> None:
        with pytest.raises(KeyError):
            tiny.subgraph(["nope"])


class TestValidation:
    def test_clean_graph_has_no_violations(self, tiny: Connectome) -> None:
        assert validate(tiny) == []

    def test_dangling_endpoint(self) -> None:
        c = make([neuron("A")], [Connection("A", "GHOST", SynapseType.CHEMICAL, 1)])
        assert any(v.code == "dangling_post" for v in validate(c))

    def test_zero_weight_rejected(self) -> None:
        c = make([neuron("A"), neuron("B")], [Connection("A", "B", SynapseType.CHEMICAL, 0)])
        assert any(v.code == "bad_weight" for v in validate(c))

    def test_non_canonical_gap_junction_rejected(self) -> None:
        """Storing B->A instead of A->B would double-count under symmetrization."""
        c = make([neuron("A"), neuron("B")], [Connection("B", "A", SynapseType.ELECTRICAL, 1)])
        assert any(v.code == "electrical_not_canonical" for v in validate(c))

    def test_duplicate_edge_rejected(self) -> None:
        c = make(
            [neuron("A"), neuron("B")],
            [
                Connection("A", "B", SynapseType.CHEMICAL, 1),
                Connection("A", "B", SynapseType.CHEMICAL, 2),
            ],
        )
        assert any(v.code == "duplicate_edge" for v in validate(c))

    def test_a_measured_synaptic_sign_is_rejected(self) -> None:
        """No connectome measures polarity. Claiming otherwise is a bug, not a value."""
        c = make(
            [neuron("A"), neuron("B")],
            [
                Connection(
                    "A", "B", SynapseType.CHEMICAL, 1,
                    sign=Sign.EXCITATORY, sign_confidence=Confidence.MEASURED,
                )
            ],
        )
        assert any(v.code == "sign_claimed_measured" for v in validate(c))

    def test_unresolvable_field_source_rejected(self) -> None:
        c = make([Cell(id="A", category=CellCategory.NEURON, field_sources={"roles": "nope"})], [])
        assert any(v.code == "unresolved_source" for v in validate(c))


class TestCellSemantics:
    def test_roles_and_neurotransmitters_are_tuples(self) -> None:
        """Multi-valued by design: polymodal cells and co-transmission are real."""
        c = Cell(
            id="X",
            category=CellCategory.NEURON,
            # Deliberately a list: __post_init__ must coerce it, so that an
            # importer passing a mutable sequence cannot leave the Cell mutable.
            roles=[SIMRole.SENSORY, SIMRole.INTER],  # type: ignore[arg-type]
        )
        assert c.roles == (SIMRole.SENSORY, SIMRole.INTER)
        assert c.has_role(SIMRole.INTER)

    def test_field_sources_is_immutable(self) -> None:
        c = Cell(id="X", category=CellCategory.NEURON, field_sources={"roles": "s"})
        with pytest.raises(TypeError):
            c.field_sources["roles"] = "other"  # type: ignore[index]

    def test_provenance_round_trips(self) -> None:
        assert Provenance.from_dict(PROV.to_dict()) == PROV
