"""Verification of the Witvliet et al. 2021 importer against published totals.

This is the milestone's core gate. The reference numbers in
``worm/data/reference_totals.toml`` come from an independent analysis of the same
source spreadsheets (the OpenWorm Connectome Toolbox factsheets), so agreement is
meaningful: the oracle did not come from our parser.

Tests that need the raw ``.xlsx`` are marked ``requires_raw`` and skip on a fresh
clone. Everything else runs against the committed normalized data, so a checkout
with no network still verifies the format, the naming rules and the graph
invariants.
"""

from __future__ import annotations

import pytest

from common.data.registry import list_datasets
from common.data.schemas import CellCategory, Connectome, SynapseType
from common.data.validate import validate
from worm.importers.sources import (
    compare_to_reference,
    normalized_path,
    raw_path,
    reference_totals,
)
from worm.loader import load_anatomy, parse_raw

ADULTS = ["witvliet_2021_7", "witvliet_2021_8"]
ALL_WITVLIET = [f"witvliet_2021_{n}" for n in range(1, 9)]

#: Witvliet et al. reconstructed the brain, so the ventral cord is absent.
VNC_MOTOR_NEURON_PREFIXES = ("VA", "VB", "VC", "VD", "DA", "DB", "DD", "AS")

#: Exactly the anterior eight segments of each of the four muscle quadrants.
EXPECTED_MUSCLES_7 = frozenset(
    f"M{q}{i:02d}" for q in ("DL", "DR", "VL", "VR") for i in range(1, 9)
)


def committed(dataset_id: str) -> Connectome:
    if not (normalized_path(dataset_id) / "meta.json").exists():
        pytest.skip(f"{dataset_id} has not been built; run tools/build_datasets.py")
    return load_anatomy(dataset_id, prefer="normalized")


@pytest.fixture(scope="module")
def w7() -> Connectome:
    return committed("witvliet_2021_7")


class TestRegistration:
    def test_all_eight_datasets_registered(self) -> None:
        assert set(ALL_WITVLIET) <= set(list_datasets())

    def test_every_dataset_has_a_published_oracle(self) -> None:
        """A dataset we cannot check against anything is a dataset we cannot trust."""
        assert set(ALL_WITVLIET) <= set(reference_totals())


@pytest.mark.parametrize("dataset_id", ALL_WITVLIET)
class TestPublishedTotals:
    def test_committed_data_matches_published_figures(self, dataset_id: str) -> None:
        c = committed(dataset_id)
        mismatches = compare_to_reference(dataset_id, c.totals().to_dict())
        assert not mismatches, "\n".join(str(m) for m in mismatches)

    def test_committed_data_is_structurally_valid(self, dataset_id: str) -> None:
        assert validate(committed(dataset_id)) == []

    @pytest.mark.requires_raw
    def test_reparsing_raw_reproduces_committed_data(self, dataset_id: str) -> None:
        """Committed derived data must still be what the importer produces."""
        if not raw_path(dataset_id).exists():
            pytest.skip("raw source not fetched; run tools/fetch_datasets.py")
        fresh = parse_raw(dataset_id)
        stored = committed(dataset_id)
        assert fresh.totals() == stored.totals()
        assert fresh.cell_ids() == stored.cell_ids()
        assert {(e.pre, e.post, e.synapse_type, e.weight) for e in fresh.connections} == {
            (e.pre, e.post, e.synapse_type, e.weight) for e in stored.connections
        }


class TestElectricalConvention:
    """The single most error-prone convention in the project.

    Gap junctions are undirected and stored once. Published adjacency-matrix
    figures count each non-self junction twice and each self junction once, so:

        directed_entries = 2 * (stored - self_loops) + self_loops
        directed_weight  = 2 *  stored_weight - self_loop_weight

    Note the asymmetry in the second formula. It is why Witvliet #8's published
    electrical weight is 851 — an odd number that a naive doubling can never
    produce. That dataset contains a weight-2 gap junction from AIZR to itself.
    """

    def test_stored_once_in_canonical_orientation(self, w7: Connectome) -> None:
        assert all(e.pre <= e.post for e in w7.electrical())

    def test_directed_view_reproduces_published_entry_count(self, w7: Connectome) -> None:
        t = w7.totals()
        assert t.electrical_edges_undirected == 291
        assert t.electrical_self_loops == 6
        assert t.electrical_edges_directed == 2 * (291 - 6) + 6 == 576

    def test_directed_view_reproduces_published_weight(self, w7: Connectome) -> None:
        t = w7.totals()
        assert t.electrical_weight_undirected == 400
        assert t.electrical_weight_directed == 2 * 400 - 6 == 794

    def test_dataset_8_odd_weight_is_explained_by_a_weight_2_self_junction(self) -> None:
        c = committed("witvliet_2021_8")
        selfies = [e for e in c.electrical() if e.pre == e.post]
        assert len(selfies) == 8
        assert sum(e.weight for e in selfies) == 9
        assert [e.pre for e in selfies if e.weight == 2] == ["AIZR"]
        assert c.totals().electrical_weight_directed == 851

    def test_every_dataset_contains_self_junctions(self) -> None:
        """Systematic feature of the reconstruction, not a one-off artefact."""
        for dataset_id in ALL_WITVLIET:
            assert committed(dataset_id).totals().electrical_self_loops > 0, dataset_id


class TestAnatomicalScope:
    """Guards against silently using a head-only dataset for locomotion work."""

    def test_scope_is_recorded_as_head(self, w7: Connectome) -> None:
        assert str(w7.scope) == "head"

    def test_no_ventral_cord_motor_neurons(self, w7: Connectome) -> None:
        import re

        pattern = re.compile(rf"^({'|'.join(VNC_MOTOR_NEURON_PREFIXES)})\d+$")
        offenders = [c.id for c in w7.cells if pattern.match(c.id)]
        assert offenders == [], (
            f"Witvliet #7 should contain no ventral-cord motor neurons, found {offenders}. "
            "Either the wrong file was parsed, or this dataset is not what we think it is."
        )

    def test_muscles_are_exactly_the_anterior_eight_segments(self, w7: Connectome) -> None:
        muscles = {c.id for c in w7.cells if c.category is CellCategory.MUSCLE}
        assert muscles == EXPECTED_MUSCLES_7

    def test_provenance_states_the_limitation(self, w7: Connectome) -> None:
        assert "ventral-cord" in w7.provenance.notes


class TestGraphInvariants:
    def test_every_edge_endpoint_is_a_known_cell(self, w7: Connectome) -> None:
        ids = set(w7.cell_ids())
        for e in w7.connections:
            assert e.pre in ids and e.post in ids

    def test_no_duplicate_edges(self, w7: Connectome) -> None:
        keys = [e.key for e in w7.connections]
        assert len(keys) == len(set(keys))

    def test_all_weights_positive(self, w7: Connectome) -> None:
        assert all(e.weight >= 1 for e in w7.connections)

    def test_chemical_self_synapses_are_present_and_preserved(self, w7: Connectome) -> None:
        """11 neurons synapse onto themselves in the source. We do not drop them."""
        autapses = [e for e in w7.chemical() if e.pre == e.post]
        assert len(autapses) == 11

    def test_no_isolated_cells(self, w7: Connectome) -> None:
        touched = {e.pre for e in w7.connections} | {e.post for e in w7.connections}
        assert set(w7.cell_ids()) == touched

    def test_category_is_attributed_to_the_registry_not_the_connectome(
        self, w7: Connectome
    ) -> None:
        """EM shows a cell exists; it does not show that the cell is a neuron."""
        assert w7.cell("AVAL").field_sources["category"] != w7.provenance.source_id

    def test_anatomy_carries_no_synaptic_sign(self, w7: Connectome) -> None:
        assert all(str(e.sign) == "unknown" for e in w7.connections)

    def test_chemical_adjacency_is_not_symmetric(self, w7: Connectome) -> None:
        _, m = w7.adjacency(SynapseType.CHEMICAL)
        assert not (m == m.T).all()
