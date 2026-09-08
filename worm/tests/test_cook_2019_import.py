"""Verification of the Cook et al. 2019 importer against published totals.

The whole-animal connectome, and the first dataset here that could support a
locomotion model. Two things get particular attention:

* the **ventral nerve cord is present** — the thing Witvliet lacks;
* the **gap-junction sheet choice**, because Cook ships two redundant sheets and
  reading the wrong one would double every gap junction weight.
"""

from __future__ import annotations

import re

import pytest

from common.data.registry import list_datasets
from common.data.schemas import CellCategory, Connectome, Scope, SIMRole
from common.data.validate import validate
from worm.importers.naming import body_wall_muscle_ids, canonical_cell_id, cell_registry
from worm.importers.sources import compare_to_reference, normalized_path, raw_path
from worm.loader import load, load_anatomy, parse_raw

DATASET = "cook_2019_herm"

#: The ventral-cord motor neuron classes that make the worm crawl.
VNC_CLASSES = ("VA", "VB", "VC", "VD", "DA", "DB", "DD", "AS")


def committed() -> Connectome:
    if not (normalized_path(DATASET) / "meta.json").exists():
        pytest.skip(f"{DATASET} has not been built; run tools/build_datasets.py")
    return load_anatomy(DATASET, prefer="normalized")


@pytest.fixture(scope="module")
def cook() -> Connectome:
    return committed()


class TestRegistration:
    def test_registered(self) -> None:
        assert DATASET in list_datasets()


class TestPublishedTotals:
    def test_matches_published_figures(self, cook: Connectome) -> None:
        mismatches = compare_to_reference(DATASET, cook.totals().to_dict())
        assert not mismatches, "\n".join(str(m) for m in mismatches)

    def test_structurally_valid(self, cook: Connectome) -> None:
        assert validate(cook) == []

    @pytest.mark.requires_raw
    def test_reparsing_raw_reproduces_committed_data(self, cook: Connectome) -> None:
        if not raw_path(DATASET).exists():
            pytest.skip("raw source not fetched; run tools/fetch_datasets.py")
        fresh = parse_raw(DATASET)
        assert fresh.totals() == cook.totals()
        assert fresh.cell_ids() == cook.cell_ids()
        assert {(e.pre, e.post, e.synapse_type, e.weight) for e in fresh.connections} == {
            (e.pre, e.post, e.synapse_type, e.weight) for e in cook.connections
        }


class TestGapJunctionSheetChoice:
    """Cook ships a mirrored sheet and an upper-triangle sheet describing one network.

    Reading the mirrored one into our once-per-pair storage would double every
    weight and still look plausible, so these assertions are the guard.
    """

    def test_stored_once_in_canonical_orientation(self, cook: Connectome) -> None:
        assert all(e.pre <= e.post for e in cook.electrical())

    def test_upper_triangle_counts(self, cook: Connectome) -> None:
        t = cook.totals()
        assert t.electrical_edges_undirected == 1450
        assert t.electrical_weight_undirected == 11680

    def test_odd_published_totals_are_explained_by_self_junctions(
        self, cook: Connectome
    ) -> None:
        """2883 and 23313 are both odd; a naive doubling can never produce them."""
        t = cook.totals()
        assert t.electrical_self_loops == 17
        assert sum(e.weight for e in cook.electrical() if e.pre == e.post) == 47
        assert t.electrical_edges_directed == 2 * (1450 - 17) + 17 == 2883
        assert t.electrical_weight_directed == 2 * 11680 - 47 == 23313

    def test_gap_junction_weights_are_not_doubled(self, cook: Connectome) -> None:
        """A canary against silently reading the mirrored sheet instead."""
        assert cook.totals().electrical_weight_undirected < 12000


class TestWholeAnimalScope:
    """The reason this dataset exists in the repository at all."""

    def test_scope_is_whole_animal(self, cook: Connectome) -> None:
        assert cook.scope is Scope.WHOLE_ANIMAL

    @pytest.mark.parametrize("klass", VNC_CLASSES)
    def test_ventral_cord_motor_neurons_are_present(
        self, cook: Connectome, klass: str
    ) -> None:
        pattern = re.compile(rf"^{klass}\d+$")
        members = [c.id for c in cook.cells if pattern.match(c.id)]
        assert members, f"no {klass} motor neurons found; this dataset cannot drive locomotion"

    def test_all_95_body_wall_muscles_are_present(self, cook: Connectome) -> None:
        muscles = {c.id for c in cook.cells if c.category is CellCategory.MUSCLE}
        assert body_wall_muscle_ids() <= muscles

    def test_all_302_hermaphrodite_neurons_are_present(self, cook: Connectome) -> None:
        assert len(cook.neurons()) == 302

    def test_motor_neurons_actually_reach_muscle(self, cook: Connectome) -> None:
        """Whole-animal means the neuromuscular junctions are in the data."""
        annotated, _ = load(DATASET)
        muscles = {c.id for c in annotated.cells if c.category is CellCategory.MUSCLE}
        motor = {c.id for c in annotated.cells_with_role(SIMRole.MOTOR)}
        nmj = [e for e in annotated.chemical() if e.pre in motor and e.post in muscles]
        assert len(nmj) > 500, f"only {len(nmj)} motor-to-muscle synapses"


class TestNaming:
    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            ("dBWML1", "MDL01"),
            ("dBWMR24", "MDR24"),
            ("vBWML7", "MVL07"),
            ("vBWMR23", "MVR23"),
            ("g1p", "g1P"),
        ],
    )
    def test_cook_specific_conventions(self, raw: str, expected: str) -> None:
        """A third muscle naming convention for the same 95 cells."""
        assert canonical_cell_id(raw) == expected

    def test_source_labels_record_the_original_spelling(self, cook: Connectome) -> None:
        assert "dBWML1" in cook.cell("MDL01").source_labels

    def test_no_two_source_labels_collapsed_onto_one_cell(self, cook: Connectome) -> None:
        """Silent merging would sum two cells' synapse counts and look fine."""
        for cell in cook.cells:
            assert len(cell.source_labels) == 1, (cell.id, cell.source_labels)


class TestGraphInvariants:
    def test_every_edge_endpoint_exists(self, cook: Connectome) -> None:
        ids = set(cook.cell_ids())
        assert all(e.pre in ids and e.post in ids for e in cook.connections)

    def test_no_duplicate_edges(self, cook: Connectome) -> None:
        keys = [e.key for e in cook.connections]
        assert len(keys) == len(set(keys))

    def test_chemical_autapses_are_preserved(self, cook: Connectome) -> None:
        assert len([e for e in cook.chemical() if e.pre == e.post]) == 38

    def test_anatomy_carries_no_synaptic_sign(self, cook: Connectome) -> None:
        assert all(str(e.sign) == "unknown" for e in cook.connections)


class TestAnnotationCoverage:
    def test_every_neuron_gets_annotated(self) -> None:
        """The overlays are whole-animal tables, so a whole-animal dataset is fully covered."""
        _, reports = load(DATASET)
        for r in reports:
            assert r.coverage == 1.0, f"{r.overlay_id}: unmatched {r.unmatched}"

    def test_leftover_table_entries_are_exactly_the_male_specific_cells(self) -> None:
        """A whole-animal *hermaphrodite* dataset should leave only male cells unused.

        The role table covers both sexes, so ray neurons (``R8AR``), spicule neurons
        (``SPCL``) and the male ventral-cord ``CA``/``CP`` classes are expected to go
        unmatched here. Anything else left over would mean a naming mismatch quietly
        preventing annotations from attaching.
        """
        _, reports = load(DATASET)
        registry = cell_registry()
        for r in reports:
            non_male = [
                cell_id
                for cell_id in r.unknown_in_source
                if cell_id in registry and not registry[cell_id].is_male_specific
            ]
            assert not non_male, (
                f"{r.overlay_id}: hermaphrodite cells in the annotation table that did "
                f"not attach to this dataset: {non_male}"
            )
