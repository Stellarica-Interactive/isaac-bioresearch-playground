"""Annotation overlays: coverage, provenance, and the absence of silent defaults.

The purpose of these tests is not only "does the overlay work" but "does it stay
honest": an annotation must never be attributed to the connectome, and a cell the
table does not cover must be *reported* rather than quietly filled in.
"""

from __future__ import annotations

import csv
import re
from importlib.resources import files

import pytest

from common.data.schemas import (
    CellCategory,
    Confidence,
    Connectome,
    Neurotransmitter,
    Sign,
    SIMRole,
)
from worm.annotations.overlays import DEFAULT_OVERLAYS, OVERLAYS, get_overlays
from worm.importers.naming import body_wall_muscle_ids
from worm.loader import load, load_anatomy

DATASET = "witvliet_2021_7"

#: The neuromuscular tests need the ventral cord and all 95 body wall muscles, so
#: they use the whole-animal dataset rather than the head-only one.
WHOLE_ANIMAL = "cook_2019_herm"


def _table(name: str) -> list[dict[str, str]]:
    text = (files("worm.data") / "annotations" / name).read_text(encoding="utf-8")
    return list(csv.DictReader(text.splitlines()))


@pytest.fixture(scope="module")
def annotated() -> Connectome:
    c, _ = load(DATASET)
    return c


@pytest.fixture(scope="module")
def anatomy() -> Connectome:
    return load_anatomy(DATASET)


class TestSeparationOfConcerns:
    def test_stored_anatomy_has_no_annotations(self, anatomy: Connectome) -> None:
        """The committed connectome files must contain measurement only."""
        assert all(not c.roles for c in anatomy.cells)
        assert all(not c.neurotransmitters for c in anatomy.cells)
        assert all(c.class_name is None for c in anatomy.cells)

    def test_polarity_is_available_but_never_applied_by_default(self) -> None:
        """A predicted synapse sign must never arrive unless it was asked for.

        The overlay exists -- the neural runtime needs it -- but selecting it is
        always an explicit act, and the runtime reports its coverage when it does.
        """
        assert "polarity" in OVERLAYS
        assert "polarity" not in DEFAULT_OVERLAYS

    def test_anatomy_loaded_by_default_carries_no_sign(self) -> None:
        c, _ = load(DATASET)
        assert all(str(e.sign) == "unknown" for e in c.connections)

    def test_annotations_are_attributed_to_their_own_sources(
        self, annotated: Connectome
    ) -> None:
        cell = annotated.cell("AVAL")
        for field in ("roles", "neurotransmitters", "class_name"):
            assert cell.field_sources[field] != annotated.provenance.source_id

    def test_every_cited_source_resolves(self, annotated: Connectome) -> None:
        for cell in annotated.cells:
            for src in cell.field_sources.values():
                assert src in annotated.sources
                assert annotated.sources[src].citation

    def test_unknown_overlay_name_raises(self) -> None:
        with pytest.raises(KeyError):
            get_overlays(("nope",))


class TestCoverage:
    def test_all_reports_are_returned(self) -> None:
        _, reports = load(DATASET)
        assert [r.overlay_id for r in reports] == list(DEFAULT_OVERLAYS)

    def test_every_neuron_is_annotated(self) -> None:
        _, reports = load(DATASET)
        for r in reports:
            assert r.coverage == 1.0, f"{r.overlay_id}: unmatched {r.unmatched}"

    def test_whole_animal_tables_over_a_head_dataset_report_the_excess(self) -> None:
        """Expected and harmless in this direction; a red flag in the other."""
        _, reports = load(DATASET)
        assert all(r.unknown_in_source for r in reports)

    def test_non_neurons_are_left_alone(self, annotated: Connectome) -> None:
        muscles = [c for c in annotated.cells if c.category is CellCategory.MUSCLE]
        assert muscles
        assert all(not m.roles and not m.neurotransmitters for m in muscles)


class TestSIMRoles:
    def test_known_assignments(self, annotated: Connectome) -> None:
        assert SIMRole.SENSORY in annotated.cell("ASHL").roles
        assert SIMRole.INTER in annotated.cell("AVAL").roles
        assert SIMRole.MOTOR in annotated.cell("RMDL").roles

    def test_multi_role_cells_keep_both_roles(self) -> None:
        """The published label for some cells records disagreement between studies."""
        rows = _table("sim_roles.csv")
        by_cell: dict[str, list[str]] = {}
        for r in rows:
            by_cell.setdefault(r["cell_id"], []).append(r["role"])
        multi = {k: v for k, v in by_cell.items() if len(v) > 1}
        assert multi, "expected some cells to carry two roles"
        assert all("interneuron in White" in r["type_label"] or
                   "motorneuron in White" in r["type_label"]
                   for r in rows if r["cell_id"] in multi)

    def test_cells_we_decline_to_classify_are_marked_unknown_not_guessed(self) -> None:
        rows = _table("sim_roles.csv")
        unknown = sorted({r["cell_id"] for r in rows if r["role"] == "unknown"})
        assert unknown == ["CANL", "CANR", "MCL", "MCR", "MI", "NSML", "NSMR"]

    def test_every_row_cites_a_source(self) -> None:
        assert all(r["source_ref"] for r in _table("sim_roles.csv"))


class TestNeurotransmitters:
    def test_known_assignments(self, annotated: Connectome) -> None:
        assert annotated.cell("ASHL").neurotransmitters == (Neurotransmitter.GLUTAMATE,)
        assert annotated.cell("AVAL").neurotransmitters == (Neurotransmitter.ACETYLCHOLINE,)

    def test_all_302_hermaphrodite_neurons_are_covered(self) -> None:
        assert len({r["cell_id"] for r in _table("neurotransmitters.csv")}) == 302

    def test_orphan_neurons_are_a_result_not_a_blank(self) -> None:
        """16 neurons express no known transmitter pathway gene. That is a finding."""
        rows = _table("neurotransmitters.csv")
        orphans = [r for r in rows if r["neurotransmitter"] == "unknown"]
        assert len(orphans) == 16
        assert all(r["evidence"] == "orphan_no_pathway_gene_detected" for r in orphans)

    def test_uptake_is_distinguished_from_synthesis(self) -> None:
        """AVFL/R take up GABA from neighbours rather than making it."""
        rows = {r["cell_id"]: r for r in _table("neurotransmitters.csv")}
        assert rows["AVFL"]["evidence"] == "uptake_not_synthesis"
        assert rows["AVFL"]["neurotransmitter"] == "gaba"
        assert rows["CANL"]["evidence"] == "uptake_not_synthesis"

    def test_dim_variable_expression_is_distinguished(self) -> None:
        """Weaker evidence than a clean positive, and marked as such."""
        rows = {r["cell_id"]: r for r in _table("neurotransmitters.csv")}
        assert rows["AWAL"]["evidence"] == "reporter_expression_dim_variable"
        assert rows["AWAL"]["neurotransmitter"] == "acetylcholine"

    def test_every_row_keeps_its_source_row_number(self) -> None:
        """So a disputed assignment can be checked in the original spreadsheet."""
        assert all(r["source_row"].isdigit() for r in _table("neurotransmitters.csv"))


class TestNeuronClasses:
    def test_bilateral_pairs_share_a_class(self, annotated: Connectome) -> None:
        assert annotated.cell("AVAL").class_name == annotated.cell("AVAR").class_name == "AVA"

    def test_class_is_not_derivable_by_stripping_letters(self, annotated: Connectome) -> None:
        """RMDDL and RMDL are both class RMD - a rule no string operation would find."""
        assert annotated.cell("RMDDL").class_name == "RMD"
        assert annotated.cell("RMDL").class_name == "RMD"


@pytest.fixture(scope="module")
def signed() -> Connectome:
    c, _ = load(WHOLE_ANIMAL, annotations=("classes", "sim", "nt", "nmj"))
    return c


class TestNeuromuscularPolarity:
    """The overlay that unblocks locomotion, and the limits it keeps.

    These synapses are what actually move the animal. No connectome or prediction
    dataset signs them, so this overlay is the only thing standing between the
    connectome and a motor output.
    """

    def test_it_is_opt_in(self) -> None:
        assert "nmj" in OVERLAYS
        assert "nmj" not in DEFAULT_OVERLAYS

    def test_requesting_it_without_transmitters_is_refused(self) -> None:
        """It reads the presynaptic transmitter, so 'nt' is not optional."""
        with pytest.raises(ValueError, match="presynaptic transmitter"):
            get_overlays(("nmj",))

    def test_dependency_order_is_enforced_not_assumed(self) -> None:
        """Asking in the wrong order must still run 'nt' first, not annotate nothing."""
        ids = [o.overlay_id for o in get_overlays(("nmj", "nt"))]
        assert ids.index("nt") < ids.index("nmj")

    def test_cholinergic_motor_neurons_excite_muscle(self, signed: Connectome) -> None:
        edges = [e for e in signed.chemical() if e.pre == "VB7" and e.post.startswith("M")]
        assert edges
        assert all(e.sign is Sign.EXCITATORY for e in edges)

    def test_gabaergic_motor_neurons_inhibit_muscle(self, signed: Connectome) -> None:
        for cell in ("DD3", "VD5"):
            edges = [e for e in signed.chemical() if e.pre == cell and e.post.startswith("M")]
            assert edges, cell
            assert all(e.sign is Sign.INHIBITORY for e in edges), cell

    def test_signs_are_measured_physiology_not_prediction(self, signed: Connectome) -> None:
        """This is the one sign source in the project that is not a guess."""
        nmj = [
            e
            for e in signed.chemical()
            if e.post in body_wall_muscle_ids() and e.sign is not Sign.UNKNOWN
        ]
        assert nmj
        assert all(e.sign_confidence is Confidence.PUBLISHED_ANNOTATION for e in nmj)
        assert {e.field_sources["sign"] for e in nmj} == {
            "richmond_1999_nmj_receptors",
            "mcintire_1993_gaba_inhibitory",
        }

    def test_coverage_matches_the_documented_figure(self) -> None:
        _, reports = load(WHOLE_ANIMAL, annotations=("classes", "sim", "nt", "nmj"))
        report = next(r for r in reports if r.overlay_id == "nmj")
        assert (report.matched, report.eligible) == (892, 956)

    def test_non_body_wall_muscle_is_left_alone(self, signed: Connectome) -> None:
        """Pharyngeal pharmacology differs: glutamate is inhibitory there."""
        muscles = {x.id for x in signed.cells if x.category is CellCategory.MUSCLE}
        other = muscles - body_wall_muscle_ids()
        assert other
        edges = [e for e in signed.chemical() if e.post in other]
        assert edges
        assert all(e.sign is Sign.UNKNOWN for e in edges)

    def test_modulatory_transmitters_are_left_unsigned(self, signed: Connectome) -> None:
        """Dopamine acts through GPCRs; a reversal potential is the wrong model."""
        dopaminergic = {
            x.id for x in signed.cells if Neurotransmitter.DOPAMINE in x.neurotransmitters
        }
        edges = [
            e
            for e in signed.chemical()
            if e.pre in dopaminergic and e.post in body_wall_muscle_ids()
        ]
        assert edges
        assert all(e.sign is Sign.UNKNOWN for e in edges)

    def test_the_motor_classes_recovered_match_textbook_biology(
        self, signed: Connectome
    ) -> None:
        """A cross-check nothing in our code arranges.

        Which cells excite muscle and which inhibit it falls out of two independent
        datasets -- Cook's connectome and Wang's transmitter atlas -- that were never
        reconciled against each other. It lands on the textbook division exactly.
        """
        bwm = body_wall_muscle_ids()
        vnc = re.compile(r"^(VA|VB|VC|VD|DA|DB|DD|AS)\d+$")

        def classes(sign: Sign) -> set[str]:
            return {
                vnc.match(e.pre).group(1)  # type: ignore[union-attr]
                for e in signed.chemical()
                if e.post in bwm and e.sign is sign and vnc.match(e.pre)
            }

        assert classes(Sign.EXCITATORY) == {"AS", "DA", "DB", "VA", "VB", "VC"}
        assert classes(Sign.INHIBITORY) == {"DD", "VD"}


def test_polarity_table_names_match_the_connectome() -> None:
    """Every Fenyves row must name cells this dataset actually has.

    This is the test that was missing. Fenyves zero-pads the ventral cord motor
    neurons -- ``DB01`` where Cook writes ``DB1`` -- and because the extractor did
    not canonicalise, 634 of 3638 rows named cells that did not exist. 582 of them
    carried a real predicted sign, so a quarter of all polarity data was silently
    discarded, and the loss fell on the entire locomotor circuit: every AS, DA, DB,
    VA, VB, DD and VD cell.

    Nothing failed. The overlay reported it faithfully as "739 table entries not in
    this dataset", which is exactly what a table describing a *different animal*
    would also produce -- so the summary line was read as expected noise for months.

    Hence an assertion rather than a report. A naming mismatch and a genuine
    dataset difference look identical in prose and completely different in numbers.
    """
    connectome, _ = load("cook_2019_herm")
    known = {c.id for c in connectome.cells}

    rows = list(
        csv.DictReader(
            (files("worm.data.annotations") / "polarity_fenyves2020.csv").open(
                encoding="utf-8"
            )
        )
    )
    assert rows, "polarity table is empty"

    unknown = {
        name
        for row in rows
        for name in (row["pre"], row["post"])
        if name not in known
    }
    # Fenyves covers a slightly different cell set, so a few genuine absences are
    # expected; a systematic naming mismatch is not.
    assert len(unknown) <= 12, (
        f"{len(unknown)} cell names in the polarity table are not in the connectome. "
        f"That is too many to be a dataset difference -- check canonicalisation. "
        f"Sample: {sorted(unknown)[:10]}"
    )


def test_polarity_overlay_covers_most_of_the_connectome() -> None:
    """A floor on how much of the connectome carries a predicted sign.

    Coverage is 72.1%. It was 59.4% before the canonicalisation fix above, and a
    regression there would quietly shrink the simulated network again rather than
    fail anything.
    """
    _, reports = load("cook_2019_herm", annotations=("classes", "sim", "nt", "polarity"))
    polarity = next(r for r in reports if r.overlay_id == "polarity")
    assert polarity.matched / (polarity.matched + len(polarity.unmatched)) > 0.70
