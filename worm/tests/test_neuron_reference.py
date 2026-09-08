"""The generated neuron reference must stay in sync and stay honest.

``docs/neurons.md`` is generated, so the important properties are that it is up to
date and that every curated claim in it carries a citation. Nothing in the
generator may paraphrase biology without a source.
"""

from __future__ import annotations

import csv
import re
from importlib.resources import files
from pathlib import Path

import pytest

from tools.generate_neuron_reference import CURATED, main, render
from worm.importers.naming import hermaphrodite_neuron_ids
from worm.loader import load

DOCS = Path(__file__).resolve().parents[2] / "docs"
NEURONS_MD = DOCS / "neurons.md"


def _descriptions() -> list[dict[str, str]]:
    text = (files("worm.data") / "annotations" / "cell_descriptions.csv").read_text(
        encoding="utf-8"
    )
    return list(csv.DictReader(text.splitlines()))


class TestGeneratedFileIsCurrent:
    def test_docs_neurons_md_is_up_to_date(self) -> None:
        assert main(["--check"]) == 0, (
            "docs/neurons.md is stale; run python tools/generate_neuron_reference.py"
        )

    def test_it_says_it_is_generated(self) -> None:
        assert "This file is generated" in NEURONS_MD.read_text(encoding="utf-8")


class TestCuratedNotes:
    def test_every_curated_note_cites_a_source(self) -> None:
        """No hand-written biological claim without a reference."""
        for klass, note in CURATED.items():
            assert note.source.strip(), f"{klass} has no source"
            assert re.search(r"\(\d{4}\)", note.source), (
                f"{klass}: source {note.source!r} has no year, so it is not a real citation"
            )
            assert note.summary.strip()

    def test_every_curated_class_is_a_real_neuron_class(self) -> None:
        """Against the whole animal, every curated note must name a cell that exists."""
        c, _ = load("cook_2019_herm")
        known = {x.class_name for x in c.cells if x.class_name}
        unknown = set(CURATED) - known
        assert not unknown, f"curated notes for classes that do not exist: {sorted(unknown)}"

    def test_absent_classes_are_flagged_rather_than_implied_present(self) -> None:
        """Rendered against a head-only dataset, missing cells must be marked missing.

        The committed reference uses the whole animal, where nothing is missing — so
        this exercises the mechanism directly rather than relying on it firing there.
        """
        head_only, _ = load("witvliet_2021_7")
        text = render(head_only, "witvliet_2021_7")
        assert "_not in this dataset_" in text
        assert "Not present in `witvliet_2021_7`" in text
        # PLM and PVD are posterior cells outside the Witvliet reconstruction volume.
        assert "`PLM`" in text.split("Not present in")[1][:200]


class TestDescriptionTable:
    def test_all_302_hermaphrodite_neurons_have_a_name_expansion(self) -> None:
        rows = {r["cell_id"]: r for r in _descriptions()}
        missing = sorted(
            n for n in hermaphrodite_neuron_ids() if not rows.get(n, {}).get("name_expansion")
        )
        assert not missing, f"neurons with no recorded name expansion: {missing}"

    def test_placeholder_text_is_not_stored(self) -> None:
        """WormAtlas writes '- To be added... -' for unrecorded fields; that is a gap."""
        for r in _descriptions():
            for field in ("name_expansion", "lineage", "classification"):
                assert "To be added" not in r[field], (r["cell_id"], field)

    def test_every_row_cites_a_source(self) -> None:
        assert all(r["source_ref"] for r in _descriptions())

    @pytest.mark.parametrize(
        ("cell", "fragment"),
        [
            ("ADEL", "Anterior DEirid"),
            ("AVAL", "Anterior Ventral Process A Left"),
            ("ASHL", "Amphid Single Cilium H Left"),
        ],
    )
    def test_known_name_expansions(self, cell: str, fragment: str) -> None:
        rows = {r["cell_id"]: r for r in _descriptions()}
        assert fragment in rows[cell]["name_expansion"]


@pytest.fixture(scope="module")
def text() -> str:
    return NEURONS_MD.read_text(encoding="utf-8")


class TestOutputContent:
    def test_every_neuron_in_the_dataset_appears(self, text: str) -> None:
        c, _ = load("witvliet_2021_7")
        missing = [n.id for n in c.neurons() if f"| `{n.id}` |" not in text]
        assert not missing, f"neurons absent from the reference: {missing}"

    def test_muscles_and_glia_are_not_filtered_out(self, text: str) -> None:
        assert "| `MDL01` |" in text
        assert "| `CEPshDL` |" in text

    def test_it_records_that_sign_is_unavailable(self, text: str) -> None:
        assert "measured" in text.lower()

    def test_evidence_marks_are_explained(self, text: str) -> None:
        assert "dim and variable" in text
        assert "taken up from neighbouring cells" in text
