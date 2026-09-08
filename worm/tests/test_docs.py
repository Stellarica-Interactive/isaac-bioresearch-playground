"""Documentation cannot silently drift away from the data.

Free-text citations rot. These tests bind the prose to the manifests so that
adding a dataset without citing it, or citing a source id that no longer exists,
fails the build.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from common.data.registry import list_datasets
from worm.annotations.overlays import OVERLAYS
from worm.importers.sources import reference_totals, source_manifest

DOCS = Path(__file__).resolve().parents[2] / "docs"


@pytest.fixture(scope="module")
def references() -> str:
    return (DOCS / "references.md").read_text(encoding="utf-8")


def test_every_source_id_is_cited(references: str) -> None:
    missing = [sid for sid in source_manifest() if sid not in references]
    assert not missing, f"source ids in sources.toml but absent from references.md: {missing}"


def test_every_source_has_a_citation_and_licence() -> None:
    for sid, entry in source_manifest().items():
        assert entry.get("citation"), f"{sid} has no citation"
        assert entry.get("license"), f"{sid} has no licence"
        assert entry.get("url"), f"{sid} has no url"


def test_every_connectome_source_pins_a_hash() -> None:
    """An unpinned connectome could change under us without notice."""
    for sid, entry in source_manifest().items():
        if entry.get("kind") == "connectome":
            assert len(str(entry.get("sha256") or "")) == 64, f"{sid} has no SHA-256"


def test_every_registered_dataset_is_in_the_manifest() -> None:
    manifest = source_manifest()
    assert all(d in manifest for d in list_datasets())


def test_every_reference_total_names_its_oracle() -> None:
    for dataset_id, entry in reference_totals().items():
        assert entry.get("source"), f"{dataset_id} reference totals cite no source"


def test_every_overlay_source_resolves() -> None:
    manifest = source_manifest()
    for name, overlay in OVERLAYS.items():
        assert overlay.provenance().source_id in manifest, name


@pytest.mark.parametrize(
    "name",
    ["architecture.md", "biology.md", "model_assumptions.md", "references.md", "datasets.md"],
)
def test_expected_docs_exist_and_are_substantial(name: str) -> None:
    p = DOCS / name
    assert p.exists(), f"{name} is missing"
    assert len(p.read_text(encoding="utf-8")) > 2000, f"{name} is a stub"


def test_head_only_limitation_is_documented() -> None:
    """The single most consequential fact about the current data."""
    text = (DOCS / "datasets.md").read_text(encoding="utf-8")
    assert "ventral" in text.lower()
    assert "cannot" in text.lower()
