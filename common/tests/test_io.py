"""Serialization determinism and round-tripping.

These are locked down before any real dataset is parsed, so that committed golden
files never need regenerating for a mere formatting change.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from common.data.io import (
    CELLS_FILE,
    CONNECTIONS_FILE,
    META_FILE,
    load_connectome,
    render,
    save_connectome,
)
from common.data.schemas import (
    Cell,
    CellCategory,
    Connection,
    Connectome,
    Provenance,
    Scope,
    SynapseType,
)

PROV = Provenance(
    source_id="test",
    kind="connectome",
    citation="Test fixture",
    url="https://example.invalid/test",
    license="CC0",
    sha256="deadbeef",
    retrieved="2026-01-01",
)


@pytest.fixture
def sample() -> Connectome:
    cells = [
        Cell(id="ZZZ", category=CellCategory.NEURON, field_sources={"category": "test"}),
        Cell(id="AAA", category=CellCategory.NEURON, field_sources={"category": "test"}),
        Cell(
            id="M1",
            category=CellCategory.MUSCLE,
            source_labels=("BWM-VL01", "MVL01"),
            field_sources={"category": "test"},
        ),
        Cell(id="G1", category=CellCategory.GLIA, field_sources={"category": "test"}),
        Cell(id="O1", category=CellCategory.OTHER, field_sources={"category": "test"}),
    ]
    conns = [
        Connection("ZZZ", "AAA", SynapseType.CHEMICAL, 2),
        Connection("AAA", "M1", SynapseType.CHEMICAL, 7),
        Connection("AAA", "ZZZ", SynapseType.ELECTRICAL, 3),
        Connection("G1", "G1", SynapseType.ELECTRICAL, 1),
    ]
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


class TestDeterminism:
    def test_rendering_twice_is_identical(self, sample: Connectome) -> None:
        assert render(sample) == render(sample)

    def test_only_lf_newlines(self, sample: Connectome) -> None:
        """A CRLF checkout would break every recorded hash. See .gitattributes."""
        for name, text in render(sample).items():
            assert "\r" not in text, name

    def test_rows_are_sorted_regardless_of_input_order(self, sample: Connectome) -> None:
        cells_csv = render(sample)[CELLS_FILE]
        ids = [line.split(",")[0] for line in cells_csv.splitlines()[1:]]
        assert ids == sorted(ids)

    def test_connections_sorted_by_type_then_endpoints(self, sample: Connectome) -> None:
        rows = render(sample)[CONNECTIONS_FILE].splitlines()[1:]
        keys = [(r.split(",")[2], r.split(",")[0], r.split(",")[1]) for r in rows]
        assert keys == sorted(keys)

    def test_meta_records_csv_hashes(self, sample: Connectome) -> None:
        files = render(sample)
        meta = json.loads(files[META_FILE])
        assert set(meta["files"]) == {CELLS_FILE, CONNECTIONS_FILE}
        assert all(len(v["sha256"]) == 64 for v in meta["files"].values())


class TestRoundTrip:
    def test_save_load_preserves_structure(self, sample: Connectome, tmp_path: Path) -> None:
        save_connectome(sample, tmp_path / "ds")
        back = load_connectome(tmp_path / "ds")
        assert back.cell_ids() == sample.cell_ids()
        assert {(e.pre, e.post, e.synapse_type, e.weight) for e in back.connections} == {
            (e.pre, e.post, e.synapse_type, e.weight) for e in sample.connections
        }
        assert back.totals() == sample.totals()
        assert back.provenance == sample.provenance
        assert back.scope is sample.scope

    def test_serialization_is_idempotent(self, sample: Connectome, tmp_path: Path) -> None:
        first = save_connectome(sample, tmp_path / "ds")
        second = render(load_connectome(tmp_path / "ds"))
        assert first == second

    def test_category_provenance_survives(self, sample: Connectome, tmp_path: Path) -> None:
        """Cell category is an annotation, not a measurement. Its citation must persist."""
        save_connectome(sample, tmp_path / "ds")
        back = load_connectome(tmp_path / "ds")
        assert back.cell("AAA").field_sources["category"] == "test"

    def test_source_labels_survive(self, sample: Connectome, tmp_path: Path) -> None:
        save_connectome(sample, tmp_path / "ds")
        assert load_connectome(tmp_path / "ds").cell("M1").source_labels == ("BWM-VL01", "MVL01")


class TestCorruptionDetection:
    def test_edited_csv_is_rejected(self, sample: Connectome, tmp_path: Path) -> None:
        dest = tmp_path / "ds"
        save_connectome(sample, dest)
        p = dest / CONNECTIONS_FILE
        p.write_text(p.read_text().replace(",7", ",99"), encoding="utf-8", newline="")
        with pytest.raises(ValueError, match="does not match the hash"):
            load_connectome(dest)

    def test_missing_meta_is_reported_clearly(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="not a normalized dataset"):
            load_connectome(tmp_path)

    def test_wrong_columns_rejected(self, sample: Connectome, tmp_path: Path) -> None:
        dest = tmp_path / "ds"
        save_connectome(sample, dest)
        p = dest / CELLS_FILE
        lines = p.read_text().splitlines()
        lines[0] = "id,category,bogus,source_labels"
        p.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="")
        with pytest.raises(ValueError):
            load_connectome(dest)
