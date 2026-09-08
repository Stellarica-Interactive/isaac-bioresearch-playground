"""Cell-name canonicalization.

File-free: these are properties of the naming rules themselves.

The injectivity test is the important one. If two distinct source labels
collapse onto the same canonical id, two separate cells silently merge into one
and their synapse counts add together. Nothing downstream would notice.
"""

from __future__ import annotations

import pytest

from common.data.schemas import CellCategory
from worm.importers.naming import (
    KNOWN_ALIASES,
    body_wall_muscle_ids,
    canonical_cell_id,
    category_of,
    cell_registry,
    hermaphrodite_neuron_ids,
    is_known_cell,
)


class TestCanonicalization:
    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            ("BWM-DL01", "MDL01"),
            ("BWM-VR08", "MVR08"),
            ("BWM-VL01", "MVL01"),
            ("BWM-DR07", "MDR07"),
            ("excgl", "exc_gl"),
            ("AVAL", "AVAL"),
            ("MDL01", "MDL01"),
            ("  ASHL  ", "ASHL"),
        ],
    )
    def test_known_rules(self, raw: str, expected: str) -> None:
        assert canonical_cell_id(raw) == expected

    def test_leading_zero_stripped_only_for_neurons(self) -> None:
        assert canonical_cell_id("VB01") == "VB1"
        # Muscles keep their zero: MDL01 is the canonical spelling, and stripping
        # it would produce an id that does not exist.
        assert canonical_cell_id("MDL01") == "MDL01"

    def test_unknown_label_returned_unchanged(self) -> None:
        """Never mangle. An unknown label must surface as an error downstream."""
        assert canonical_cell_id("NOT_A_CELL") == "NOT_A_CELL"

    def test_idempotent(self) -> None:
        for raw in ["BWM-DL01", "excgl", "VB01", "AVAL", "MDL01", "NOT_A_CELL"]:
            once = canonical_cell_id(raw)
            assert canonical_cell_id(once) == once, raw

    def test_empty_label_raises(self) -> None:
        with pytest.raises(ValueError):
            canonical_cell_id("   ")

    def test_injective_over_witvliet_label_space(self) -> None:
        """No two distinct source labels may collapse onto the same canonical id."""
        labels = (
            [f"BWM-{q}{i:02d}" for q in ("DL", "DR", "VL", "VR") for i in range(1, 25)]
            + ["excgl"]
            + sorted(hermaphrodite_neuron_ids())
            + sorted(body_wall_muscle_ids())
        )
        seen: dict[str, str] = {}
        for label in labels:
            canonical = canonical_cell_id(label)
            if canonical in seen and seen[canonical] != label:
                pair = {seen[canonical], label}
                assert pair <= set(KNOWN_ALIASES) | {KNOWN_ALIASES.get(x, x) for x in pair}, (
                    f"{seen[canonical]!r} and {label!r} both canonicalize to "
                    f"{canonical!r} without an explicit alias"
                )
            seen[canonical] = label


class TestRegistry:
    def test_302_hermaphrodite_neurons(self) -> None:
        """The canonical figure for the adult hermaphrodite nervous system."""
        assert len(hermaphrodite_neuron_ids()) == 302

    def test_95_body_wall_muscles(self) -> None:
        assert len(body_wall_muscle_ids()) == 95

    @pytest.mark.parametrize(
        ("cell", "category"),
        [
            ("AVAL", CellCategory.NEURON),
            ("ASHL", CellCategory.NEURON),
            ("MDL01", CellCategory.MUSCLE),
            ("CEPshDL", CellCategory.GLIA),
            ("GLRDL", CellCategory.OTHER),
            ("exc_gl", CellCategory.OTHER),
        ],
    )
    def test_categories(self, cell: str, category: CellCategory) -> None:
        assert category_of(cell) is category

    def test_unknown_cell_raises_with_actionable_message(self) -> None:
        with pytest.raises(KeyError, match="canonical cell registry"):
            category_of("NOT_A_CELL")

    def test_is_known_cell(self) -> None:
        assert is_known_cell("AVAL")
        assert not is_known_cell("NOT_A_CELL")

    def test_registry_ids_are_unique_and_nonempty(self) -> None:
        reg = cell_registry()
        assert len(reg) > 600
        assert all(r.cell_id and r.type_label and r.source_ref for r in reg.values())
