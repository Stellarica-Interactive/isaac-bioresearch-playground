"""Canonical *C. elegans* cell naming.

Published datasets do not agree on how to spell a cell. Witvliet et al. write a
body wall muscle as ``BWM-VL01``; WormAtlas and Cook et al. write the same cell
as ``MVL01``. Witvliet write the excretory gland as ``excgl``; elsewhere it is
``exc_gl``. Ventral-cord motor neurons appear as both ``VB1`` and ``VB01``.

Every such rule lives here and nowhere else, because importers *and* annotation
overlays must agree on identity. If they drift, an annotation silently fails to
attach and the graph quietly loses a whole class of cells.

The rules below are not guesses. They were derived by enumerating every cell
label appearing in all eight Witvliet datasets and checking which ones failed to
resolve against the canonical registry in ``worm/data/cells.csv``; exactly two
patterns did (``BWM-*`` and ``excgl``). See ``docs/model_assumptions.md``.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files

from common.data.schemas import CellCategory

#: Cell identifiers with no systematic rule, mapped one by one.
KNOWN_ALIASES: dict[str, str] = {
    # Witvliet et al. 2021 spelling of the excretory gland cell.
    "excgl": "exc_gl",
    # Cook et al. 2019 lower-cases the posterior pharyngeal glial cell.
    "g1p": "g1P",
}

#: Witvliet body wall muscle labels: ``BWM-<quadrant><index>`` -> ``M<quadrant><index>``.
#: The two-digit index is preserved: ``MDL01`` is the canonical spelling, not ``MDL1``.
_BWM_RE = re.compile(r"^BWM-([DV][LR])(\d{2})$")

#: Cook et al. 2019 body wall muscle labels: ``dBWML1`` -> ``MDL01``, ``vBWMR23`` -> ``MVR23``.
#: A third convention for the same 95 cells; the index is one or two digits and is
#: zero-padded on the way in.
_COOK_BWM_RE = re.compile(r"^([dv])BWM([LR])(\d{1,2})$")

#: Optional leading zero on a neuron index, e.g. ``VB01`` for ``VB1``. Applied only
#: when stripping it produces a *known* cell and the original is unknown, so that
#: muscle names such as ``MDL01`` — where the zero is canonical — are never touched.
_LEADING_ZERO_RE = re.compile(r"^([A-Za-z]+)0(\d)$")


@dataclass(frozen=True, slots=True)
class CellRecord:
    """One row of the canonical cell registry."""

    cell_id: str
    category: CellCategory
    type_label: str
    sex_specific: str
    source_ref: str

    @property
    def is_male_specific(self) -> bool:
        return self.sex_specific == "male"


@lru_cache(maxsize=1)
def cell_registry() -> dict[str, CellRecord]:
    """Canonical cell registry, loaded from ``worm/data/cells.csv``.

    Regenerate with ``python tools/extract_annotations.py``.
    """
    text = (files("worm.data") / "cells.csv").read_text(encoding="utf-8")
    out: dict[str, CellRecord] = {}
    for row in csv.DictReader(text.splitlines()):
        rec = CellRecord(
            cell_id=row["cell_id"],
            category=CellCategory(row["category"]),
            type_label=row["type_label"],
            sex_specific=row["sex_specific"],
            source_ref=row["source_ref"],
        )
        out[rec.cell_id] = rec
    return out


def canonical_cell_id(raw: str) -> str:
    """Map a source-file cell label onto our canonical identifier.

    Idempotent: ``canonical_cell_id(canonical_cell_id(x)) == canonical_cell_id(x)``.
    An unrecognised label is returned unchanged rather than mangled, so that an
    importer reports "unknown cell" instead of silently merging it into another.
    """
    name = raw.strip()
    if not name:
        raise ValueError("empty cell label")

    if name in KNOWN_ALIASES:
        return KNOWN_ALIASES[name]

    m = _BWM_RE.match(name)
    if m:
        return f"M{m.group(1)}{m.group(2)}"

    m = _COOK_BWM_RE.match(name)
    if m:
        side = "D" if m.group(1) == "d" else "V"
        return f"M{side}{m.group(2)}{int(m.group(3)):02d}"

    registry = cell_registry()
    if name in registry:
        return name

    m = _LEADING_ZERO_RE.match(name)
    if m:
        stripped = f"{m.group(1)}{m.group(2)}"
        if stripped in registry:
            return stripped

    return name


def is_known_cell(cell_id: str) -> bool:
    return cell_id in cell_registry()


def category_of(cell_id: str) -> CellCategory:
    """Category of a *canonical* cell id.

    Raises :class:`KeyError` for unknown cells. Importers must handle that
    explicitly rather than defaulting to ``OTHER``: an unexpected cell name means
    either a new naming convention or the wrong file, and both deserve a failure.
    """
    try:
        return cell_registry()[cell_id].category
    except KeyError:
        raise KeyError(
            f"{cell_id!r} is not in the canonical cell registry (worm/data/cells.csv). "
            "Either it needs a canonicalization rule in worm/importers/naming.py, "
            "or the registry needs regenerating with tools/extract_annotations.py."
        ) from None


def type_label_of(cell_id: str) -> str:
    return cell_registry()[cell_id].type_label


def known_cell_ids() -> frozenset[str]:
    return frozenset(cell_registry())


def hermaphrodite_neuron_ids() -> frozenset[str]:
    """The 302 neurons of the adult hermaphrodite."""
    return frozenset(
        r.cell_id
        for r in cell_registry().values()
        if r.category is CellCategory.NEURON and not r.is_male_specific
    )


def body_wall_muscle_ids() -> frozenset[str]:
    """The 95 body wall muscles (head muscles plus main body muscles)."""
    return frozenset(
        r.cell_id
        for r in cell_registry().values()
        if r.type_label in ("Head muscle", "Main body muscle")
    )
