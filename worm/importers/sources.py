"""Access to the source manifest and published reference totals.

``worm/data/sources.toml`` is the single place that says where a published file
came from, under what licence, and what its bytes hash to. ``fetch``, the
importers, and the provenance recorded in every derived file all read from it, so
a citation cannot drift away from the data it describes.

``worm/data/reference_totals.toml`` holds numbers quoted from an *independent*
analysis of the same source files. Keeping them in data rather than in test code
means the same oracle backs both ``pytest`` and
``inspect_connectome.py ... verify``.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from pathlib import Path

from common.data.registry import SourceFile

RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"
NORMALIZED_DIR = Path(__file__).resolve().parent.parent / "data" / "normalized"


@lru_cache(maxsize=1)
def source_manifest() -> dict[str, dict[str, object]]:
    text = (files("worm.data") / "sources.toml").read_text(encoding="utf-8")
    return tomllib.loads(text)


@lru_cache(maxsize=1)
def reference_totals() -> dict[str, dict[str, object]]:
    text = (files("worm.data") / "reference_totals.toml").read_text(encoding="utf-8")
    return tomllib.loads(text)


def source_file(source_id: str, *, key: str = "connectome") -> SourceFile:
    try:
        entry = source_manifest()[source_id]
    except KeyError:
        raise KeyError(
            f"{source_id!r} is not in worm/data/sources.toml; known: "
            f"{sorted(source_manifest())}"
        ) from None
    return SourceFile(
        key=key,
        url=str(entry["url"]),
        filename=str(entry["filename"]),
        sha256=str(entry.get("sha256") or ""),
        license=str(entry["license"]),
        citation=str(entry["citation"]),
        doi=_opt(entry.get("doi")),
        notes=str(entry.get("mirror_note") or ""),
    )


def raw_path(source_id: str) -> Path:
    """Where :mod:`tools.fetch_datasets` puts the file for ``source_id``."""
    return RAW_DIR / source_file(source_id).filename


def normalized_path(dataset_id: str) -> Path:
    return NORMALIZED_DIR / dataset_id


@dataclass(frozen=True, slots=True)
class TotalsMismatch:
    field: str
    expected: int
    actual: int

    def __str__(self) -> str:  # pragma: no cover - trivial
        return f"{self.field}: published {self.expected}, parsed {self.actual}"


def compare_to_reference(dataset_id: str, observed: dict[str, int]) -> list[TotalsMismatch]:
    """Compare parsed totals against the published oracle.

    Only fields present in ``reference_totals.toml`` are checked; a dataset with
    no published figures returns an empty list, which the caller should treat as
    "unverified", not "verified".
    """
    ref = reference_totals().get(dataset_id)
    if not ref:
        return []
    out: list[TotalsMismatch] = []
    for field, expected in ref.items():
        if not isinstance(expected, int):
            continue  # e.g. the `source` URL
        if field not in observed:
            continue
        if observed[field] != expected:
            out.append(TotalsMismatch(field, expected, observed[field]))
    return out


def has_reference(dataset_id: str) -> bool:
    return dataset_id in reference_totals()


def _opt(v: object) -> str | None:
    return None if v is None or v == "" else str(v)
