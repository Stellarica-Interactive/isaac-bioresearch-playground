"""Dataset registry.

An importer turns one published source file into an unannotated
:class:`~common.data.schemas.Connectome`. Registering it under a stable
``dataset_id`` is what lets every downstream tool — and eventually the simulator —
switch datasets without any code change.

Importers produce **anatomy only**. They must not attach functional roles,
neurotransmitters or synaptic signs; those come from separate publications and
are applied by :mod:`common.data.overlay`.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable

from common.data.schemas import Connectome


@dataclass(frozen=True, slots=True)
class SourceFile:
    """One downloadable published file an importer needs."""

    key: str
    """Key used to look the file up in the mapping passed to :meth:`ConnectomeImporter.parse`."""

    url: str
    filename: str
    sha256: str
    license: str
    citation: str
    doi: str | None = None
    notes: str = ""


@runtime_checkable
class ConnectomeImporter(Protocol):
    """Parses one published dataset into our normalized schema."""

    @property
    def dataset_id(self) -> str:
        """Stable key this importer is registered under."""
        ...

    def source_files(self) -> list[SourceFile]:
        """Files this importer needs, in ``sources.toml`` terms."""
        ...

    def parse(self, files: Mapping[str, Path]) -> Connectome:
        """Build an unannotated connectome from the given local file paths."""
        ...


_REGISTRY: dict[str, ConnectomeImporter] = {}


def register_importer(importer: ConnectomeImporter) -> ConnectomeImporter:
    if importer.dataset_id in _REGISTRY:
        raise ValueError(f"dataset id already registered: {importer.dataset_id!r}")
    _REGISTRY[importer.dataset_id] = importer
    return importer


def get_importer(dataset_id: str) -> ConnectomeImporter:
    try:
        return _REGISTRY[dataset_id]
    except KeyError:
        raise KeyError(
            f"unknown dataset {dataset_id!r}; known: {sorted(_REGISTRY)}"
        ) from None


def list_datasets() -> list[str]:
    return sorted(_REGISTRY)


def clear_registry() -> None:
    """Test helper. Not used in normal operation."""
    _REGISTRY.clear()
