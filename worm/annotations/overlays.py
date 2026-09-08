"""Annotation overlays for *C. elegans*.

Each overlay attaches one kind of published annotation to a connectome that was
parsed from electron-microscopy data alone. See :mod:`common.data.overlay` for
why these are kept separate from the anatomy.

Three overlays exist today:

``classes``
    Anatomical neuron class (``AVAL`` -> ``AVA``), from Wang et al. 2024.

``sim``
    Coarse sensory / interneuron / motor role, from WormAtlas cell listings and
    the Cook et al. 2019 groupings.

``nt``
    Neurotransmitter usage, from the Wang et al. 2024 reporter-allele atlas.

A polarity overlay (which chemical synapses are excitatory or inhibitory) is
**not** included yet. That value is not measured anywhere; it can only be
predicted from receptor gene expression, and letting a prediction reach a
dynamics model without a deliberate opt-in is the most likely way this project
would end up fooling itself. See ``docs/model_assumptions.md``.
"""

from __future__ import annotations

import csv
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, replace
from functools import lru_cache
from importlib.resources import files

from common.data.overlay import AnnotationOverlay, OverlayReport
from common.data.schemas import (
    Cell,
    CellCategory,
    Connectome,
    Neurotransmitter,
    Provenance,
    SIMRole,
)
from worm.importers.sources import source_file


@lru_cache(maxsize=8)
def _rows(name: str) -> tuple[Mapping[str, str], ...]:
    text = (files("worm.data") / "annotations" / name).read_text(encoding="utf-8")
    return tuple(csv.DictReader(text.splitlines()))


def _provenance(source_id: str, note: str) -> Provenance:
    spec = source_file(source_id)
    return Provenance(
        source_id=source_id,
        kind="annotation",
        citation=spec.citation,
        url=spec.url,
        license=spec.license,
        doi=spec.doi,
        file_name=spec.filename or None,
        sha256=spec.sha256 or None,
        notes=note,
    )


def _report(
    overlay_id: str, annotated: set[str], connectome: Connectome, table_ids: set[str]
) -> OverlayReport:
    eligible = {c.id for c in connectome.cells if c.category is CellCategory.NEURON}
    return OverlayReport(
        overlay_id=overlay_id,
        matched=len(annotated),
        eligible=len(eligible),
        unmatched=tuple(sorted(eligible - annotated)),
        unknown_in_source=tuple(sorted(table_ids - {c.id for c in connectome.cells})),
    )


def _rebuild(connectome: Connectome, updated: Mapping[str, Cell]) -> Connectome:
    return replace(
        connectome, cells=tuple(updated.get(c.id, c) for c in connectome.cells)
    )


# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NeuronClassOverlay:
    """Anatomical neuron class: the group a bilateral or radial set belongs to.

    ``AVAL`` and ``AVAR`` are both class ``AVA``; ``RMDDL``, ``RMDDR``, ``RMDVL``,
    ``RMDVR``, ``RMDL`` and ``RMDR`` are all class ``RMD``. This grouping is
    anatomical, published, and *not* derivable by stripping letters from the name,
    which is why it is read from a table rather than computed.
    """

    overlay_id: str = "classes"
    source_id: str = "wang_2024_nt_atlas"

    def provenance(self) -> Provenance:
        return _provenance(self.source_id, "Neuron class column of Supplementary File 2.")

    def apply(self, c: Connectome) -> tuple[Connectome, OverlayReport]:
        table = {r["cell_id"]: r["class_name"] for r in _rows("neuron_classes.csv")}
        updated: dict[str, Cell] = {}
        for cell in c.cells:
            klass = table.get(cell.id)
            if klass is None:
                continue
            updated[cell.id] = replace(
                cell,
                class_name=klass,
                field_sources={**cell.field_sources, "class_name": self.source_id},
            )
        return _rebuild(c, updated), _report(self.overlay_id, set(updated), c, set(table))


@dataclass(frozen=True)
class SIMRoleOverlay:
    """Coarse sensory / interneuron / motor classification.

    Not a measurement made by any connectome. It is a curated summary of decades
    of ablation, imaging and anatomy work. A cell may carry several roles, and
    seven cells (CANL/R, MI, MCL/R, NSML/R) carry
    :attr:`~common.data.schemas.SIMRole.UNKNOWN` because the published label for
    them does not resolve to a single role and we decline to invent one.
    """

    overlay_id: str = "sim"
    source_id: str = "wormatlas_cook_2019_via_cect"

    def provenance(self) -> Provenance:
        return _provenance(self.source_id, "Coarse functional role per cell.")

    def apply(self, c: Connectome) -> tuple[Connectome, OverlayReport]:
        table: dict[str, list[SIMRole]] = {}
        for r in _rows("sim_roles.csv"):
            table.setdefault(r["cell_id"], []).append(SIMRole(r["role"]))
        updated: dict[str, Cell] = {}
        for cell in c.cells:
            roles = table.get(cell.id)
            if not roles:
                continue
            updated[cell.id] = replace(
                cell,
                roles=tuple(roles),
                field_sources={**cell.field_sources, "roles": self.source_id},
            )
        return _rebuild(c, updated), _report(self.overlay_id, set(updated), c, set(table))


@dataclass(frozen=True)
class NeurotransmitterOverlay:
    """Neurotransmitter usage from CRISPR reporter-allele expression.

    The source distinguishes three things that a naive import would flatten, and
    which are preserved in ``worm/data/annotations/neurotransmitters.csv``:

    * clean reporter expression;
    * *dim and variable* expression — real but weaker evidence, marked ``*`` in
      the source;
    * **uptake rather than synthesis** — the cell does not make the transmitter
      but takes it up from neighbours (AVFL/R for GABA; CANL/R and ASIL/R for
      betaine). Such a cell can still release the transmitter, but the genetics
      and the manipulability differ.

    Sixteen neurons express no known transmitter pathway gene at all. That is a
    published *result*, not a missing measurement, and it is recorded as
    :attr:`~common.data.schemas.Neurotransmitter.UNKNOWN` with the evidence string
    ``orphan_no_pathway_gene_detected`` rather than left blank.
    """

    overlay_id: str = "nt"
    source_id: str = "wang_2024_nt_atlas"

    def provenance(self) -> Provenance:
        return _provenance(
            self.source_id,
            "Reporter-allele neurotransmitter assignments, Supplementary File 2 "
            "(hermaphrodite).",
        )

    def apply(self, c: Connectome) -> tuple[Connectome, OverlayReport]:
        table: dict[str, list[Neurotransmitter]] = {}
        for r in _rows("neurotransmitters.csv"):
            table.setdefault(r["cell_id"], []).append(Neurotransmitter(r["neurotransmitter"]))
        updated: dict[str, Cell] = {}
        for cell in c.cells:
            nts = table.get(cell.id)
            if not nts:
                continue
            updated[cell.id] = replace(
                cell,
                neurotransmitters=tuple(nts),
                field_sources={**cell.field_sources, "neurotransmitters": self.source_id},
            )
        return _rebuild(c, updated), _report(self.overlay_id, set(updated), c, set(table))


def neurotransmitter_evidence() -> Iterator[Mapping[str, str]]:
    """Full neurotransmitter table including evidence strings and source rows."""
    yield from _rows("neurotransmitters.csv")


OVERLAYS: dict[str, AnnotationOverlay] = {
    "classes": NeuronClassOverlay(),
    "sim": SIMRoleOverlay(),
    "nt": NeurotransmitterOverlay(),
}

DEFAULT_OVERLAYS = ("classes", "sim", "nt")


def get_overlays(names: tuple[str, ...] = DEFAULT_OVERLAYS) -> list[AnnotationOverlay]:
    unknown = [n for n in names if n not in OVERLAYS]
    if unknown:
        raise KeyError(f"unknown overlay(s) {unknown}; known: {sorted(OVERLAYS)}")
    return [OVERLAYS[n] for n in names]
