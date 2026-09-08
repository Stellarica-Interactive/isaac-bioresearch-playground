"""Annotation overlays.

An electron-microscopy connectome records contacts. It does not record what a
cell *does* (sensory? motor?), what transmitter it releases, or whether a synapse
excites or inhibits. Each of those comes from a different body of work by a
different group, on a different set of animals, with a different confidence.

So they are kept out of the connectome files entirely and merged in at load time.
Three consequences, all deliberate:

* a correction to a neurotransmitter atlas never dirties a connectome's committed
  golden files, and a parse-verification test cannot accidentally also be
  asserting that the annotation table is intact;
* the same annotation table applies to every dataset — AVAL is an interneuron in
  Witvliet and in Cook — so the tables cannot drift between importers;
* cells an overlay does not cover are **reported**, not silently defaulted.
  "Which of these 181 neurons has no neurotransmitter assignment" is a real
  scientific fact about the state of the field, and it should be visible.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import Protocol, runtime_checkable

from common.data.schemas import Connectome, Provenance


@dataclass(frozen=True, slots=True)
class OverlayReport:
    """What an overlay actually managed to annotate."""

    overlay_id: str
    matched: int
    """Number of cells or edges that received at least one value."""

    eligible: int
    """Number that could in principle have received one."""

    unmatched: tuple[str, ...]
    """Identifiers left un-annotated. Not an error — usually a gap in the data."""

    unknown_in_source: tuple[str, ...]
    """Identifiers present in the annotation table but absent from this dataset.

    Expected and harmless when annotating a head-only dataset with a
    whole-animal table; a red flag if it is large in the other direction.
    """

    conflicts: tuple[str, ...] = ()

    @property
    def coverage(self) -> float:
        return self.matched / self.eligible if self.eligible else 0.0

    def summary(self) -> str:
        return (
            f"{self.overlay_id}: {self.matched}/{self.eligible} "
            f"({self.coverage:.1%}) annotated, {len(self.unmatched)} unmatched, "
            f"{len(self.unknown_in_source)} table entries not in this dataset"
        )


@runtime_checkable
class AnnotationOverlay(Protocol):
    """Attaches one kind of published annotation to a connectome."""

    @property
    def overlay_id(self) -> str:
        """Short name used to select this overlay, e.g. ``"nt"``."""
        ...

    def provenance(self) -> Provenance: ...

    def apply(self, c: Connectome) -> tuple[Connectome, OverlayReport]: ...


def apply_overlays(
    c: Connectome, overlays: Sequence[AnnotationOverlay]
) -> tuple[Connectome, list[OverlayReport]]:
    """Apply overlays in order, threading provenance through."""
    reports: list[OverlayReport] = []
    for ov in overlays:
        prov = ov.provenance()
        sources = dict(c.sources)
        sources[prov.source_id] = prov
        c = replace(c, sources=sources)
        c, report = ov.apply(c)
        reports.append(report)
    return c, reports
