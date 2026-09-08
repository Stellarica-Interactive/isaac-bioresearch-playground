"""Structural invariants for a :class:`~common.data.schemas.Connectome`.

Shared by the test suite and by ``tools/inspect_connectome.py validate`` so the
two can never disagree about what "well-formed" means.

This checks *internal consistency* only. Agreement with the totals published in
the source paper is a different question, answered by
``tools/inspect_connectome.py verify`` against ``worm/data/reference_totals.toml``.
The two fail for entirely different reasons and are deliberately kept apart.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from common.data.schemas import CELL_ID_RE, Confidence, Connectome, Sign, SynapseType


@dataclass(frozen=True, slots=True)
class Violation:
    code: str
    message: str

    def __str__(self) -> str:  # pragma: no cover - trivial
        return f"[{self.code}] {self.message}"


def validate(c: Connectome) -> list[Violation]:
    """Return every structural problem found. Empty list means the graph is sound."""
    v: list[Violation] = []
    ids = {cell.id for cell in c.cells}

    if len(ids) != len(c.cells):
        dupes = [k for k, n in Counter(cell.id for cell in c.cells).items() if n > 1]
        v.append(Violation("duplicate_cell", f"duplicate cell ids: {sorted(dupes)}"))

    for cell in c.cells:
        if not CELL_ID_RE.match(cell.id):
            v.append(Violation("bad_cell_id", f"cell id is not canonical: {cell.id!r}"))

    seen: set[tuple[str, str, SynapseType]] = set()
    for e in c.connections:
        if e.pre not in ids:
            v.append(Violation("dangling_pre", f"edge {e.pre}->{e.post} has unknown pre {e.pre!r}"))
        if e.post not in ids:
            v.append(
                Violation("dangling_post", f"edge {e.pre}->{e.post} has unknown post {e.post!r}")
            )
        if e.weight < 1:
            v.append(
                Violation(
                    "bad_weight",
                    f"edge {e.pre}->{e.post} ({e.synapse_type}) has weight {e.weight}; "
                    "a connection that is not observed should be absent, not zero-weighted",
                )
            )
        if e.synapse_type is SynapseType.ELECTRICAL and e.pre > e.post:
            v.append(
                Violation(
                    "electrical_not_canonical",
                    f"gap junction {e.pre}->{e.post} violates the pre <= post storage rule; "
                    "gap junctions are undirected and stored exactly once",
                )
            )
        if e.key in seen:
            v.append(
                Violation("duplicate_edge", f"duplicate edge {e.pre}->{e.post} {e.synapse_type}")
            )
        seen.add(e.key)

        if e.sign is not Sign.UNKNOWN and e.sign_confidence is Confidence.MEASURED:
            v.append(
                Violation(
                    "sign_claimed_measured",
                    f"edge {e.pre}->{e.post} claims a MEASURED sign. No connectome measures "
                    "synaptic polarity; it is inferred from receptor gene expression at best",
                )
            )

    known_sources = set(c.sources)
    for cell in c.cells:
        for fname, src in cell.field_sources.items():
            if src not in known_sources:
                v.append(
                    Violation(
                        "unresolved_source",
                        f"cell {cell.id} field {fname!r} cites unknown source {src!r}",
                    )
                )
    for e in c.connections:
        for fname, src in e.field_sources.items():
            if src not in known_sources:
                v.append(
                    Violation(
                        "unresolved_source",
                        f"edge {e.pre}->{e.post} field {fname!r} cites unknown source {src!r}",
                    )
                )

    return v


def assert_valid(c: Connectome) -> None:
    """Raise :class:`ValueError` listing every violation, or return silently."""
    problems = validate(c)
    if problems:
        joined = "\n  ".join(str(p) for p in problems)
        raise ValueError(f"connectome {c.id!r} failed validation:\n  {joined}")
