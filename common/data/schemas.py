"""Normalized, organism-agnostic connectome schema.

This module defines *our* internal representation of a nervous system graph. It is
deliberately independent of any published dataset's file format, and of any
neuroscience framework.

Design principles
-----------------

1. **Anatomy and annotation are separate.** A :class:`Connectome` produced by an
   importer contains only what the source electron-microscopy reconstruction
   actually measured: which cells exist and which of them contact which, how
   often, and by what kind of synapse. Everything else — whether a cell is
   sensory or motor, which neurotransmitter it releases, whether a chemical
   synapse excites or inhibits — comes from *different* experiments by
   *different* groups and is attached later by an annotation overlay
   (:mod:`common.data.overlay`).

2. **Every field can name its source.** Each :class:`Cell` and :class:`Connection`
   carries ``field_sources``, mapping a field name to a ``source_id`` that
   resolves in :attr:`Connectome.sources`. Combined with
   :class:`Confidence`, this makes the project's central honesty requirement
   machine-checkable: you can ask a loaded connectome to list everything that is
   not directly measured.

3. **Gap junctions are stored once.** Electrical synapses are physically
   bidirectional. Storing them twice (once per direction) is a common convention
   in published adjacency matrices but makes weight sums and degree counts
   ambiguous. We store each gap junction exactly once with ``pre <= post``, and
   reproduce the published directed form on demand via
   :meth:`Connectome.electrical_directed_view`.

Biological note for the reader
------------------------------

A "connectome" here is a *wiring diagram reconstructed from serial-section
electron micrographs of one fixed, dead animal*. It records physical contacts.
It does not record synaptic strength, sign, neurotransmitter identity, the
animal's ongoing neural activity, or anything learned during its life. Those are
separate measurements — several of which do not exist for most cells.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field, replace
from enum import StrEnum
from types import MappingProxyType
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    import networkx as nx
    import numpy as np

SCHEMA_VERSION = "1.0"

# Canonical cell identifiers: uppercase letters/digits plus '_' (for cells such
# as ``exc_gl``) and lowercase for non-neuronal cells such as ``pm1`` or ``hyp``.
CELL_ID_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


# ---------------------------------------------------------------------------
# Enumerations
# ---------------------------------------------------------------------------


class CellCategory(StrEnum):
    """What kind of cell a graph node represents.

    Connectome datasets are not neuron-only. Witvliet #7 has 222 nodes of which
    only 181 are neurons; Cook 2019 has 473 nodes of which only 302 are neurons.
    The rest are muscles (the actual output of the nervous system), glia and
    other cells such as the excretory gland.
    """

    NEURON = "neuron"
    MUSCLE = "muscle"
    GLIA = "glia"
    OTHER = "other"


class SynapseType(StrEnum):
    """How two cells are coupled.

    ``CHEMICAL``
        A directed synapse. The presynaptic cell releases a neurotransmitter
        which binds receptors on the postsynaptic cell. Whether the effect is
        excitatory or inhibitory depends on the *postsynaptic receptor*, which
        electron microscopy cannot see — see :class:`Sign`.

    ``ELECTRICAL``
        A gap junction: a direct cytoplasmic channel between two cells, passing
        current in both directions. Undirected, and stored once in this schema.
    """

    CHEMICAL = "chemical"
    ELECTRICAL = "electrical"

    @property
    def is_directed(self) -> bool:
        return self is SynapseType.CHEMICAL


class SIMRole(StrEnum):
    """Coarse functional role: Sensory / Interneuron / Motor.

    This is a *published annotation*, not a measurement made by the connectome
    reconstruction itself. It comes from decades of ablation, imaging and
    anatomy work, summarized in sources such as WormAtlas and Cook et al. 2019.

    A cell may carry more than one role — many *C. elegans* neurons are
    polymodal (for example ASH is a sensory neuron that also makes substantial
    interneuron-like connections), which is why :attr:`Cell.roles` is a tuple.
    """

    SENSORY = "sensory"
    INTER = "inter"
    MOTOR = "motor"
    MODULATORY = "modulatory"
    UNKNOWN = "unknown"


class Neurotransmitter(StrEnum):
    """Signalling molecule a cell releases.

    ``NONE`` means "positively determined to release none of the assayed
    transmitters"; ``UNKNOWN`` means "not assayed / not reported". Keeping these
    distinct matters: the first is a result, the second is a gap.

    Co-transmission (a cell releasing two or more of these) is common, which is
    why :attr:`Cell.neurotransmitters` is a tuple.
    """

    ACETYLCHOLINE = "acetylcholine"
    GLUTAMATE = "glutamate"
    GABA = "gaba"
    DOPAMINE = "dopamine"
    SEROTONIN = "serotonin"
    OCTOPAMINE = "octopamine"
    TYRAMINE = "tyramine"
    BETAINE = "betaine"
    NEUROPEPTIDE = "neuropeptide"
    UNKNOWN = "unknown"
    NONE = "none"


class Sign(StrEnum):
    """Whether a chemical synapse depolarizes or hyperpolarizes its target.

    **This is never measured by a connectome.** Electron microscopy shows that a
    synapse exists and how large it is; it cannot show which receptor the
    postsynaptic cell expresses, and that receptor is what determines the sign.

    Values here therefore always arrive from a separate source — typically a
    gene-expression-based prediction (Fenyves et al. 2020) — and always carry
    :attr:`Connection.sign_confidence` of :attr:`Confidence.PREDICTED`.
    """

    EXCITATORY = "excitatory"
    INHIBITORY = "inhibitory"
    MIXED = "mixed"
    UNKNOWN = "unknown"


class Confidence(StrEnum):
    """How a value came to be, ordered from strongest to weakest evidence.

    ``MEASURED``
        Directly observed in the cited experiment. In practice this covers the
        connectome itself: cell presence, synapse presence, synapse counts.

    ``PUBLISHED_ANNOTATION``
        A curated classification from the literature, resting on many separate
        experiments (for example "AVAL is an interneuron", or a neurotransmitter
        assignment from a reporter-strain atlas).

    ``PREDICTED``
        Computationally inferred from other data, not directly observed. Synapse
        polarity predicted from receptor gene expression lives here.

    ``ASSUMED``
        Our own modelling choice or engineering simplification, with no direct
        biological backing. Anything tagged this way must also appear in
        ``docs/model_assumptions.md``.
    """

    MEASURED = "measured"
    PUBLISHED_ANNOTATION = "published_annotation"
    PREDICTED = "predicted"
    ASSUMED = "assumed"


class Scope(StrEnum):
    """Anatomical extent of a reconstruction.

    This is load-bearing, not decorative. Witvliet et al. 2021 reconstructed the
    *brain* (nerve ring and immediately adjacent ganglia), so those datasets
    contain no ventral-nerve-cord motor neurons and only the anterior eight body
    wall muscle quadrant segments. A model of crawling locomotion cannot be built
    from a ``HEAD`` dataset, and code that silently tries will produce a worm
    that cannot move.
    """

    HEAD = "head"
    WHOLE_ANIMAL = "whole_animal"
    PHARYNX = "pharynx"


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Provenance:
    """Where one body of data came from.

    One :class:`Provenance` describes a single source *file* (a published
    supplementary spreadsheet, an annotation table we derived from a paper),
    not a whole project.
    """

    source_id: str
    """Stable key referenced by ``field_sources`` (e.g. ``"witvliet_2021_7"``)."""

    kind: str
    """``"connectome"`` or ``"annotation"``."""

    citation: str
    """Full human-readable citation, suitable for pasting into a paper."""

    url: str
    """Exact location the file was retrieved from."""

    license: str
    """Licence under which the source is distributed."""

    doi: str | None = None
    file_name: str | None = None
    sha256: str | None = None
    retrieved: str | None = None
    """ISO-8601 date the file was downloaded."""

    notes: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "source_id": self.source_id,
            "kind": self.kind,
            "citation": self.citation,
            "url": self.url,
            "license": self.license,
            "doi": self.doi,
            "file_name": self.file_name,
            "sha256": self.sha256,
            "retrieved": self.retrieved,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, d: Mapping[str, object]) -> Provenance:
        return cls(
            source_id=str(d["source_id"]),
            kind=str(d["kind"]),
            citation=str(d["citation"]),
            url=str(d["url"]),
            license=str(d["license"]),
            doi=_opt_str(d.get("doi")),
            file_name=_opt_str(d.get("file_name")),
            sha256=_opt_str(d.get("sha256")),
            retrieved=_opt_str(d.get("retrieved")),
            notes=str(d.get("notes") or ""),
        )


def _opt_str(v: object) -> str | None:
    if v is None or v == "":
        return None
    return str(v)


# ---------------------------------------------------------------------------
# Graph elements
# ---------------------------------------------------------------------------


def _freeze(m: Mapping[str, str] | None) -> Mapping[str, str]:
    return MappingProxyType(dict(m or {}))


@dataclass(frozen=True, slots=True)
class Cell:
    """One node of the connectome graph.

    Named ``Cell`` rather than ``Neuron`` on purpose: a substantial minority of
    nodes in every real dataset are muscles, glia or other cells, and code that
    assumes otherwise silently drops the nervous system's actual output.
    """

    id: str
    """Canonical identifier, e.g. ``"AVAL"``, ``"MVL01"``, ``"exc_gl"``."""

    category: CellCategory

    class_name: str | None = None
    """Anatomical class shared by a bilateral pair, e.g. ``"AVA"`` for AVAL/AVAR."""

    roles: tuple[SIMRole, ...] = ()
    """Empty until a role annotation overlay is applied."""

    neurotransmitters: tuple[Neurotransmitter, ...] = ()
    """Empty until a neurotransmitter overlay is applied."""

    source_labels: tuple[str, ...] = ()
    """Raw label(s) as written in the source file, e.g. ``("BWM-VL01",)``.

    Retained so that a canonicalization bug is diagnosable after the fact.
    """

    field_sources: Mapping[str, str] = field(default_factory=dict)
    """Field name -> ``source_id``. Keys resolve in :attr:`Connectome.sources`."""

    def __post_init__(self) -> None:
        object.__setattr__(self, "field_sources", _freeze(self.field_sources))
        object.__setattr__(self, "roles", tuple(self.roles))
        object.__setattr__(self, "neurotransmitters", tuple(self.neurotransmitters))
        object.__setattr__(self, "source_labels", tuple(self.source_labels))

    @property
    def is_neuron(self) -> bool:
        return self.category is CellCategory.NEURON

    def has_role(self, role: SIMRole) -> bool:
        return role in self.roles


@dataclass(frozen=True, slots=True)
class Connection:
    """One edge of the connectome graph.

    Invariant enforced by :mod:`common.data.validate`: if
    ``synapse_type is SynapseType.ELECTRICAL`` then ``pre <= post``, because gap
    junctions are stored exactly once.
    """

    pre: str
    post: str
    synapse_type: SynapseType

    weight: int
    """Synapse count as published.

    For Witvliet this is the number of presynaptic active zones agreed by at
    least two of three independent human annotators. It is a count of physical
    structures, **not** a physiological conductance — converting one to the other
    is a modelling assumption we have not yet made.
    """

    sign: Sign = Sign.UNKNOWN
    sign_confidence: Confidence = Confidence.ASSUMED
    field_sources: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "field_sources", _freeze(self.field_sources))

    @property
    def key(self) -> tuple[str, str, SynapseType]:
        return (self.pre, self.post, self.synapse_type)

    def reversed(self) -> Connection:
        """Same edge with endpoints swapped (used for the directed gap-junction view)."""
        return replace(self, pre=self.post, post=self.pre)


# ---------------------------------------------------------------------------
# Totals
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Totals:
    """Summary counts, in the same conventions published papers use.

    The ``electrical_directed_*`` fields exist so our numbers can be compared
    against published figures without changing how we store data. See
    :meth:`Connectome.electrical_directed_view`.
    """

    cells: int
    neurons: int
    muscles: int
    glia: int
    other: int
    chemical_edges: int
    chemical_weight: int
    electrical_edges_undirected: int
    electrical_weight_undirected: int
    electrical_edges_directed: int
    electrical_weight_directed: int
    electrical_self_loops: int

    @property
    def other_cells(self) -> int:
        """Cells that are neither neurons nor muscles.

        Published factsheets report this as a single "other cells" figure; we
        keep glia and everything else apart internally because the distinction is
        real, and expose the sum for comparison against those figures.
        """
        return self.glia + self.other

    def to_dict(self) -> dict[str, int]:
        return {
            "cells": self.cells,
            "neurons": self.neurons,
            "muscles": self.muscles,
            "glia": self.glia,
            "other": self.other,
            "other_cells": self.other_cells,
            "chemical_edges": self.chemical_edges,
            "chemical_weight": self.chemical_weight,
            "electrical_edges_undirected": self.electrical_edges_undirected,
            "electrical_weight_undirected": self.electrical_weight_undirected,
            "electrical_edges_directed": self.electrical_edges_directed,
            "electrical_weight_directed": self.electrical_weight_directed,
            "electrical_self_loops": self.electrical_self_loops,
        }


# ---------------------------------------------------------------------------
# Connectome
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Connectome:
    """A complete normalized nervous-system graph plus its provenance."""

    id: str
    organism: str
    sex: str
    stage: str
    scope: Scope
    cells: tuple[Cell, ...]
    connections: tuple[Connection, ...]
    provenance: Provenance
    sources: Mapping[str, Provenance] = field(default_factory=dict)
    schema_version: str = SCHEMA_VERSION

    def __post_init__(self) -> None:
        object.__setattr__(self, "cells", tuple(self.cells))
        object.__setattr__(self, "connections", tuple(self.connections))
        sources = dict(self.sources)
        sources.setdefault(self.provenance.source_id, self.provenance)
        object.__setattr__(self, "sources", MappingProxyType(sources))

    # -- lookup ------------------------------------------------------------

    def _index(self) -> dict[str, Cell]:
        # Rebuilt on demand; frozen dataclass with a tuple field cannot cache
        # cheaply without __dict__, and these graphs are small.
        return {c.id: c for c in self.cells}

    def cell(self, cell_id: str) -> Cell:
        try:
            return self._index()[cell_id]
        except KeyError:
            raise KeyError(f"no such cell in {self.id!r}: {cell_id!r}") from None

    def has_cell(self, cell_id: str) -> bool:
        return cell_id in self._index()

    def cell_ids(self) -> list[str]:
        return sorted(c.id for c in self.cells)

    def cells_by_category(self, category: CellCategory) -> list[Cell]:
        return sorted((c for c in self.cells if c.category is category), key=lambda c: c.id)

    def neurons(self) -> list[Cell]:
        return self.cells_by_category(CellCategory.NEURON)

    def cells_with_role(self, role: SIMRole) -> list[Cell]:
        return sorted((c for c in self.cells if role in c.roles), key=lambda c: c.id)

    # -- edges -------------------------------------------------------------

    def edges(
        self,
        synapse_type: SynapseType | None = None,
        *,
        pre: str | None = None,
        post: str | None = None,
        min_weight: int = 1,
    ) -> list[Connection]:
        out = []
        for e in self.connections:
            if synapse_type is not None and e.synapse_type is not synapse_type:
                continue
            if pre is not None and e.pre != pre:
                continue
            if post is not None and e.post != post:
                continue
            if e.weight < min_weight:
                continue
            out.append(e)
        return out

    def chemical(self) -> list[Connection]:
        return self.edges(SynapseType.CHEMICAL)

    def electrical(self) -> list[Connection]:
        """Gap junctions, stored form: each pair exactly once, ``pre <= post``."""
        return self.edges(SynapseType.ELECTRICAL)

    def electrical_directed_view(self) -> Iterator[Connection]:
        """Gap junctions in the symmetrized directed form used by published tables.

        Every non-self gap junction is emitted twice (``a->b`` and ``b->a``) with
        identical weight; a self-loop is emitted once. This is the convention
        behind figures such as "576 electrical entries, weight sum 794" for
        Witvliet #7, and lets us assert against those numbers without adopting a
        storage format that double-counts.
        """
        for e in self.connections:
            if e.synapse_type is not SynapseType.ELECTRICAL:
                continue
            yield e
            if e.pre != e.post:
                yield e.reversed()

    def partners(self, cell_id: str) -> dict[str, list[Connection]]:
        """Incoming/outgoing chemical edges and gap junctions touching ``cell_id``."""
        out: dict[str, list[Connection]] = {
            "chemical_out": [],
            "chemical_in": [],
            "electrical": [],
        }
        for e in self.connections:
            if e.synapse_type is SynapseType.CHEMICAL:
                if e.pre == cell_id:
                    out["chemical_out"].append(e)
                elif e.post == cell_id:
                    out["chemical_in"].append(e)
            elif cell_id in (e.pre, e.post):
                out["electrical"].append(e)
        for v in out.values():
            v.sort(key=lambda e: (-e.weight, e.pre, e.post))
        return out

    # -- matrices and graphs -----------------------------------------------

    def adjacency(
        self,
        synapse_type: SynapseType,
        order: Sequence[str] | None = None,
    ) -> tuple[list[str], np.ndarray]:
        """Weighted adjacency matrix.

        Node order is explicit and sorted by default — never insertion order, so
        that repeated runs and separate processes agree. Electrical matrices are
        symmetrized on the way out.
        """
        import numpy as np

        nodes = list(order) if order is not None else self.cell_ids()
        idx = {n: i for i, n in enumerate(nodes)}
        m = np.zeros((len(nodes), len(nodes)), dtype=np.int64)
        edges = (
            self.electrical_directed_view()
            if synapse_type is SynapseType.ELECTRICAL
            else iter(self.edges(synapse_type))
        )
        for e in edges:
            i, j = idx.get(e.pre), idx.get(e.post)
            if i is None or j is None:
                continue
            m[i, j] += e.weight
        return nodes, m

    def to_networkx(self, synapse_type: SynapseType | None = None) -> nx.MultiDiGraph:
        """NetworkX view. Gap junctions appear in both directions."""
        import networkx as nx

        g: nx.MultiDiGraph = nx.MultiDiGraph(name=self.id)
        for c in self.cells:
            g.add_node(
                c.id,
                category=str(c.category),
                class_name=c.class_name,
                roles=[str(r) for r in c.roles],
                neurotransmitters=[str(n) for n in c.neurotransmitters],
            )
        if synapse_type is not SynapseType.ELECTRICAL:
            for e in self.chemical():
                g.add_edge(
                    e.pre, e.post, key="chemical", weight=e.weight, sign=str(e.sign),
                    synapse_type="chemical",
                )
        if synapse_type is not SynapseType.CHEMICAL:
            for e in self.electrical_directed_view():
                g.add_edge(
                    e.pre, e.post, key="electrical", weight=e.weight, synapse_type="electrical"
                )
        return g

    def subgraph(self, seeds: Iterable[str], hops: int = 1) -> Connectome:
        """Neighbourhood around ``seeds``, following edges in either direction.

        Used to produce readable circuit diagrams instead of hairballs.
        """
        keep = {s for s in seeds if self.has_cell(s)}
        missing = set(seeds) - keep
        if missing:
            raise KeyError(f"seeds not in {self.id!r}: {sorted(missing)}")
        frontier = set(keep)
        for _ in range(hops):
            nxt: set[str] = set()
            for e in self.connections:
                if e.pre in frontier:
                    nxt.add(e.post)
                if e.post in frontier:
                    nxt.add(e.pre)
            nxt -= keep
            keep |= nxt
            frontier = nxt
            if not frontier:
                break
        return replace(
            self,
            id=f"{self.id}:subgraph",
            cells=tuple(c for c in self.cells if c.id in keep),
            connections=tuple(e for e in self.connections if e.pre in keep and e.post in keep),
        )

    # -- summary -----------------------------------------------------------

    def totals(self) -> Totals:
        by_cat: dict[CellCategory, int] = dict.fromkeys(CellCategory, 0)
        for c in self.cells:
            by_cat[c.category] += 1
        chem = self.chemical()
        elec = self.electrical()
        self_loops = sum(1 for e in elec if e.pre == e.post)
        return Totals(
            cells=len(self.cells),
            neurons=by_cat[CellCategory.NEURON],
            muscles=by_cat[CellCategory.MUSCLE],
            glia=by_cat[CellCategory.GLIA],
            other=by_cat[CellCategory.OTHER],
            chemical_edges=len(chem),
            chemical_weight=sum(e.weight for e in chem),
            electrical_edges_undirected=len(elec),
            electrical_weight_undirected=sum(e.weight for e in elec),
            electrical_edges_directed=2 * (len(elec) - self_loops) + self_loops,
            electrical_weight_directed=sum(e.weight for e in self.electrical_directed_view()),
            electrical_self_loops=self_loops,
        )

    def unmeasured(self) -> list[tuple[str, str, str, str]]:
        """Every annotated field whose evidence is weaker than ``MEASURED``.

        Returns ``(kind, element_id, field_name, source_id)`` rows. This backs the
        ``inspect_connectome gaps`` command: the point of the project is that this
        list is easy to produce and impossible to lose track of.
        """
        rows: list[tuple[str, str, str, str]] = []
        for c in self.cells:
            for fname, src in sorted(c.field_sources.items()):
                if src != self.provenance.source_id:
                    rows.append(("cell", c.id, fname, src))
        for e in self.connections:
            for fname, src in sorted(e.field_sources.items()):
                if src != self.provenance.source_id:
                    rows.append(("connection", f"{e.pre}->{e.post}", fname, src))
        return rows
