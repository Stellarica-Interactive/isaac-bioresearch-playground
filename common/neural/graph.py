"""Graph analysis over a normalized connectome.

Organism-agnostic, and deliberately thin: it wraps NetworkX rather than
reimplementing graph algorithms, and exists mainly to fix the conventions that
are easy to get wrong when a graph mixes directed chemical synapses with
undirected gap junctions.

Interpretation warning
----------------------

Every metric here is a statement about a *wiring diagram*, not about
information flow. A short path from a sensory neuron to a motor neuron does not
mean a signal travels it; a high-betweenness cell is not necessarily important.
The connectome constrains function without determining it — signalling that
leaves no ultrastructural trace (neuropeptides, monoamines acting
extrasynaptically) is entirely invisible here, and in *C. elegans* that
extrasynaptic network is known to be large. Treat these numbers as descriptions
of anatomy.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import networkx as nx

from common.data.schemas import Connectome, SynapseType


@dataclass(frozen=True, slots=True)
class DegreeRow:
    cell_id: str
    in_degree: int
    out_degree: int
    in_weight: int
    out_weight: int
    gap_partners: int
    gap_weight: int

    @property
    def total_degree(self) -> int:
        return self.in_degree + self.out_degree


def degrees(c: Connectome) -> list[DegreeRow]:
    """Per-cell chemical in/out degree and gap-junction partner count.

    Gap-junction partners are counted once per neighbouring cell, not twice.
    """
    ids = c.cell_ids()
    zeros = dict.fromkeys(ids, 0)
    in_deg, out_deg = dict(zeros), dict(zeros)
    in_w, out_w = dict(zeros), dict(zeros)
    gap_n, gap_w = dict(zeros), dict(zeros)

    for e in c.chemical():
        out_deg[e.pre] += 1
        in_deg[e.post] += 1
        out_w[e.pre] += e.weight
        in_w[e.post] += e.weight
    for e in c.electrical():
        gap_n[e.pre] += 1
        gap_w[e.pre] += e.weight
        if e.pre != e.post:
            gap_n[e.post] += 1
            gap_w[e.post] += e.weight

    return [DegreeRow(i, in_deg[i], out_deg[i], in_w[i], out_w[i], gap_n[i], gap_w[i]) for i in ids]


def _simple_digraph(c: Connectome, synapse_type: SynapseType | None) -> nx.DiGraph:
    """Collapse the multigraph to one weighted edge per ordered pair."""
    g: nx.DiGraph = nx.DiGraph()
    g.add_nodes_from(c.cell_ids())
    for u, v, data in c.to_networkx(synapse_type).edges(data=True):
        if g.has_edge(u, v):
            g[u][v]["weight"] += data["weight"]
        else:
            g.add_edge(u, v, weight=data["weight"])
    return g


def betweenness(
    c: Connectome, synapse_type: SynapseType | None = None, *, weighted: bool = False
) -> dict[str, float]:
    """Betweenness centrality.

    ``weighted=False`` treats every connection as one hop. That is usually what
    you want: edge *weight* here is a synapse count, and a high synapse count means
    a stronger connection, i.e. a **shorter** effective distance — the opposite of
    what NetworkX's ``weight`` parameter assumes. When ``weighted=True`` we invert
    the counts so the metric means what a reader expects.
    """
    g = _simple_digraph(c, synapse_type)
    if not weighted:
        return nx.betweenness_centrality(g)
    for _, _, data in g.edges(data=True):
        data["distance"] = 1.0 / data["weight"]
    return nx.betweenness_centrality(g, weight="distance")


def components(c: Connectome, synapse_type: SynapseType | None = None) -> list[list[str]]:
    """Weakly connected components, largest first."""
    g = _simple_digraph(c, synapse_type)
    return sorted((sorted(s) for s in nx.weakly_connected_components(g)), key=len, reverse=True)


def reciprocity(c: Connectome) -> float:
    """Fraction of chemical connections whose reverse also exists.

    Gap junctions are excluded: they are undirected, so reciprocity is undefined
    for them and including them would trivially inflate the figure.
    """
    pairs = {(e.pre, e.post) for e in c.chemical()}
    if not pairs:
        return 0.0
    return sum(1 for u, v in pairs if (v, u) in pairs) / len(pairs)


def shortest_paths(
    c: Connectome,
    source: str,
    target: str,
    *,
    synapse_type: SynapseType | None = None,
    max_paths: int = 10,
    cutoff: int | None = None,
) -> list[list[str]]:
    """Shortest hop-count paths from ``source`` to ``target``.

    Chemical synapses are followed in their signalling direction; gap junctions in
    both. A path here is an *anatomical* possibility, nothing more.
    """
    for cell in (source, target):
        if not c.has_cell(cell):
            raise KeyError(f"no such cell in {c.id!r}: {cell!r}")
    g = _simple_digraph(c, synapse_type)
    if not nx.has_path(g, source, target):
        return []
    out: list[list[str]] = []
    for path in nx.all_shortest_paths(g, source, target):
        if cutoff is not None and len(path) - 1 > cutoff:
            break
        out.append(path)
        if len(out) >= max_paths:
            break
    return out


def summary_stats(c: Connectome) -> dict[str, float]:
    g = _simple_digraph(c, None)
    rows = degrees(c)
    n = len(rows) or 1
    return {
        "nodes": g.number_of_nodes(),
        "edges_collapsed": g.number_of_edges(),
        "density": nx.density(g),
        "mean_chemical_in_degree": sum(r.in_degree for r in rows) / n,
        "mean_chemical_out_degree": sum(r.out_degree for r in rows) / n,
        "mean_gap_partners": sum(r.gap_partners for r in rows) / n,
        "chemical_reciprocity": reciprocity(c),
        "weakly_connected_components": nx.number_weakly_connected_components(g),
    }


def top_by(rows: Sequence[DegreeRow], key: str, n: int = 20) -> list[DegreeRow]:
    return sorted(rows, key=lambda r: getattr(r, key), reverse=True)[:n]
