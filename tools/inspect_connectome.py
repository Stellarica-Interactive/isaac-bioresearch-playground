"""Inspect a normalized connectome from the command line.

    python tools/inspect_connectome.py --organism worm list
    python tools/inspect_connectome.py --dataset witvliet_2021_7 summary
    python tools/inspect_connectome.py --dataset witvliet_2021_7 cell ASHL
    python tools/inspect_connectome.py --dataset witvliet_2021_7 verify
    python tools/inspect_connectome.py --dataset witvliet_2021_7 gaps

Two commands are deliberately separate because they fail for different reasons:

``validate``
    Is the graph internally well-formed? (No dangling endpoints, no duplicate
    edges, gap junctions stored canonically, every cited source resolvable.)

``verify``
    Do our totals agree with the figures published for this dataset by someone
    else? A parser can be perfectly self-consistent and still wrong.

Exit codes: 0 success, 1 assertion failure, 2 usage error.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from typing import Any

from common.data.registry import list_datasets
from common.data.schemas import Connectome, SIMRole, SynapseType
from common.data.validate import validate
from common.neural import graph as gr
from worm.annotations.overlays import DEFAULT_OVERLAYS, neurotransmitter_evidence
from worm.importers.sources import (
    compare_to_reference,
    has_reference,
    normalized_path,
    raw_path,
    reference_totals,
)
from worm.loader import load

EXIT_OK, EXIT_FAIL, EXIT_USAGE = 0, 1, 2


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def table(rows: Sequence[Sequence[Any]], headers: Sequence[str]) -> str:
    cells = [[str(x) for x in r] for r in rows]
    widths = [
        max(len(h), *(len(r[i]) for r in cells)) if cells else len(h)
        for i, h in enumerate(headers)
    ]
    line = "  ".join(h.ljust(w) for h, w in zip(headers, widths, strict=True))
    rule = "  ".join("-" * w for w in widths)
    body = "\n".join(
        "  ".join(c.ljust(w) for c, w in zip(r, widths, strict=True)) for r in cells
    )
    return f"{line}\n{rule}\n{body}" if body else f"{line}\n{rule}\n(none)"


def emit(args: argparse.Namespace, payload: dict[str, Any], text: str) -> None:
    if args.format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True, default=str))
    else:
        print(text)


def provenance_stamp(c: Connectome) -> dict[str, Any]:
    """Included in every JSON payload so a number carries its own citation."""
    return {
        "dataset_id": c.id,
        "schema_version": c.schema_version,
        "source_sha256": c.provenance.sha256,
        "citation": c.provenance.citation,
        "scope": str(c.scope),
    }


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


def cmd_list(args: argparse.Namespace) -> int:
    rows = []
    for did in list_datasets():
        raw = raw_path(did)
        norm = normalized_path(did)
        rows.append(
            [
                did,
                "yes" if (norm / "meta.json").exists() else "-",
                "yes" if raw.exists() else "-",
                "yes" if has_reference(did) else "-",
            ]
        )
    payload = {
        "datasets": [
            {"id": r[0], "normalized": r[1] == "yes", "raw": r[2] == "yes", "oracle": r[3] == "yes"}
            for r in rows
        ]
    }
    emit(args, payload, table(rows, ["dataset", "normalized", "raw", "published totals"]))
    return EXIT_OK


def cmd_summary(args: argparse.Namespace) -> int:
    c, reports = _load(args)
    t = c.totals()
    stats = gr.summary_stats(c)

    role_counts = {
        str(r): len(c.cells_with_role(r)) for r in SIMRole if c.cells_with_role(r)
    }
    nt_counts: dict[str, int] = {}
    for cell in c.cells:
        for nt in cell.neurotransmitters:
            nt_counts[str(nt)] = nt_counts.get(str(nt), 0) + 1

    lines = [
        f"{c.id}  ({c.organism}, {c.sex}, {c.stage})",
        f"  scope:     {c.scope}",
        f"  source:    {c.provenance.citation}",
        f"  licence:   {c.provenance.license}",
    ]
    if c.provenance.notes:
        lines.append(f"  caveat:    {c.provenance.notes}")
    lines += [
        "",
        "CELLS",
        f"  neurons {t.neurons}   muscles {t.muscles}   glia {t.glia}   "
        f"other {t.other}   total {t.cells}",
        "",
        "CONNECTIONS (as measured: counts of presynaptic active zones)",
        f"  chemical    {t.chemical_edges:5d} directed edges,  weight sum {t.chemical_weight}",
        f"  electrical  {t.electrical_edges_undirected:5d} gap junctions stored once, "
        f"weight sum {t.electrical_weight_undirected}",
        f"              {t.electrical_edges_directed:5d} entries / {t.electrical_weight_directed} "
        f"weight in the symmetrized form published tables use",
        f"              ({t.electrical_self_loops} of them join a cell to itself)",
        "",
        "GRAPH",
    ] + [f"  {k:32s} {v:.4g}" for k, v in stats.items()]

    if role_counts:
        lines += ["", "FUNCTIONAL ROLE (annotation, not measured by this dataset)"]
        lines += [f"  {k:12s} {v}" for k, v in sorted(role_counts.items())]
    if nt_counts:
        lines += ["", "NEUROTRANSMITTER (annotation, not measured by this dataset)"]
        lines += [f"  {k:16s} {v}" for k, v in sorted(nt_counts.items(), key=lambda x: -x[1])]
    if reports:
        lines += ["", "ANNOTATION COVERAGE"] + [f"  {r.summary()}" for r in reports]

    payload = {
        **provenance_stamp(c),
        "totals": t.to_dict(),
        "graph": stats,
        "roles": role_counts,
        "neurotransmitters": nt_counts,
        "coverage": [
            {"overlay": r.overlay_id, "matched": r.matched, "eligible": r.eligible}
            for r in reports
        ],
    }
    emit(args, payload, "\n".join(lines))
    return EXIT_OK


def cmd_cell(args: argparse.Namespace) -> int:
    c, _ = _load(args)
    try:
        cell = c.cell(args.cell_id)
    except KeyError as exc:
        print(exc, file=sys.stderr)
        return EXIT_FAIL

    nt_detail = {
        r["cell_id"]: r for r in neurotransmitter_evidence() if r["cell_id"] == cell.id
    }
    p = c.partners(cell.id)

    def rows(edges: list[Any], other: str) -> list[list[Any]]:
        return [
            [getattr(e, other), e.weight, str(e.synapse_type)]
            for e in edges
        ]

    lines = [f"{cell.id}   class {cell.class_name or '?'}   {cell.category}"]
    lines.append(f"  roles:             {', '.join(map(str, cell.roles)) or '(not annotated)'}")
    if cell.neurotransmitters:
        nt = ", ".join(map(str, cell.neurotransmitters))
        detail = nt_detail.get(cell.id)
        extra = (
            f"   [{detail['evidence']}, source label {detail['source_label']!r}]"
            if detail
            else ""
        )
        lines.append(f"  neurotransmitter:  {nt}{extra}")
    lines.append(f"  labels in source:  {', '.join(cell.source_labels)}")
    lines.append("  field sources:")
    for k, v in sorted(cell.field_sources.items()):
        lines.append(f"    {k:20s} {v}")

    lines += ["", f"OUTGOING CHEMICAL ({len(p['chemical_out'])})",
              table(rows(p["chemical_out"], "post"), ["post", "synapses", "type"])]
    lines += ["", f"INCOMING CHEMICAL ({len(p['chemical_in'])})",
              table(rows(p["chemical_in"], "pre"), ["pre", "synapses", "type"])]
    gap = [[e.post if e.pre == cell.id else e.pre, e.weight] for e in p["electrical"]]
    lines += ["", f"GAP JUNCTIONS ({len(gap)})", table(gap, ["partner", "synapses"])]

    payload = {
        **provenance_stamp(c),
        "cell": {
            "id": cell.id,
            "class_name": cell.class_name,
            "category": str(cell.category),
            "roles": [str(r) for r in cell.roles],
            "neurotransmitters": [str(n) for n in cell.neurotransmitters],
            "field_sources": dict(cell.field_sources),
        },
        "chemical_out": [{"post": e.post, "weight": e.weight} for e in p["chemical_out"]],
        "chemical_in": [{"pre": e.pre, "weight": e.weight} for e in p["chemical_in"]],
        "electrical": [
            {"partner": e.post if e.pre == cell.id else e.pre, "weight": e.weight}
            for e in p["electrical"]
        ],
    }
    emit(args, payload, "\n".join(lines))
    return EXIT_OK


def cmd_edges(args: argparse.Namespace) -> int:
    c, _ = _load(args)
    stype = SynapseType(args.type) if args.type else None
    edges = c.edges(stype, pre=args.pre, post=args.post, min_weight=args.min_weight)
    edges.sort(key=lambda e: (-e.weight, e.pre, e.post))
    if args.limit:
        edges = edges[: args.limit]
    rows = [[e.pre, e.post, str(e.synapse_type), e.weight] for e in edges]
    payload = {
        **provenance_stamp(c),
        "edges": [
            {"pre": e.pre, "post": e.post, "type": str(e.synapse_type), "weight": e.weight}
            for e in edges
        ],
    }
    emit(args, payload, table(rows, ["pre", "post", "type", "synapses"]))
    return EXIT_OK


def cmd_role(args: argparse.Namespace) -> int:
    c, _ = _load(args)
    cells = c.cells_with_role(SIMRole(args.role))
    rows = [
        [x.id, x.class_name or "", ", ".join(map(str, x.neurotransmitters))] for x in cells
    ]
    payload = {**provenance_stamp(c), "role": args.role, "cells": [x.id for x in cells]}
    emit(args, payload, table(rows, ["cell", "class", "neurotransmitter"]))
    return EXIT_OK


def cmd_validate(args: argparse.Namespace) -> int:
    c, _ = _load(args)
    problems = validate(c)
    payload = {**provenance_stamp(c), "violations": [str(p) for p in problems]}
    text = (
        f"{c.id}: OK, no structural violations"
        if not problems
        else f"{c.id}: {len(problems)} violation(s)\n" + "\n".join(f"  {p}" for p in problems)
    )
    emit(args, payload, text)
    return EXIT_OK if not problems else EXIT_FAIL


def cmd_verify(args: argparse.Namespace) -> int:
    c, _ = _load(args)
    observed = c.totals().to_dict()
    if not has_reference(c.id):
        emit(
            args,
            {**provenance_stamp(c), "verified": False, "reason": "no published totals"},
            f"{c.id}: UNVERIFIED - no entry in worm/data/reference_totals.toml",
        )
        return EXIT_FAIL

    ref = reference_totals()[c.id]
    mismatches = compare_to_reference(c.id, observed)
    checked = [(k, v) for k, v in ref.items() if isinstance(v, int)]
    rows = [
        [k, v, observed.get(k, "-"), "ok" if observed.get(k) == v else "MISMATCH"]
        for k, v in sorted(checked)
    ]
    text = (
        f"{c.id}  oracle: {ref.get('source', '?')}\n\n"
        + table(rows, ["field", "published", "parsed", ""])
        + (
            f"\n\n{len(checked)} field(s) checked, all agree."
            if not mismatches
            else f"\n\n{len(mismatches)} MISMATCH(ES)."
        )
    )
    payload = {
        **provenance_stamp(c),
        "verified": not mismatches,
        "oracle": ref.get("source"),
        "fields": {k: {"published": v, "parsed": observed.get(k)} for k, v in checked},
    }
    emit(args, payload, text)
    return EXIT_OK if not mismatches else EXIT_FAIL


def cmd_gaps(args: argparse.Namespace) -> int:
    """Everything in the loaded graph that was not measured by this reconstruction."""
    c, _ = _load(args)
    rows_by_source: dict[str, dict[str, int]] = {}
    for kind, _elem, fieldname, src in c.unmeasured():
        rows_by_source.setdefault(src, {}).setdefault(f"{kind}.{fieldname}", 0)
        rows_by_source[src][f"{kind}.{fieldname}"] += 1

    lines = [
        f"{c.id}: values present in the loaded graph that this dataset did NOT measure.",
        "",
        f"MEASURED by {c.provenance.source_id}:",
        "  cell presence, connection presence, synapse type, synapse count",
        "",
        "NOT MEASURED - supplied by other sources:",
    ]
    rows = []
    for src, fields in sorted(rows_by_source.items()):
        for fieldname, n in sorted(fields.items()):
            rows.append([fieldname, n, src, c.sources[src].citation[:70] + "..."])
    lines.append(table(rows, ["field", "count", "source_id", "citation"]))
    lines += [
        "",
        "NOT AVAILABLE AT ALL - no source in this repository provides these:",
        "  - synaptic sign (excitatory / inhibitory): not measurable by electron",
        "    microscopy; only predictable from receptor gene expression.",
        "  - synaptic strength / conductance: the weight above is a count of",
        "    physical active zones, not a physiological quantity.",
        "  - the animal's neural activity, or anything it had learned.",
        "  - extrasynaptic signalling (neuropeptides, monoamines acting at a",
        "    distance), which leaves no ultrastructural trace.",
        "",
        "See docs/model_assumptions.md.",
    ]
    payload = {**provenance_stamp(c), "unmeasured_by_source": rows_by_source}
    emit(args, payload, "\n".join(lines))
    return EXIT_OK


def cmd_graph(args: argparse.Namespace) -> int:
    c, _ = _load(args)

    if args.metric == "components":
        comps = gr.components(c)
        comp_rows: list[Sequence[Any]] = [
            [i, len(s), ", ".join(s[:8]) + ("..." if len(s) > 8 else "")]
            for i, s in enumerate(comps)
        ]
        emit(
            args,
            {**provenance_stamp(c), "components": comps},
            table(comp_rows, ["#", "size", "members"]),
        )
        return EXIT_OK

    if args.metric == "reciprocity":
        r = gr.reciprocity(c)
        emit(
            args,
            {**provenance_stamp(c), "chemical_reciprocity": r},
            f"chemical reciprocity: {r:.4f}",
        )
        return EXIT_OK

    if args.metric == "betweenness":
        ranked = sorted(gr.betweenness(c).items(), key=lambda kv: -kv[1])[: args.top]
        emit(
            args,
            {**provenance_stamp(c), "betweenness": dict(ranked)},
            table([[k, f"{v:.5f}"] for k, v in ranked], ["cell", "betweenness"]),
        )
        return EXIT_OK

    key = {"degree": "total_degree", "in-degree": "in_degree", "out-degree": "out_degree"}[
        args.metric
    ]
    top = gr.top_by(gr.degrees(c), key, args.top)
    emit(
        args,
        {
            **provenance_stamp(c),
            "degrees": [
                {
                    "cell_id": d.cell_id,
                    "in_degree": d.in_degree,
                    "out_degree": d.out_degree,
                    "in_weight": d.in_weight,
                    "out_weight": d.out_weight,
                    "gap_partners": d.gap_partners,
                    "gap_weight": d.gap_weight,
                }
                for d in top
            ],
        },
        table(
            [
                [d.cell_id, d.in_degree, d.out_degree, d.in_weight, d.out_weight, d.gap_partners]
                for d in top
            ],
            ["cell", "chem in", "chem out", "in wt", "out wt", "gap partners"],
        ),
    )
    return EXIT_OK


def cmd_compare(args: argparse.Namespace) -> int:
    a, _ = load(args.a, annotations=(), prefer=args.prefer)
    b, _ = load(args.b, annotations=(), prefer=args.prefer)
    a_cells, b_cells = set(a.cell_ids()), set(b.cell_ids())
    a_edges = {(e.pre, e.post, str(e.synapse_type)): e.weight for e in a.connections}
    b_edges = {(e.pre, e.post, str(e.synapse_type)): e.weight for e in b.connections}
    shared = set(a_edges) & set(b_edges)

    lines = [
        f"{a.id}  vs  {b.id}",
        "",
        f"  cells:  {len(a_cells)} vs {len(b_cells)};  shared {len(a_cells & b_cells)}",
        f"    only in {a.id}: {sorted(a_cells - b_cells)}",
        f"    only in {b.id}: {sorted(b_cells - a_cells)}",
        "",
        f"  connections: {len(a_edges)} vs {len(b_edges)};  shared {len(shared)}",
        f"    only in {a.id}: {len(set(a_edges) - set(b_edges))}",
        f"    only in {b.id}: {len(set(b_edges) - set(a_edges))}",
        "",
        "  NOTE: these are two different animals reconstructed by the same method.",
        "  A connection present in one and not the other may reflect genuine",
        "  individual variability, developmental stage, or reconstruction",
        "  uncertainty. This comparison cannot distinguish between those.",
    ]
    payload = {
        "a": a.id,
        "b": b.id,
        "cells_only_a": sorted(a_cells - b_cells),
        "cells_only_b": sorted(b_cells - a_cells),
        "edges_a": len(a_edges),
        "edges_b": len(b_edges),
        "edges_shared": len(shared),
    }
    emit(args, payload, "\n".join(lines))
    return EXIT_OK


def cmd_coverage(args: argparse.Namespace) -> int:
    c, reports = _load(args)
    lines = [f"{c.id} annotation coverage", ""]
    for r in reports:
        lines.append(f"  {r.summary()}")
        if r.unmatched:
            lines.append(f"    unmatched: {', '.join(r.unmatched)}")
    payload = {
        **provenance_stamp(c),
        "coverage": [
            {
                "overlay": r.overlay_id,
                "matched": r.matched,
                "eligible": r.eligible,
                "unmatched": list(r.unmatched),
            }
            for r in reports
        ],
    }
    emit(args, payload, "\n".join(lines))
    return EXIT_OK


# ---------------------------------------------------------------------------


def _load(args: argparse.Namespace) -> tuple[Connectome, list[Any]]:
    ann = () if args.annotations == "none" else tuple(args.annotations.split(","))
    return load(args.dataset, annotations=ann, prefer=args.prefer)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="inspect_connectome",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--organism", default="worm", choices=["worm"])
    p.add_argument("--dataset", default="witvliet_2021_7")
    p.add_argument("--format", default="table", choices=["table", "json"])
    p.add_argument(
        "--annotations",
        default=",".join(DEFAULT_OVERLAYS),
        help="comma-separated overlay names, or 'none'. Synaptic polarity is never "
        "in the default set: it is predicted, not measured.",
    )
    p.add_argument("--prefer", default="normalized", choices=["normalized", "raw"])

    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list").set_defaults(func=cmd_list)
    sub.add_parser("summary").set_defaults(func=cmd_summary)
    sub.add_parser("validate").set_defaults(func=cmd_validate)
    sub.add_parser("verify").set_defaults(func=cmd_verify)
    sub.add_parser("gaps").set_defaults(func=cmd_gaps)
    sub.add_parser("coverage").set_defaults(func=cmd_coverage)

    c = sub.add_parser("cell")
    c.add_argument("cell_id")
    c.set_defaults(func=cmd_cell)

    e = sub.add_parser("edges")
    e.add_argument("--pre")
    e.add_argument("--post")
    e.add_argument("--type", choices=["chemical", "electrical"])
    e.add_argument("--min-weight", type=int, default=1)
    e.add_argument("--limit", type=int, default=50)
    e.set_defaults(func=cmd_edges)

    r = sub.add_parser("role")
    r.add_argument("role", choices=[str(x) for x in SIMRole])
    r.set_defaults(func=cmd_role)

    g = sub.add_parser("graph")
    g.add_argument(
        "--metric",
        default="degree",
        choices=["degree", "in-degree", "out-degree", "betweenness", "components", "reciprocity"],
    )
    g.add_argument("--top", type=int, default=20)
    g.set_defaults(func=cmd_graph)

    cm = sub.add_parser("compare")
    cm.add_argument("a")
    cm.add_argument("b")
    cm.set_defaults(func=cmd_compare)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args))
    except (KeyError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_FAIL


if __name__ == "__main__":
    sys.exit(main())
