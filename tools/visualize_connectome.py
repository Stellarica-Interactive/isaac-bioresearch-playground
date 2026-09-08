"""Draw a readable subnetwork of a connectome.

    python tools/visualize_connectome.py --circuit gentle-touch --hops 1 -o out.png
    python tools/visualize_connectome.py --cells ASHL,ASHR --hops 2 -o out.png
    python tools/visualize_connectome.py --list-circuits

Drawing all 222 cells at once produces a hairball that conveys nothing, so this
tool always works on a neighbourhood: seed cells plus everything within ``--hops``
of them. Named circuits live in ``worm/data/circuits/circuits.json`` as data, with
a citation each, so the biology stays reviewable and the tool stays generic.

Encoding:

* node colour  -- functional role (sensory / interneuron / motor)
* node shape   -- cell category (neuron = circle, muscle = square, other = diamond)
* solid arrow  -- chemical synapse, pointing from pre- to postsynaptic
* dashed line  -- gap junction, undirected
* line width   -- synapse count

Note what is **not** encoded: whether a chemical synapse excites or inhibits. That
is not in the data, so it is not in the picture.

Requires the ``viz`` extra: ``pip install -e ".[viz]"``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from importlib.resources import files
from pathlib import Path
from typing import TYPE_CHECKING, Any

from common.data.schemas import Cell, CellCategory, Connectome, SIMRole
from worm.loader import load

if TYPE_CHECKING:  # pragma: no cover - typing only
    import networkx as nx

ROLE_COLOUR = {
    SIMRole.SENSORY: "#4C9F70",
    SIMRole.INTER: "#3E6DA8",
    SIMRole.MOTOR: "#C4553B",
    SIMRole.MODULATORY: "#8E6BB0",
    SIMRole.UNKNOWN: "#9A9A9A",
}
CATEGORY_SHAPE = {
    CellCategory.NEURON: "o",
    CellCategory.MUSCLE: "s",
    CellCategory.GLIA: "^",
    CellCategory.OTHER: "D",
}
NO_ROLE_COLOUR = "#D8D8D8"


def circuits() -> dict[str, dict]:
    text = (files("worm.data") / "circuits" / "circuits.json").read_text(encoding="utf-8")
    return json.loads(text)["circuits"]


def _colour(cell: Cell) -> str:
    for role in (SIMRole.SENSORY, SIMRole.MOTOR, SIMRole.INTER, SIMRole.MODULATORY):
        if role in cell.roles:
            return ROLE_COLOUR[role]
    return NO_ROLE_COLOUR


def _layout(g: nx.DiGraph, nx: Any) -> dict[str, tuple[float, float]]:
    """Position nodes, then push apart any pair that would overlap.

    Force-directed layouts routinely place two nodes close enough that their
    labels collide, which is the single thing most likely to make a circuit
    diagram useless. A plain spring layout plus this relaxation pass keeps the
    tool free of a scipy dependency (which ``kamada_kawai_layout`` would need).
    """
    n = max(len(g), 1)
    pos = nx.spring_layout(g, seed=0, k=2.4 / math.sqrt(n), iterations=500)

    min_sep = max(0.18, 1.25 / math.sqrt(n))
    nodes = list(pos)
    for _ in range(80):
        moved = False
        for i, a in enumerate(nodes):
            for b in nodes[i + 1 :]:
                ax_, ay = pos[a]
                bx, by = pos[b]
                dx, dy = bx - ax_, by - ay
                d = math.hypot(dx, dy)
                if d >= min_sep:
                    continue
                if d < 1e-9:
                    dx, dy, d = 1e-3, 0.0, 1e-3
                push = (min_sep - d) / 2.0
                ux, uy = dx / d, dy / d
                pos[a] = (ax_ - ux * push, ay - uy * push)
                pos[b] = (bx + ux * push, by + uy * push)
                moved = True
        if not moved:
            break
    return pos


def draw(sub: Connectome, out: Path, title: str, *, seeds: set[str]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import networkx as nx

    g = nx.DiGraph()
    g.add_nodes_from(sub.cell_ids())
    for e in sub.chemical():
        g.add_edge(e.pre, e.post, weight=e.weight)

    pos = _layout(g, nx)

    # Capped: beyond about 90 cells the picture stops being readable anyway, and
    # an uncapped figure size produces multi-megabyte PNGs.
    size = min(16.0, max(9.0, 0.34 * len(g)))
    fig, ax = plt.subplots(figsize=(size, size * 0.78))

    for category, shape in CATEGORY_SHAPE.items():
        nodes = [c.id for c in sub.cells if c.category is category]
        if not nodes:
            continue
        nx.draw_networkx_nodes(
            g,
            pos,
            nodelist=nodes,
            node_shape=shape,
            node_color=[_colour(sub.cell(n)) for n in nodes],
            node_size=[900 if n in seeds else 480 for n in nodes],
            edgecolors=["#111111" if n in seeds else "#666666" for n in nodes],
            linewidths=[2.0 if n in seeds else 0.7 for n in nodes],
            ax=ax,
        )

    chem = sub.chemical()
    if chem:
        wmax = max(e.weight for e in chem)
        nx.draw_networkx_edges(
            g,
            pos,
            edgelist=[(e.pre, e.post) for e in chem],
            width=[0.4 + 2.6 * e.weight / wmax for e in chem],
            edge_color="#555555",
            alpha=0.55,
            arrows=True,
            arrowsize=9,
            connectionstyle="arc3,rad=0.08",
            node_size=520,
            ax=ax,
        )

    elec = [e for e in sub.electrical() if e.pre != e.post]
    if elec:
        wmax = max(e.weight for e in elec)
        gj = nx.Graph()
        gj.add_nodes_from(g.nodes)
        for e in elec:
            gj.add_edge(e.pre, e.post)
        nx.draw_networkx_edges(
            gj,
            pos,
            edgelist=[(e.pre, e.post) for e in elec],
            width=[0.8 + 3.0 * e.weight / wmax for e in elec],
            edge_color="#C9A227",
            style="dashed",
            alpha=0.85,
            arrows=False,
            ax=ax,
        )

    nx.draw_networkx_labels(g, pos, font_size=7.0, font_family="DejaVu Sans", ax=ax)

    legend = [
        mpatches.Patch(color=ROLE_COLOUR[SIMRole.SENSORY], label="sensory"),
        mpatches.Patch(color=ROLE_COLOUR[SIMRole.INTER], label="interneuron"),
        mpatches.Patch(color=ROLE_COLOUR[SIMRole.MOTOR], label="motor"),
        mpatches.Patch(color=NO_ROLE_COLOUR, label="not a neuron / unannotated"),
    ]
    ax.legend(handles=legend, loc="lower left", fontsize=8, frameon=False)

    ax.set_title(title, fontsize=11)
    ax.text(
        0.5,
        -0.02,
        "solid arrow = chemical synapse (pre -> post)   dashed = gap junction   "
        "width = synapse count\nSynaptic sign is not shown because it is not measured.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8,
        color="#444444",
    )
    ax.axis("off")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--organism", default="worm", choices=["worm"])
    p.add_argument("--dataset", default="witvliet_2021_7")
    p.add_argument("--circuit", help="named circuit from worm/data/circuits/circuits.json")
    p.add_argument("--cells", help="comma-separated seed cell ids")
    p.add_argument("--hops", type=int, default=1)
    p.add_argument("--max-cells", type=int, default=90, help="refuse to draw a hairball")
    p.add_argument("-o", "--out", type=Path, default=Path("docs/img/circuit.png"))
    p.add_argument("--list-circuits", action="store_true")
    args = p.parse_args(argv)

    known = circuits()
    if args.list_circuits:
        for name, spec in sorted(known.items()):
            print(f"{name}\n  {spec['description']}\n  seeds: {', '.join(spec['seeds'])}")
            print(f"  source: {spec['source']}")
            if spec.get("note"):
                print(f"  note: {spec['note']}")
            print()
        return 0

    if bool(args.circuit) == bool(args.cells):
        p.error("give exactly one of --circuit or --cells (or --list-circuits)")

    if args.circuit:
        if args.circuit not in known:
            print(f"unknown circuit {args.circuit!r}; known: {sorted(known)}", file=sys.stderr)
            return 2
        spec = known[args.circuit]
        seeds = list(spec["seeds"])
        title = f"{args.circuit} - {args.dataset}"
    else:
        seeds = [s.strip() for s in args.cells.split(",") if s.strip()]
        title = f"{', '.join(seeds)} +{args.hops} hop(s) - {args.dataset}"

    c, _ = load(args.dataset)
    present = [s for s in seeds if c.has_cell(s)]
    absent = [s for s in seeds if not c.has_cell(s)]
    if absent:
        print(f"not present in {args.dataset} (skipped): {absent}", file=sys.stderr)
    if not present:
        print("no seed cells present in this dataset", file=sys.stderr)
        return 1

    sub = c.subgraph(present, hops=args.hops)
    if len(sub.cells) > args.max_cells:
        print(
            f"{len(sub.cells)} cells at {args.hops} hop(s) exceeds --max-cells "
            f"({args.max_cells}). Reduce --hops or raise the limit; a denser graph than "
            "this is not readable.",
            file=sys.stderr,
        )
        return 1

    try:
        import matplotlib  # noqa: F401
    except ImportError:
        print('matplotlib not installed. Run: pip install -e ".[viz]"', file=sys.stderr)
        return 1
    draw(sub, args.out, title, seeds=set(present))

    t = sub.totals()
    print(
        f"wrote {args.out}  ({t.cells} cells, {t.chemical_edges} chemical, "
        f"{t.electrical_edges_undirected} gap junctions)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
