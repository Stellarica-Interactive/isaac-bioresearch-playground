"""Draw a readable view of a connectome.

    python tools/visualize_connectome.py --circuit gentle-touch --hops 1 -o out.png
    python tools/visualize_connectome.py --cells ASHL,ASHR --hops 2 -o out.png
    python tools/visualize_connectome.py --whole -o out.png
    python tools/visualize_connectome.py --list-circuits

Two modes, because "show me the connectome" and "show me this circuit" want
completely different pictures.

**Circuit mode** (``--circuit`` / ``--cells``). Drawing all 222 cells with a force
layout produces a hairball, so this mode works on a neighbourhood: seed cells plus
everything within ``--hops`` of them. Named circuits live in
``worm/data/circuits/circuits.json`` as data, with a citation each, so the biology
stays reviewable and the tool stays generic.

**Whole-network mode** (``--whole``). Concentric rings by functional role -- motor
innermost, then interneurons, sensory, and body wall muscles outermost -- so the
sensory-to-muscle organisation is visible at a glance. Radius encodes role;
*angular* position is taken from a force layout of the real graph, so neighbours on
a ring really are wiring neighbours rather than alphabetical accidents.

Encoding, both modes:

* node colour  -- functional role (sensory / interneuron / motor / muscle)
* solid line   -- chemical synapse (circuit mode draws the arrowhead)
* gold line    -- gap junction, undirected
* line width   -- synapse count
* node size    -- number of partners

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


def _has_matplotlib() -> bool:
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return False
    return True


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


# --- whole-network ring view ------------------------------------------------
# Ring 0 sensory, 1 interneuron, 2 motor, 3 muscle, 4 glia/other. A cell that is
# both sensory and motor lands on the sensory ring; the ambiguity is real and is
# reported in docs/neurons.md rather than resolved here.
RING_OF_ROLE = {SIMRole.SENSORY: 0, SIMRole.MOTOR: 2}
RING_RADIUS = {0: 1.00, 1: 0.64, 2: 0.30, 3: 1.36, 4: 1.36}
RING_LABEL = {0: "SENSORY NEURONS", 1: "INTERNEURONS", 2: "MOTOR NEURONS", 3: "BODY WALL MUSCLES"}
RING_COLOUR = {
    0: "#3FBF7F",
    1: "#4C8FD6",
    2: "#E4643C",
    3: "#A879D0",
    4: "#9AA3AE",
}


def _ring_of(cell: Cell) -> int:
    if cell.category is CellCategory.MUSCLE:
        return 3
    if cell.category is not CellCategory.NEURON:
        return 4
    for role, ring in RING_OF_ROLE.items():
        if role in cell.roles:
            return ring
    return 1


def _ring_positions(
    c: Connectome, nx: Any
) -> tuple[dict[str, tuple[float, float]], dict[int, list[str]]]:
    """Radius from functional role, angle from a force layout of the real graph.

    Placing cells alphabetically around each ring would scatter every circuit; taking
    the angle from a force layout keeps wiring neighbours adjacent. Spacing each ring
    evenly afterwards is what stops the whole thing collapsing onto one side, which a
    naive "point each node at the mean angle of its partners" relaxation does.
    """
    g = nx.Graph()
    g.add_nodes_from(c.cell_ids())
    for e in c.connections:
        if e.pre == e.post:
            continue
        prev = g.get_edge_data(e.pre, e.post, {}).get("weight", 0)
        g.add_edge(e.pre, e.post, weight=prev + e.weight)
    force = nx.spring_layout(g, seed=11, iterations=600, weight="weight")
    angle0 = {n: math.atan2(force[n][1], force[n][0]) for n in g}

    members: dict[int, list[str]] = {r: [] for r in RING_RADIUS}
    for cell in c.cells:
        members[_ring_of(cell)].append(cell.id)

    pos: dict[str, tuple[float, float]] = {}
    for ring, ids in members.items():
        ids.sort(key=lambda n: angle0[n])
        for k, n in enumerate(ids):
            # Inner rings get a small offset so bilateral pairs do not line up
            # radially with the ring outside them and hide each other's labels.
            a = 2 * math.pi * k / max(len(ids), 1) + (0.35 if ring in (1, 2) else 0.0)
            pos[n] = (RING_RADIUS[ring] * math.cos(a), RING_RADIUS[ring] * math.sin(a))
    return pos, members


def _curve(
    p: tuple[float, float], q: tuple[float, float], bend: float = 0.22
) -> tuple[list[float], list[float]]:
    """Bow an edge toward the centre so long chords do not slice across the figure."""
    mx, my = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
    cx, cy = mx * (1 - bend), my * (1 - bend)
    ts = [i / 18 for i in range(19)]
    xs = [(1 - t) ** 2 * p[0] + 2 * (1 - t) * t * cx + t * t * q[0] for t in ts]
    ys = [(1 - t) ** 2 * p[1] + 2 * (1 - t) * t * cy + t * t * q[1] for t in ts]
    return xs, ys


def draw_whole(c: Connectome, out: Path, *, label_top: int = 34) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import networkx as nx

    pos, members = _ring_positions(c, nx)

    partners: dict[str, set[str]] = {i: set() for i in c.cell_ids()}
    for e in c.connections:
        partners[e.pre].add(e.post)
        partners[e.post].add(e.pre)
    deg = {n: len(p) for n, p in partners.items()}

    fig, ax = plt.subplots(figsize=(13, 13), facecolor="#0b0c10")
    ax.set_facecolor("#0b0c10")
    for r in (3, 0, 1, 2):
        ax.add_patch(plt.Circle((0, 0), RING_RADIUS[r], fill=False, color="#20232d", lw=1.0))

    chem = c.chemical()
    wmax = max((e.weight for e in chem), default=1)
    for e in chem:
        if e.pre == e.post:
            continue
        w = e.weight / wmax
        xs, ys = _curve(pos[e.pre], pos[e.post])
        ax.plot(xs, ys, color="#7396c4", lw=0.12 + 2.0 * w**0.55,
                alpha=0.09 + 0.5 * w**0.5, solid_capstyle="round", zorder=1)

    elec = [e for e in c.electrical() if e.pre != e.post]
    ewmax = max((e.weight for e in elec), default=1)
    for e in elec:
        xs, ys = _curve(pos[e.pre], pos[e.post])
        ax.plot(xs, ys, color="#F0C04A", lw=0.5 + 2.6 * (e.weight / ewmax) ** 0.55,
                alpha=0.55, solid_capstyle="round", zorder=2)

    for ring, ids in members.items():
        if not ids:
            continue
        ax.scatter([pos[n][0] for n in ids], [pos[n][1] for n in ids],
                   s=[16 + 8.0 * deg[n] for n in ids], c=RING_COLOUR[ring],
                   edgecolors="#0b0c10", linewidths=0.8, zorder=3)

    for n in sorted(c.cell_ids(), key=lambda n: -deg[n])[:label_top]:
        x, y = pos[n]
        norm = math.hypot(x, y) or 1.0
        ax.text(x + 0.055 * x / norm, y + 0.055 * y / norm, n, fontsize=7.0,
                color="#e8ebf0", ha="center", va="center", zorder=4)

    for ring, lbl in RING_LABEL.items():
        if members.get(ring):
            ax.text(0, RING_RADIUS[ring] + 0.05, lbl, color=RING_COLOUR[ring],
                    fontsize=10, ha="center", va="bottom", alpha=0.9, zorder=5)

    t = c.totals()
    ax.set_title(
        f"{c.organism} nervous system, arranged by function\n"
        f"{c.id}  ·  {t.cells} cells  ·  {t.chemical_edges} chemical synapses  ·  "
        f"{t.electrical_edges_undirected} gap junctions",
        color="#f4f5f7", fontsize=15, pad=22,
    )
    counts = {r: len(ids) for r, ids in members.items()}
    ax.legend(
        handles=[
            mpatches.Patch(color=RING_COLOUR[0], label=f"sensory neuron  ({counts.get(0, 0)})"),
            mpatches.Patch(color=RING_COLOUR[1], label=f"interneuron  ({counts.get(1, 0)})"),
            mpatches.Patch(color=RING_COLOUR[2], label=f"motor neuron  ({counts.get(2, 0)})"),
            mpatches.Patch(color=RING_COLOUR[3], label=f"body wall muscle  ({counts.get(3, 0)})"),
            mpatches.Patch(color=RING_COLOUR[4], label=f"glia / other  ({counts.get(4, 0)})"),
            mpatches.Patch(color="#7396c4", label="chemical synapse"),
            mpatches.Patch(color="#F0C04A", label="gap junction"),
        ],
        loc="lower left", fontsize=9.5, frameon=False, labelcolor="#cfd4dc",
        bbox_to_anchor=(-0.03, -0.02),
    )
    ax.text(
        0.5, -0.035,
        "Rings are functional role, not physical position.  Angular order preserves "
        "wiring neighbourhood.  Node size = number of partners.\n"
        "Synaptic sign is not drawn, because no connectome measures it.",
        transform=ax.transAxes, ha="center", va="top", fontsize=9.5, color="#7d8695",
    )
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_aspect("equal")
    ax.axis("off")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out,
        dpi=100,
        facecolor=fig.get_facecolor(),
        bbox_inches="tight",
        pil_kwargs={"optimize": True},
    )
    plt.close(fig)


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
    p.add_argument(
        "--whole",
        action="store_true",
        help="draw the entire network as concentric rings by functional role",
    )
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

    if not _has_matplotlib():
        print('matplotlib not installed. Run: pip install -e ".[viz]"', file=sys.stderr)
        return 1

    if args.whole:
        if args.circuit or args.cells:
            p.error("--whole draws the entire network; do not also give --circuit or --cells")
        c, _ = load(args.dataset)
        draw_whole(c, args.out)
        t = c.totals()
        print(
            f"wrote {args.out}  ({t.cells} cells, {t.chemical_edges} chemical, "
            f"{t.electrical_edges_undirected} gap junctions)"
        )
        return 0

    if bool(args.circuit) == bool(args.cells):
        p.error("give exactly one of --circuit, --cells or --whole (or --list-circuits)")

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

    draw(sub, args.out, title, seeds=set(present))

    t = sub.totals()
    print(
        f"wrote {args.out}  ({t.cells} cells, {t.chemical_edges} chemical, "
        f"{t.electrical_edges_undirected} gap junctions)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
