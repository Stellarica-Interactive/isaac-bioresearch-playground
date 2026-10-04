"""Which connections damp the forward gait mode, and which drive the instability?

    .venv/Scripts/python tools/gait_mode_sensitivity.py --mode forward
    .venv/Scripts/python tools/gait_mode_sensitivity.py --mode unstable
    .venv/Scripts/python tools/gait_mode_sensitivity.py --check "DB>BWM-D"

At its equilibrium the committed loop has two modes that matter (model
assumptions 5AN.8): an *unstable* slow standing mode, 0.065 Hz, which grows into
the 22-second coiling cycle the loop actually runs, and a damped *forward*
travelling mode, 0.17 Hz at zeta 0.94, which is the gait. The loop needs the
second to go unstable instead of the first. A global knob scales hundreds of
pathways at once, some of which may push one way and some the other, so a
global sweep can only report their net effect.

This asks the connectome pathway by pathway instead. A *pathway* is every
chemical synapse from one cell class onto another (``DB>VD``), or every gap
junction between two classes (``AVB~VB``). Body-wall muscles are grouped by side,
``BWM-D`` and ``BWM-V``. For each pathway it computes how the forward mode's
eigenvalue moves when that pathway's conductance is scaled, from first-order
eigenvalue perturbation theory:

    d mu / d eps = w^T (dJ/d eps) v / (w^T v)

with ``v`` and ``w`` the mode's right and left eigenvectors and ``eps`` the
fractional change in the pathway's conductance, so the result is a derivative
with respect to log conductance and pathways of different size compare fairly.

``dJ/d eps`` is the *total* derivative: scaling a pathway also moves the
equilibrium, and that shift can matter as much as the direct effect. The
equilibrium's shift comes from the fixed-point
condition, ``dx0/d eps = (I - J)^-1 dM/d eps``, and ``J v`` is differentiated
along the path ``(x0(eps), eps)``. Ten evaluations of the one-step map per
pathway, against two thousand to re-linearise.

**Thresholds are held fixed.** Each cell's sigmoid midpoint is solved once from
the network at build time; a ``--param`` sweep re-solves it, a pathway scaling
here does not. That is the cleaner question -- what does this pathway's strength
do -- but it means the numbers are not the same derivative as a ``--param``
sweep. ``--check`` re-linearises the loop with one pathway actually scaled,
thresholds still fixed, and compares the prediction with what happens.

First-order: valid for small changes. The linear extrapolation to zero damping
it prints is a ranking aid, not a prediction that the loop would crawl there.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass, replace

import numpy as np

sys.path.insert(0, ".")
from common.data.schemas import CellCategory  # noqa: E402
from tools.analyse_loop import Loop, fixed_point  # noqa: E402
from tools.track_gait_mode import FORWARD_STEP_DEG, classify  # noqa: E402
from worm.loader import load  # noqa: E402
from worm.neural.config import RUNTIME_OVERLAYS  # noqa: E402


@dataclass(frozen=True)
class Pathway:
    """Every synapse of one kind between two groups of cells."""

    label: str
    kind: str  # "chemical" or "gap"
    rows: np.ndarray
    cols: np.ndarray
    synapses: int
    sign: str


@dataclass(frozen=True)
class Mode:
    mu: complex
    v: np.ndarray
    w: np.ndarray
    lam: complex
    hz: float
    zeta: float
    step: float


def cell_groups() -> dict[str, str]:
    """Cell id -> group label: the cell class, with body-wall muscles by side."""
    connectome, _ = load("cook_2019_herm", annotations=RUNTIME_OVERLAYS)
    groups = {}
    for cell in connectome.cells:
        if cell.category is CellCategory.MUSCLE:
            if cell.id.startswith("MD"):
                groups[cell.id] = "BWM-D"
            elif cell.id.startswith("MV"):
                groups[cell.id] = "BWM-V"
            else:
                groups[cell.id] = cell.id.rstrip("0123456789LR") or cell.id
        else:
            groups[cell.id] = cell.class_name or cell.id
    return groups


def pathways(loop: Loop, groups: dict[str, str], e_exc_mv: float) -> list[Pathway]:
    net = loop.runtime.network
    label = [groups.get(c, c) for c in net.cell_ids]

    chem: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for i, j in zip(*np.nonzero(net.g_syn), strict=True):
        chem[f"{label[j]}>{label[i]}"].append((i, j))
    gap: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for i, j in zip(*np.nonzero(net.g_gap), strict=True):
        a, b = sorted((label[i], label[j]))
        gap[f"{a}~{b}"].append((i, j))  # both orientations: the matrix is symmetric

    out = []
    for name, pairs in chem.items():
        rows, cols = (np.array(x) for x in zip(*pairs, strict=True))
        exc = np.isclose(net.e_rev[rows, cols], e_exc_mv)
        sign = "exc" if exc.all() else "inh" if not exc.any() else "mixed"
        out.append(Pathway(name, "chemical", rows, cols, len(pairs), sign))
    for name, pairs in gap.items():
        rows, cols = (np.array(x) for x in zip(*pairs, strict=True))
        out.append(Pathway(name, "gap", rows, cols, len(pairs), "gap"))
    return out


def scaled(base, pathway: Pathway, factor: float):
    if pathway.kind == "chemical":
        g = base.g_syn.copy()
        g[pathway.rows, pathway.cols] *= factor
        return replace(base, g_syn=g)
    g = base.g_gap.copy()
    g[pathway.rows, pathway.cols] *= factor
    return replace(base, g_gap=g)


def pick_mode(loop: Loop, jac: np.ndarray, kind: str, band: float) -> Mode:
    """The least-damped forward-travelling mode, or the fastest-growing oscillation."""
    mu_all, vec = np.linalg.eig(jac)
    # The posture history is a shift register: some multipliers are exactly 0.
    with np.errstate(divide="ignore", invalid="ignore"):
        lam_all = np.log(mu_all.astype(complex)) / loop.dt
    found = classify(loop, lam_all, vec)
    # Fractions of a hertz; see track_gait_mode for the fast modes excluded.
    found = [f for f in found if f[0] <= 5.0]
    if kind == "forward":
        found = [f for f in found if abs(f[2] - FORWARD_STEP_DEG) < band]
    if not found:
        raise SystemExit(f"no {kind} mode at this equilibrium")
    hz, zeta, step, _, k = min(found, key=lambda f: f[1])
    mu = complex(mu_all[k])
    left_mu, left = np.linalg.eig(jac.T)
    w = left[:, int(np.argmin(np.abs(left_mu - mu)))]
    v = vec[:, k] / np.abs(vec[:, k]).max()
    return Mode(mu, v, w, complex(lam_all[k]), hz, zeta, step)


def jv(loop: Loop, x: np.ndarray, v: np.ndarray, h: float) -> np.ndarray:
    """``J(x) v`` for complex ``v``, by central differences on each part."""
    out = np.zeros(x.size, dtype=complex)
    for part, unit in ((v.real, 1.0), (v.imag, 1.0j)):
        out += unit * (loop.map(x + h * part) - loop.map(x - h * part)) / (2.0 * h)
    return out


def sensitivity(loop, x0, solve, mode: Mode, perturb, delta: float, h: float):
    """d lambda / d eps for one perturbation, operating-point shift included."""
    perturb(+delta)
    m_plus = loop.map(x0)
    perturb(-delta)
    m_minus = loop.map(x0)
    dx = solve((m_plus - m_minus) / (2.0 * delta))

    perturb(+delta)
    jv_plus = jv(loop, x0 + delta * dx, mode.v, h)
    perturb(-delta)
    jv_minus = jv(loop, x0 - delta * dx, mode.v, h)
    perturb(0.0)

    dmu = (mode.w @ ((jv_plus - jv_minus) / (2.0 * delta))) / (mode.w @ mode.v)
    return dmu / (mode.mu * loop.dt)


def d_zeta(lam: complex, dlam: complex) -> float:
    a, b = lam.real, lam.imag
    r = abs(lam)
    return float(b * (-b * dlam.real + a * dlam.imag) / r**3)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--mode", choices=("forward", "unstable"), default="forward")
    parser.add_argument("--settle-s", type=float, default=10.0)
    parser.add_argument("--band", type=float, default=12.0)
    parser.add_argument("--delta", type=float, default=1.0e-2)
    parser.add_argument("--top", type=int, default=15)
    parser.add_argument("--min-synapses", type=int, default=1)
    parser.add_argument("--proprioceptive-mv", type=float, default=20.0)
    parser.add_argument("--proprio-rate", type=float, default=0.0)
    parser.add_argument(
        "--start",
        metavar="NPY",
        help="an equilibrium to continue from, where the trajectory's average is a poor guess",
    )
    parser.add_argument(
        "--check",
        action="append",
        default=[],
        metavar="PATHWAY",
        help="re-linearise with this pathway scaled by +-10%% and compare with the "
        "first-order prediction",
    )
    args = parser.parse_args()

    ns = argparse.Namespace(
        unknown_sign="exclude",
        physics_hz=240.0,
        neural_dt_ms=1.0,
        torque_scale=3.0e-3,
        proprioceptive_mv=args.proprioceptive_mv,
        proprio_rate=args.proprio_rate,
        command_mv=20.0,
        settle_s=args.settle_s,
        open_loop=False,
        param=[],
    )
    loop = Loop(ns)
    start = np.load(args.start) if args.start else None
    x0, jac = fixed_point(loop, settle_s=args.settle_s, start=start)
    residual = float(np.abs(loop.map(x0) - x0).max())
    bend = float(np.degrees(np.abs(loop.joints(x0)).max()))
    mode = pick_mode(loop, jac, args.mode, args.band)
    print(f"  equilibrium: bend {bend:.1f} deg, residual {residual:.2e}")
    print(
        f"  {args.mode} mode: {mode.hz:.3f} Hz, damping {mode.zeta:+.4f}, "
        f"phase step {mode.step:+.2f} deg/joint"
    )

    i_minus_j = np.eye(jac.shape[0]) - jac
    lu = np.linalg.inv(i_minus_j)

    def solve(b: np.ndarray) -> np.ndarray:
        return lu @ b

    base = loop.runtime.network
    e_exc = loop.runtime.model.params.e_exc_mv
    groups = cell_groups()
    paths = [p for p in pathways(loop, groups, e_exc) if p.synapses >= args.min_synapses]

    def network_perturb(pathway: Pathway):
        def apply(eps: float) -> None:
            loop.runtime.network = base if eps == 0.0 else scaled(base, pathway, 1.0 + eps)

        return apply

    proprio_base = loop.proprio

    def proprio_perturb(eps: float) -> None:
        loop.proprio = replace(
            proprio_base, gain_pa_per_rad=proprio_base.gain_pa_per_rad * (1 + eps)
        )

    command_base = dict(loop.command)

    def command_perturb(eps: float) -> None:
        loop.command = {c: v * (1.0 + eps) for c, v in command_base.items()}

    h = 1.0e-5
    rows = []
    every = Pathway("ALL chemical", "chemical", *np.nonzero(base.g_syn), 0, "-")
    all_gap = Pathway("ALL gap", "gap", *np.nonzero(base.g_gap), 0, "-")
    for p in [every, all_gap, *paths]:
        dlam = sensitivity(loop, x0, solve, mode, network_perturb(p), args.delta, h)
        rows.append((p, dlam))
    extras = []
    for name, fn in (("proprioceptive gain", proprio_perturb), ("command drive", command_perturb)):
        extras.append((name, sensitivity(loop, x0, solve, mode, fn, args.delta, h)))

    def line(label, kind, n, sign, dlam):
        dz = d_zeta(mode.lam, dlam)
        # The pathway scale at which damping would extrapolate to zero, if any.
        factor = 1.0 - mode.zeta / dz if dz != 0.0 else -1.0
        reach = f"{factor:7.2f}x" if factor > 0.0 else f"{'-':>8}"
        return (
            f"  {label:<22} {kind:>8} {n:>5} {sign:>5} {dz:+10.4f} "
            f"{dlam.real:+10.4f} {dlam.imag / (2 * np.pi):+9.4f} {reach}"
        )

    head = (
        f"  {'pathway':<22} {'kind':>8} {'syn':>5} {'sign':>5} {'dzeta':>10} "
        f"{'dgrowth/s':>10} {'dfreq Hz':>9} {'to zero':>8}"
    )
    print()
    print("  derivatives with respect to log conductance (a fractional change)")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for p, dlam in rows[:2]:
        print(line(p.label, p.kind, "", "", dlam))
    for name, dlam in extras:
        print(line(name, "input", "", "", dlam))
    total = sum((d for _, d in rows[2:]), 0j)
    print(line("sum over pathways", "", "", "", total))

    # No single pathway reaching zero does not rule out many changing together.
    # The smallest coordinated change, in log conductance, that extrapolates to
    # zero damping is along the gradient, with length |zeta| / |grad|.
    grad = np.array([d_zeta(mode.lam, d) for _, d in rows[2:]])
    norm = float(np.linalg.norm(grad))
    step = abs(mode.zeta) / norm
    biggest = abs(mode.zeta) * float(np.abs(grad).max()) / norm**2
    print(
        f"\n  smallest coordinated change reaching zero damping, linearly extrapolated:\n"
        f"    {step:.2f} in log conductance over {grad.size} pathways, "
        f"largest single pathway x{np.exp(biggest):.2f} or x{np.exp(-biggest):.2f}"
    )

    body = sorted(rows[2:], key=lambda r: d_zeta(mode.lam, r[1]))
    print()
    print(f"  pathways whose strengthening most LOWERS this mode's damping (top {args.top})")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for p, dlam in body[: args.top]:
        print(line(p.label, p.kind, p.synapses, p.sign, dlam))
    print()
    print(f"  pathways whose strengthening most RAISES it (top {args.top})")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for p, dlam in body[::-1][: args.top]:
        print(line(p.label, p.kind, p.synapses, p.sign, dlam))

    by_label = {p.label: (p, d) for p, d in rows}
    for name in args.check:
        if name not in by_label:
            print(f"\n  --check {name!r}: no such pathway")
            continue
        p, dlam = by_label[name]
        predicted = d_zeta(mode.lam, dlam)
        print(f"\n  check {name}: predicted dzeta/d eps {predicted:+.4f}")
        for eps in (-0.1, 0.1):
            loop.runtime.network = scaled(base, p, 1.0 + eps)
            # Continue *this* equilibrium: from a trajectory average Newton can
            # land on a different one (a saddle), which is not the comparison.
            _, moved_jac = fixed_point(loop, start=x0)
            moved = pick_mode(loop, moved_jac, args.mode, args.band)
            print(
                f"    scaled {1 + eps:.1f}x: damping {moved.zeta:+.4f} "
                f"(predicted {mode.zeta + predicted * eps:+.4f}), {moved.hz:.3f} Hz"
            )
        loop.runtime.network = base
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
