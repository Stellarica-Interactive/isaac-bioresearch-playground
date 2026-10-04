"""How much negative conductance, in which cells, would sustain the gait?

    .venv/Scripts/python tools/gait_mode_requirement.py
    .venv/Scripts/python tools/gait_mode_requirement.py --validate B-type BWM

model_assumptions 5AN.9 to 5AN.12: the forward gait mode is damped (zeta 0.94),
no connection weight changes that, no measured cell amplifies in the gait band,
and the mode lives mostly in the body-wall muscle cells -- the B-type motor
neurons carry 3% of it. Whatever sustains a gait has to inject energy into that
mode, and this asks how much it would take, cell group by cell group.

The probe is a conductance ``dG`` added to every cell of a group, with its
reversal at that cell's own equilibrium voltage:

    I = dG (V_eq - V)

so the equilibrium does not move at all and the Jacobian changes only through
the new term. A negative ``dG`` is a source. It is applied by raising the cell's
gap-junction row sum by ``dG`` -- which adds ``-dG V`` to its current and is
seen by the exponential integrator's conductance as well -- and injecting the
constant ``dG V_eq``.

A constant negative conductance is regenerative at every frequency, and that is
the point of reporting *what goes unstable first*: a real cell that sustained
the gait would have to be negative near the gait frequency and not at 0 Hz, or
it would latch instead. To first order a band-limited element tuned to the mode
shifts its eigenvalue by the same amount, so the threshold for the gait mode
read off here is the requirement on such an element; whether the other modes
would also go is the thing a constant conductance answers pessimistically.

``--validate`` re-linearises with the conductance actually added, at 0.5x and
1.0x the predicted threshold, and prints the gait mode's damping and the fastest
growth of any other mode -- the honest check that a first-order threshold is
real and that nothing else goes unstable first.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace

import numpy as np

sys.path.insert(0, ".")
from tools.analyse_loop import Loop, fixed_point, jacobian  # noqa: E402
from tools.gait_mode_sensitivity import cell_groups, pick_mode, sensitivity  # noqa: E402

#: Cell groups to probe: label -> the class names (or muscle sides) it covers.
GROUPS: dict[str, tuple[str, ...]] = {
    "B-type": ("VB", "DB"),
    "A-type": ("VA", "DA"),
    "D-type": ("VD", "DD"),
    "AS": ("AS",),
    "VNC motor": ("VB", "DB", "VA", "DA", "VD", "DD", "AS", "VC"),
    "AVB": ("AVB",),
    "head motor": ("RMD", "SMD", "SMB", "RIV"),
    "BWM-D": ("BWM-D",),
    "BWM-V": ("BWM-V",),
    "BWM": ("BWM-D", "BWM-V"),
}


def with_conductance(base, rows: np.ndarray, dg: float):
    """``base`` with ``dg`` nS added to each listed cell's own conductance."""
    net = replace(base)
    row_sum = base.gap_row_sum.copy()
    row_sum[rows] += dg
    object.__setattr__(net, "gap_row_sum", row_sum)
    return net


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--settle-s", type=float, default=10.0)
    parser.add_argument("--band", type=float, default=12.0)
    parser.add_argument("--delta", type=float, default=1.0e-2, help="nS, for the derivative")
    parser.add_argument("--validate", nargs="*", default=[], metavar="GROUP")
    args = parser.parse_args()

    ns = argparse.Namespace(
        unknown_sign="exclude",
        physics_hz=240.0,
        neural_dt_ms=1.0,
        torque_scale=3.0e-3,
        proprioceptive_mv=20.0,
        command_mv=20.0,
        settle_s=args.settle_s,
        open_loop=False,
        param=[],
    )
    loop = Loop(ns)
    x0, jac = fixed_point(loop, settle_s=args.settle_s)
    mode = pick_mode(loop, jac, "forward", args.band)
    n = loop.sizes[0]
    v_eq = x0[:n]
    print(
        f"  gait mode at the equilibrium: {mode.hz:.3f} Hz, damping {mode.zeta:+.4f}, "
        f"growth {mode.lam.real:+.3f} /s"
    )

    base = loop.runtime.network
    label = cell_groups()
    ids = loop.runtime.network.cell_ids

    def members(group: str) -> np.ndarray:
        wanted = set(GROUPS[group])
        return np.array([i for i, c in enumerate(ids) if label.get(c, c) in wanted], dtype=int)

    def apply(rows: np.ndarray, dg: float) -> None:
        loop.runtime.network = base if dg == 0.0 else with_conductance(base, rows, dg)
        # Added after proprioception, not merged into the command: ``inject``
        # sets rather than adds, and B-type cells are proprioceptive targets.
        loop.extra = {ids[i]: dg * v_eq[i] for i in rows} if dg != 0.0 else {}

    def solve(b: np.ndarray) -> np.ndarray:
        # The equilibrium does not move under this probe by construction, so the
        # shift term is identically zero; returned as such rather than computed.
        return np.zeros_like(b)

    print()
    head = (
        f"  {'group':<12} {'cells':>5} {'dgrowth per nS':>15} {'dfreq Hz/nS':>12} "
        f"{'threshold, per cell':>20}"
    )
    print(head)
    print("  " + "-" * (len(head) - 2))
    thresholds = {}
    for group in GROUPS:
        rows = members(group)
        if rows.size == 0:
            print(f"  {group:<12} {0:>5}  no cells")
            continue

        def perturb(eps: float, rows: np.ndarray = rows) -> None:
            apply(rows, eps)

        dlam = sensitivity(loop, x0, solve, mode, perturb, args.delta, 1.0e-5)
        apply(rows, 0.0)
        # Growth is Re(lambda); the gait mode sustains itself when it reaches 0.
        threshold = -mode.lam.real / dlam.real if dlam.real != 0.0 else float("nan")
        thresholds[group] = (rows, threshold)
        print(
            f"  {group:<12} {rows.size:>5} {dlam.real:+15.4f} "
            f"{dlam.imag / (2 * np.pi):+12.4f} {threshold:+17.3f} nS"
        )
    print()
    print("  threshold: the added conductance per cell at which the gait mode's growth")
    print("  rate reaches zero, first order. Negative is a source. For scale, at this")
    print("  equilibrium the B-type cells are loaded by 1.2 to 4.5 nS of gap junctions")
    print("  and VB6's measured membrane is +0.084 nS (model_assumptions 5AN.10).")

    for group in args.validate:
        if group not in thresholds:
            print(f"\n  --validate {group!r}: unknown group")
            continue
        rows, threshold = thresholds[group]
        print(f"\n  validate {group}: threshold {threshold:+.3f} nS per cell")
        for fraction in (0.5, 1.0):
            dg = fraction * threshold
            apply(rows, dg)
            residual = float(np.abs(loop.map(x0) - x0).max())
            moved_jac = jacobian(loop, x0)
            mu, vec = np.linalg.eig(moved_jac)
            lam = np.log(mu.astype(complex)) / loop.dt
            # Follow *the same* mode, by eigenvector overlap, rather than
            # re-picking by phase band -- which can switch to a different mode.
            overlap = np.abs(vec.conj().T @ mode.v) / (
                np.linalg.norm(vec, axis=0) * np.linalg.norm(mode.v)
            )
            overlap[lam.imag < 0] = 0.0
            k = int(np.argmax(overlap))
            same = lam[k]
            zeta = -same.real / abs(same)
            others = np.delete(lam, [k, int(np.argmin(np.abs(lam - same.conj())))])
            worst = others[np.argmax(others.real)]
            print(
                f"    {fraction:.1f}x ({dg:+.3f} nS): gait mode growth {same.real:+.3f} /s "
                f"(predicted {mode.lam.real + fraction * (-mode.lam.real):+.3f}), damping "
                f"{zeta:+.4f} at {same.imag / (2 * np.pi):.3f} Hz, overlap {overlap[k]:.2f}; "
                f"fastest other mode {worst.real:+.3f} /s at "
                f"{abs(worst.imag) / (2 * np.pi):.3f} Hz; residual {residual:.1e}"
            )
        apply(rows, 0.0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
