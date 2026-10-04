"""Does the loop pass a bend head-to-tail more than tail-to-head, and at what frequency?

    .venv/Scripts/python tools/loop_reciprocity.py
    .venv/Scripts/python tools/loop_reciprocity.py --open-loop
    .venv/Scripts/python tools/loop_reciprocity.py --proprio-rate 0.75

A travelling wave needs coupling along the body that is *non-reciprocal*: a bend
at one joint must drive the joints behind it more than the joints ahead, at the
frequency the wave runs at. Gap junctions and a passive body are reciprocal by
construction. model_assumptions 5AN.12 proposed that this loop is dominated by
reciprocal coupling at the gait frequency, which would explain why every source
of energy it was given fed a standing mode.

This measures it. At the loop's equilibrium, the joint-angle block of the
resolvent of the one-step map,

    R(z) = [(z I - J)^-1]_{joints, joints},   z = exp(i 2 pi f dt),

is the bend that appears at joint ``j`` when joint ``i`` is nudged, for a nudge
repeating at frequency ``f``. For each distance ``d`` it compares the response
``d`` joints *posterior* to the nudge with the response ``d`` joints *anterior*.
A ratio of 1 is reciprocal; above 1, the loop passes bends tailward, which is
the direction a forward-crawling wave travels.

``--open-loop`` cuts proprioception: the control, which a reciprocal model
predicts to be exactly 1 at every distance and frequency.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, ".")
from tools.analyse_loop import Loop, fixed_point  # noqa: E402


def transmission(loop: Loop, jac: np.ndarray, hz: float, max_d: int):
    """Mean |posterior| and |anterior| response at each distance, and their ratio."""
    n, _, m, j = loop.sizes
    q = np.arange(2 * n + m, 2 * n + m + j)
    z = np.exp(2j * np.pi * hz * loop.dt)
    rhs = np.zeros((jac.shape[0], j), dtype=complex)
    rhs[q, np.arange(j)] = 1.0
    r = np.linalg.solve(z * np.eye(jac.shape[0]) - jac, rhs)[q, :]  # [response, nudged]
    out = []
    for d in range(1, max_d + 1):
        post_c = np.array([r[i + d, i] / r[i, i] for i in range(j - d)])
        ant_c = np.array([r[i - d, i] / r[i, i] for i in range(d, j)])
        post = np.abs([r[i + d, i] for i in range(j - d)]).mean()
        ant = np.abs([r[i - d, i] for i in range(d, j)]).mean()
        # Phase of the response d joints behind the nudge minus d joints ahead,
        # each relative to the nudged joint's own response. Negative: the bend
        # arrives behind later than ahead, so bends move tailward -- the
        # direction of a forward-crawling wave. Positive: headward, backward.
        shift = float(np.degrees(np.angle(post_c.mean()) - np.angle(ant_c.mean())))
        shift = (shift + 180.0) % 360.0 - 180.0
        out.append((d, post, ant, post / ant, shift))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--frequencies", type=float, nargs="+", default=[0.0, 0.17, 0.5])
    parser.add_argument("--distances", type=int, default=6)
    parser.add_argument("--open-loop", action="store_true")
    parser.add_argument("--proprioceptive-mv", type=float, default=20.0)
    parser.add_argument("--proprio-rate", type=float, default=0.0)
    parser.add_argument("--settle-s", type=float, default=10.0)
    parser.add_argument("--param", action="append", default=[], metavar="NAME=VALUE")
    parser.add_argument(
        "--start",
        metavar="NPY",
        help="an equilibrium to continue from (a saved state vector), for settings "
        "where the trajectory swings too far for its average to be a good guess",
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
        open_loop=args.open_loop,
        param=args.param,
    )
    loop = Loop(ns)
    start = np.load(args.start) if args.start else None
    x0, jac = fixed_point(loop, settle_s=args.settle_s, start=start)
    bend = float(np.degrees(np.abs(loop.joints(x0)).max()))
    print(
        f"  {'proprioception cut' if args.open_loop else 'closed loop'}, "
        f"rate fraction {args.proprio_rate:g}; equilibrium bend {bend:.1f} deg"
    )
    print()
    head = "  " + f"{'distance':>9}" + "".join(f"{f'{f:g} Hz':>10}" for f in args.frequencies)
    print("  posterior / anterior response to a nudged joint")
    print(head)
    print("  " + "-" * (len(head) - 2))
    tables = [transmission(loop, jac, f, args.distances) for f in args.frequencies]
    for k in range(args.distances):
        print(f"  {k + 1:>9}" + "".join(f"{t[k][3]:10.3f}" for t in tables))
    print()
    print("  phase behind minus phase ahead, degrees (negative: bends move tailward)")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for k in range(args.distances):
        print(f"  {k + 1:>9}" + "".join(f"{t[k][4]:+10.1f}" for t in tables))
    print()
    print("  1.000 is reciprocal. Above 1 the loop passes bends tailward, the")
    print("  direction a forward wave travels; a gait needs that at its own frequency.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
