"""Does a gentle touch change what the crawling worm does? Paired, standalone.

    .venv/Scripts/python tools/touch_response.py --at 0.15
    .venv/Scripts/python tools/touch_response.py --at 0.85 --command AVAL,AVAR

model_assumptions 5D found the touch circuit correctly directed but tiny -- a
quarter of a millivolt at the command interneurons, half a percent at the
muscles -- in a body that could not crawl, so there was no escape to look for.
The loop now crawls (5AN.22-5AN.24). This asks the behavioural question.

Two runs of the runner's chain, without Isaac, identical in every respect --
same build, same calibration, same start -- except that one is touched. Every
difference between them is the touch's, so no sham window is needed and a small
effect is not lost in the body's own motion, which is what defeated 5D.1b.

The touch is the runner's: each receptor's current is solved for a target
receptor potential (``--touch-mv``, 20 by default) and weighted by where the
contact falls in its receptive field (``worm/body/touch.py``). It is added after
every other input for the touch window.

Reported per window, touched minus untouched: the centroid's velocity along the
body's own tail-to-head axis (positive is head first, so a reversal is a sign
change), the mean voltage of the backward (AVA, AVD) and forward (AVB, PVC)
command interneurons, and the touched run's own axial velocity so the size of
the change can be judged against it.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, ".")
from common.neural.stimulus import solve_for_depolarisation  # noqa: E402
from tools.analyse_loop import Loop  # noqa: E402
from worm.body.touch import TouchField  # noqa: E402

COMMANDS = {
    "AVA": ("AVAL", "AVAR"),
    "AVD": ("AVDL", "AVDR"),
    "AVB": ("AVBL", "AVBR"),
    "PVC": ("PVCL", "PVCR"),
}


def build(args: argparse.Namespace) -> tuple[Loop, TouchField]:
    loop = Loop(
        argparse.Namespace(
            unknown_sign="exclude",
            physics_hz=240.0,
            neural_dt_ms=1.0,
            torque_scale=args.torque_scale,
            proprioceptive_mv=args.proprioceptive_mv,
            proprio_rate=args.proprio_rate,
            command_mv=args.command_mv,
            command_cells=args.command,
            settle_s=0.0,
            open_loop=False,
            param=args.param,
            gap_scale=args.gap_scale,
            self_contact=True,
        )
    )
    receptors = [
        c
        for c in ("ALML", "ALMR", "AVM", "PLML", "PLMR", "PVM")
        if c in loop.runtime.network.cell_ids
    ]
    per_cell = solve_for_depolarisation(loop.runtime, receptors, args.touch_mv)
    field = TouchField.build(loop.plan, cells=loop.runtime.network.cell_ids, per_cell_pa=per_cell)
    return loop, field


def run(loop: Loop, field: TouchField | None, args: argparse.Namespace):
    """Sample axial velocity and command voltages once a second."""
    ids = loop.runtime.network.cell_ids
    idx = {k: [ids.index(c) for c in v if c in ids] for k, v in COMMANDS.items()}
    per = int(round(1.0 / loop.dt))
    t_end = args.touch_start + args.touch_s + args.after_s
    rows = []
    prev = None
    for second in range(int(t_end)):
        touching = (
            field is not None and args.touch_start <= second < args.touch_start + args.touch_s
        )
        loop.extra = (
            field.currents(contact_fraction=args.at, ventral=args.ventral) if touching else {}
        )
        for _ in range(per):
            loop.step()
        c = loop.body.segment_centres()
        centroid = c.mean(axis=0)
        axis = c[0] - c[-1]
        axis /= np.linalg.norm(axis)
        axial = 0.0 if prev is None else float((centroid - prev) @ axis) / loop.plan.total_length_m
        prev = centroid.copy()
        v = loop.runtime.state[0]
        rows.append((second + 1, axial, {k: float(v[i].mean()) for k, i in idx.items()}))
    loop.extra = {}
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--at", type=float, default=0.15, help="contact position, 0 head, 1 tail")
    parser.add_argument("--ventral", action="store_true", help="include AVM / PVM")
    parser.add_argument("--touch-mv", type=float, default=20.0)
    parser.add_argument("--touch-start", type=float, default=30.0)
    parser.add_argument("--touch-s", type=float, default=2.0)
    parser.add_argument("--after-s", type=float, default=20.0)
    parser.add_argument("--command", default="AVBL,AVBR")
    parser.add_argument("--command-mv", type=float, default=20.0)
    parser.add_argument("--proprio-rate", type=float, default=0.7)
    parser.add_argument("--proprioceptive-mv", type=float, default=400.0)
    parser.add_argument("--torque-scale", type=float, default=1.0e-3)
    parser.add_argument("--param", action="append", default=[], metavar="NAME=VALUE")
    parser.add_argument("--gap-scale", action="append", default=None, metavar="CLASS=FACTOR")
    args = parser.parse_args()
    if args.gap_scale is None:
        # The measured muscle coupling is committed since 5AO.5; nothing to add.
        args.gap_scale = []

    loop_a, _ = build(args)
    loop_b, field = build(args)
    calm = run(loop_a, None, args)
    touched = run(loop_b, field, args)
    print(
        f"  touch at {args.at:g} of body length{' (ventral)' if args.ventral else ''}, "
        f"{args.touch_mv:g} mV, t = {args.touch_start:g}-{args.touch_start + args.touch_s:g} s; "
        f"command {args.command}; currents "
        f"{field.currents(contact_fraction=args.at, ventral=args.ventral)}"
    )
    head = f"  {'t':>4} {'axial BL/s':>11} {'change':>9}  " + "".join(
        f"{f'd{k} mV':>9}" for k in COMMANDS
    )
    print(head)
    print("  " + "-" * (len(head) - 2))
    for (t, ax_c, v_c), (_, ax_t, v_t) in zip(calm, touched, strict=True):
        if t < args.touch_start - 3:
            continue
        mark = " *" if args.touch_start < t <= args.touch_start + args.touch_s else "  "
        print(
            f"{mark}{t:4d} {ax_t:+11.4f} {ax_t - ax_c:+9.4f}  "
            + "".join(f"{v_t[k] - v_c[k]:+9.3f}" for k in COMMANDS)
        )
    pre = [ax for t, ax, _ in calm if args.touch_start - 10 < t <= args.touch_start]
    print(f"\n  untouched axial speed over the 10 s before: {np.mean(pre):+.4f} BL/s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
