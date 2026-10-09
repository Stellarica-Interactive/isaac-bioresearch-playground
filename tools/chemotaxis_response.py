"""Does a food gradient steer the crawling worm? Paired, standalone.

    .venv/Scripts/python tools/chemotaxis_response.py
    .venv/Scripts/python tools/chemotaxis_response.py --food-mv 200 --modality volatile

Real chemotaxis has two strategies. The pirouette strategy -- reverse and turn
more often when concentration is falling (Pierce-Shimomura, Morse & Lockery
1999) -- needs reversals, which model_assumptions 5AO finds the model cannot
choose to make. The weathervane strategy -- curving toward the side where each
head swing finds more odour -- needs only forward crawling, which the loop now
has. So that is what this tests.

Three runs of the runner's chain, identical until the food appears: crawl for
``--warmup-s``, then place a food source ``--distance`` body lengths from the
head, ``--angle`` degrees to the left in one run and to the right in another,
and nowhere in the third. The chemosensors are ``worm/body/chemotaxis.py``'s,
calibrated as the runner calibrates them. A weathervane turn shows as the
heading change of the left run exceeding the control's and the right run's
falling short of it; anything the two share is the gait's own drift.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, ".")
from common.neural.stimulus import solve_for_depolarisation  # noqa: E402
from tools.analyse_loop import Loop  # noqa: E402
from worm.body.chemotaxis import (  # noqa: E402
    CHEMOSENSORS,
    DEFAULT_DECAY_BODY_LENGTHS,
    Chemosensation,
    FoodSource,
)


def build(args):
    loop = Loop(
        argparse.Namespace(
            unknown_sign="exclude",
            physics_hz=240.0,
            neural_dt_ms=1.0,
            torque_scale=args.torque_scale,
            proprioceptive_mv=args.proprioceptive_mv,
            proprio_rate=args.proprio_rate,
            command_mv=20.0,
            settle_s=0.0,
            open_loop=False,
            param=args.param,
            gap_scale=args.gap_scale,
            self_contact=True,
        )
    )
    cells = [s.cell for s in CHEMOSENSORS if s.cell in loop.runtime.network.cell_ids]
    per_cell = solve_for_depolarisation(loop.runtime, cells, args.food_mv)
    sense = Chemosensation.build(
        cells=loop.runtime.network.cell_ids, modality=args.modality, per_cell_pa=per_cell
    )
    return loop, sense


def heading(loop: Loop) -> float:
    c = loop.body.segment_centres()
    axis = c[0] - c[-1]
    return float(np.arctan2(axis[1], axis[0]))


def run(args, side: float | None):
    """Heading (unwrapped, degrees) every second after the food appears."""
    loop, sense = build(args)
    per = int(round(1.0 / loop.dt))
    for _ in range(int(args.warmup_s) * per):
        loop.step()
    food = None
    if side is not None:
        h = heading(loop)
        bearing = h + np.radians(side * args.angle)
        nose = loop.body.segment_centres()[0]
        d = args.distance * loop.plan.total_length_m
        food = FoodSource(
            x_m=float(nose[0] + d * np.cos(bearing)),
            y_m=float(nose[1] + d * np.sin(bearing)),
            decay_m=DEFAULT_DECAY_BODY_LENGTHS * loop.plan.total_length_m,
        )
    out = [heading(loop)]
    deviations = []
    for _ in range(int(args.run_s)):
        for _ in range(per):
            if food is not None:
                sense.step(food.concentration(loop.body.segment_centres()[0]), loop.dt * 1000.0)
                loop.extra = sense.currents()
                deviations.append(sense.deviation)
            loop.step()
        out.append(heading(loop))
    loop.extra = {}
    unwrapped = np.degrees(np.unwrap(np.array(out)))
    return unwrapped - unwrapped[0], (np.abs(deviations).max() if deviations else 0.0)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--distance", type=float, default=3.0, help="body lengths")
    parser.add_argument("--angle", type=float, default=45.0, help="degrees off the heading")
    parser.add_argument("--food-mv", type=float, default=20.0)
    parser.add_argument("--modality", default="awc")
    parser.add_argument("--warmup-s", type=float, default=20.0)
    parser.add_argument("--run-s", type=float, default=60.0)
    parser.add_argument("--proprio-rate", type=float, default=0.7)
    parser.add_argument("--proprioceptive-mv", type=float, default=400.0)
    parser.add_argument("--torque-scale", type=float, default=1.0e-3)
    parser.add_argument("--param", action="append", default=[], metavar="NAME=VALUE")
    parser.add_argument("--gap-scale", action="append", default=None, metavar="CLASS=FACTOR")
    args = parser.parse_args()
    if args.gap_scale is None:
        # The measured muscle coupling is committed since 5AO.5; nothing to add.
        args.gap_scale = []

    control, _ = run(args, None)
    left, dev_l = run(args, +1.0)
    right, dev_r = run(args, -1.0)
    print(
        f"  food {args.distance:g} BL away, {args.angle:g} deg left / right of the heading; "
        f"{args.modality} pathway at {args.food_mv:g} mV full scale; "
        f"largest concentration change sensed {max(dev_l, dev_r):.2e}"
    )
    print(f"  {'t s':>4} {'control':>9} {'left':>9} {'right':>9} {'left-right':>11}")
    for k in range(0, len(control), 10):
        print(
            f"  {k:4d} {control[k]:+9.2f} {left[k]:+9.2f} {right[k]:+9.2f} "
            f"{left[k] - right[k]:+11.2f}"
        )
    print("  heading change since the food appeared, degrees; positive is to the left.")
    print("  A weathervane turn makes left-right grow positive; zero means the food did nothing.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
