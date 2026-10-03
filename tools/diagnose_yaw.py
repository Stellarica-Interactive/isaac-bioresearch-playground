"""Does the quasi-static body yaw as it crawls, and is that physical?

    .venv/Scripts/python tools/diagnose_yaw.py

Watching the scripted gait in the viewport, the whole animal appeared to rock
back and forth in rotation rather than hold a course. It does: the body frame
swings tens of degrees every gait cycle. The question this script answers is
whether that is a real consequence of the mechanics or a defect of the solver.

The mechanism is not in doubt. Nothing external torques the animal, so as the
bending wave carries body mass around, the frame must counter-rotate to keep the
net torque zero -- the same reason a cat in free fall turns by swinging its tail.
A body carrying exactly one wavelength does that asymmetrically, because the
front and back of the wave are not mirror images across the mid-point. More
wavelengths should cancel more completely.

So the test is a sweep over wavelength. If the yaw falls as the body carries more
waves, it is counter-rotation and physical; if it does not, it is the solver.

Reported per wavelength, over the cycles AFTER the start transient has settled
(the body begins straight, and winding the wave up turns it once, permanently --
that initial turn is not yaw and must not be counted):

* yaw amplitude -- peak-to-peak swing of the rigid frame within a cycle
* course drift  -- how far the cycle-mean heading moves per cycle
* straightness  -- net displacement over path length
* speed         -- body lengths per second
"""

from __future__ import annotations

import argparse

import numpy as np

from worm.body.drag import DragParameters
from worm.body.geometry import BodyPlan
from worm.body.quasistatic import QuasiStaticBody

#: Cycles discarded before anything is measured. The body starts straight, and
#: winding the wave up onto it turns the frame once by tens of degrees. That is a
#: one-off reorientation, not a per-cycle swing, and including it was how an
#: earlier reading of this measurement reported 81 degrees against a straight
#: line fitted through a curve that jumps and then flattens.
SETTLE_CYCLES = 6


def _wave(
    body: QuasiStaticBody, t: float, *, amplitude: float, waves: float, hz: float
) -> np.ndarray:
    """A travelling wave of joint torque carrying ``waves`` wavelengths."""
    joints = np.arange(body.plan.n_joints)
    phase = 2.0 * np.pi * (hz * t - waves * joints / body.plan.n_joints)
    return amplitude * np.sin(phase)


def _run(*, waves: float, seconds: float, rate_hz: float, gait_hz: float, amplitude: float):
    body = QuasiStaticBody(BodyPlan(), drag=DragParameters())
    dt = 1.0 / rate_hz
    headings, centroids, axes, means_dir = [], [], [], []
    for step in range(int(seconds * rate_hz)):
        body.step(_wave(body, step * dt, amplitude=amplitude, waves=waves, hz=gait_hz), dt_s=dt)
        nodes = body.nodes()
        headings.append(body.heading)
        centroids.append(body.segment_centres().mean(axis=0))
        axis = nodes[-1] - nodes[0]
        axes.append(np.arctan2(axis[1], axis[0]))
        # The body's mean orientation: the honest answer to "is the whole animal
        # rotating". `heading` is segment 0's direction and segment 0 is the
        # head, so a large `heading` swing is a head swing -- the same confusion
        # that cost four wrong diagnoses in 5AE.
        mean_dir = body._directions().mean(axis=0)
        means_dir.append(np.arctan2(mean_dir[1], mean_dir[0]))
    return (
        body,
        np.array(headings),
        np.array(centroids),
        np.unwrap(np.array(axes)),
        np.unwrap(np.array(means_dir)),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--seconds", type=float, default=40.0)
    parser.add_argument("--physics-hz", type=float, default=240.0)
    parser.add_argument("--gait-hz", type=float, default=0.5, help="undulation frequency")
    parser.add_argument("--amplitude", type=float, default=2.0e-5, help="joint torque, N m")
    parser.add_argument(
        "--waves",
        type=float,
        nargs="+",
        default=[0.75, 1.0, 1.25, 1.5, 2.0],
        help="wavelengths carried on the body. A crawling C. elegans holds about "
        "1.5 (Gray and Lissmann 1964); the committed scripted gait uses 1.0.",
    )
    args = parser.parse_args()

    plan = BodyPlan()
    cycle = 1.0 / args.gait_hz
    per_cycle = int(cycle * args.physics_hz)
    print(
        f"  {plan.n_segments} segments, {plan.total_length_m * 1000:.0f} mm, "
        f"gait {args.gait_hz:g} Hz ({cycle:g} s/cycle), {args.physics_hz:g} Hz steps"
    )
    print(
        f"  discarding the first {SETTLE_CYCLES} cycles: the start transient is a "
        f"one-off reorientation, not yaw"
    )
    print()
    header = (
        f"{'waves':>6} {'head p-p':>9} {'mean dir':>9} {'end-end':>8} "
        f"{'drift/cyc':>10} {'straight':>9} {'speed':>9} {'travel':>8}"
    )
    print(header)
    print("-" * len(header))

    for waves in args.waves:
        body, headings, centroids, axes, means_dir = _run(
            waves=waves,
            seconds=args.seconds,
            rate_hz=args.physics_hz,
            gait_hz=args.gait_hz,
            amplitude=args.amplitude,
        )
        first = SETTLE_CYCLES * per_cycle
        if first + per_cycle > len(headings):
            raise SystemExit(
                f"--seconds {args.seconds} is too short to discard "
                f"{SETTLE_CYCLES} cycles and still measure one"
            )
        # Per-cycle statistics, each computed inside one cycle so a slow course
        # change cannot leak into the swing. A statistic is only meaningful if
        # the signal has the structure it assumes: the swing assumes an
        # oscillation about a level, so the level is removed cycle by cycle.
        swings, means, body_swings, mean_swings = [], [], [], []
        for start in range(first, len(headings) - per_cycle + 1, per_cycle):
            window = headings[start : start + per_cycle]
            swings.append(window.max() - window.min())
            means.append(window.mean())
            axis_window = axes[start : start + per_cycle]
            body_swings.append(axis_window.max() - axis_window.min())
            dir_window = means_dir[start : start + per_cycle]
            mean_swings.append(dir_window.max() - dir_window.min())
        drift = float(np.mean(np.diff(means))) if len(means) > 1 else 0.0

        path = float(np.sum(np.linalg.norm(np.diff(centroids[first:], axis=0), axis=1)))
        net = float(np.linalg.norm(centroids[-1] - centroids[first]))
        straightness = net / path if path > 0 else 0.0
        measured_s = (len(headings) - first) / args.physics_hz
        speed = net / plan.total_length_m / measured_s
        travel = float(np.linalg.norm(centroids[-1] - centroids[0])) / plan.total_length_m

        print(
            f"{waves:6.2f} {np.degrees(np.mean(swings)):8.1f}d "
            f"{np.degrees(np.mean(mean_swings)):8.1f}d "
            f"{np.degrees(np.mean(body_swings)):7.1f}d "
            f"{np.degrees(drift):+9.2f}d {straightness:9.3f} "
            f"{speed:8.3f}L/s {travel:8.3f}"
        )

    print()
    print("head p-p is segment 0's direction; mean dir averages all 24 segments.")
    print("A large head swing with a small mean swing is a worm waving its head,")
    print("not a worm spinning. Falling with wavelength confirms counter-rotation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
