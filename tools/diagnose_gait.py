"""Is the gait the right KIND of motion, not just the right distance?

    .venv/Scripts/python tools/diagnose_gait.py
    .venv/Scripts/python tools/diagnose_gait.py --sweep wavelength

Reported from the viewport, about the scripted gait: "besides the fact that the
animal moves, it happens incorrectly... The whole rotation is also sinusoidal.
For this scripted movement I expect full straight movement."

Every aggregate said the gait was fine -- path straightness 0.997, a clean
head-to-tail wave, a biological speed. They were all true and none of them could
see what was wrong, because none of them is about the *shape* of the motion.
Two that are:

**Track following.** Real undulatory crawling on agar slithers: each body point
passes through roughly where the point ahead of it was, so the animal travels
along a groove its own head cut. A body that sweeps sideways through the medium
covers the same ground with the same aggregates. Measured as the distance from
each tail position to the nearest *earlier* head position -- earlier matters, and
comparing against the whole track regardless of time was how this measurement
first reported a spurious 52 mm.

**The rotation, decomposed.** The nose-to-tail chord of an undulating body swings
even when the body does not rotate at all, because its endpoints are riding a
wave. So the chord is measured three ways: with the body's real heading, with the
heading pinned at zero (shape only), and as the mean of every segment's
direction. The shape contributes about 119 degrees, the rigid counter-rotation
that zero net torque demands removes about 115 of it, and some 25 degrees are
left -- which the mean direction confirms, so the animal really does rotate.

See model_assumptions 5AI.
"""

from __future__ import annotations

import argparse

import numpy as np

from worm.body.drag import DragParameters
from worm.body.geometry import BodyPlan
from worm.body.muscles import MuscleModel, MuscleParameters, sine_wave_drive
from worm.body.quasistatic import QuasiStaticBody

#: What the runner applies to the scripted wave, cancelling ``--torque-scale`` so
#: the control always runs at the muscle model's default peak torque.
SEED_GAIN = 5.0e-5 / 3.0e-3

#: Cycles discarded before measuring. The body starts straight and winding the
#: wave onto it reorients the frame once; that is not part of the gait.
SETTLE_CYCLES = 6


def _chord_angle(joint_angles: np.ndarray, heading: float) -> float:
    """Nose-to-tail direction of a body with these angles and this heading."""
    angles = heading + np.concatenate([[0.0], np.cumsum(joint_angles)])
    chord = np.stack([np.cos(angles), np.sin(angles)], axis=1).sum(axis=0)
    return float(np.arctan2(chord[1], chord[0]))


def _per_cycle_swing(series: list[float], per_cycle: int) -> float:
    """Mean within-cycle peak-to-peak, degrees.

    Within a cycle, so a slow course change cannot leak into the swing: the
    statistic assumes an oscillation about a level, so the level is removed
    cycle by cycle.
    """
    values = np.unwrap(np.asarray(series, dtype=np.float64))
    windows = range(0, len(values) - per_cycle + 1, per_cycle)
    return float(
        np.degrees(
            np.mean(
                [values[i : i + per_cycle].max() - values[i : i + per_cycle].min() for i in windows]
            )
        )
    )


def measure(
    *,
    drag: DragParameters,
    taper: float,
    wavelength: float,
    gait_hz: float,
    seconds: float,
    rate_hz: float,
    groove: bool = False,
    self_contact: bool = False,
) -> dict[str, float]:
    """Run the runner's own drive chain and report the shape of the motion.

    ``groove`` and ``self_contact`` are passed to the constructor rather than set
    on the class: they are dataclass fields, so their defaults are baked into
    ``__init__`` when the class is created and assigning to the class attribute
    afterwards changes nothing. A sweep that did that reported four identical
    rows and looked like a null result.
    """
    plan = BodyPlan()
    body = QuasiStaticBody(plan, drag=drag, groove=groove, self_contact=self_contact)
    muscle = MuscleModel(plan, MuscleParameters(peak_torque_scale=3.0e-3))
    dt = 1.0 / rate_hz
    dt_ms = 1000.0 / rate_hz
    per_cycle = int(rate_hz / gait_hz)
    first = SETTLE_CYCLES * per_cycle

    heads, tails, centroids = [], [], []
    total, shape, rigid, mean_dir, bends = [], [], [], [], []
    for step in range(int(seconds * rate_hz)):
        drive = sine_wave_drive(
            plan,
            step * dt_ms,
            frequency_hz=gait_hz,
            wavelength_fraction=wavelength,
            taper=taper,
        )
        muscle.step(SEED_GAIN * drive, dt_ms=dt_ms)
        body.step(muscle.joint_torques(), dt_s=dt)
        if step < first:
            continue
        nodes = body.nodes()
        heads.append(nodes[0].copy())
        tails.append(nodes[-1].copy())
        centroids.append(body.segment_centres().mean(axis=0))
        total.append(_chord_angle(body.joint_angles, body.heading))
        shape.append(_chord_angle(body.joint_angles, 0.0))
        rigid.append(body.heading)
        direction = body._directions().mean(axis=0)
        mean_dir.append(float(np.arctan2(direction[1], direction[0])))
        bends.append(float(np.abs(body.joint_angles).max()))

    heads_a, tails_a = np.array(heads), np.array(tails)
    centroids_a = np.array(centroids)

    # Only head positions that precede the tail sample. The head is roughly a
    # body length ahead of the tail, so a generous lead is dropped.
    lead = min(int(8.0 * rate_hz), len(tails_a) // 2)
    stride = max(len(heads_a) // 400, 1)
    track = [
        float(np.min(np.linalg.norm(heads_a[: i + 1 : stride] - tails_a[i], axis=1)))
        for i in range(lead, len(tails_a), stride)
    ]

    relative = centroids_a - centroids_a[0]
    direction = float(np.arctan2(relative[-1][1], relative[-1][0]))
    across = np.array([-np.sin(direction), np.cos(direction)])
    lateral = relative @ across
    head_lateral = (heads_a - heads_a[0]) @ across
    measured_s = len(centroids_a) / rate_hz
    plan_length = plan.total_length_m

    return {
        "track": float(np.median(track)) * 1e3,
        "track95": float(np.percentile(track, 95)) * 1e3,
        "shape": _per_cycle_swing(shape, per_cycle),
        "rigid": _per_cycle_swing(rigid, per_cycle),
        "net": _per_cycle_swing(total, per_cycle),
        "mean": _per_cycle_swing(mean_dir, per_cycle),
        "head": float(np.ptp(head_lateral)) * 1e3,
        "wobble": float(lateral.max() - lateral.min()) * 1e3,
        "speed": float(np.linalg.norm(relative[-1])) / plan_length / measured_s,
        "bend": float(np.degrees(np.mean(bends))),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--seconds", type=float, default=28.0)
    parser.add_argument("--physics-hz", type=float, default=120.0)
    parser.add_argument("--gait-hz", type=float, default=0.5)
    parser.add_argument("--wavelength", type=float, default=0.65)
    parser.add_argument("--taper", type=float, default=0.0)
    parser.add_argument("--drag-exponent", type=float, default=1.0)
    parser.add_argument("--drag-yield", type=float, default=0.0)
    parser.add_argument(
        "--sweep",
        choices=("none", "wavelength", "taper"),
        default="none",
        help="sweep one parameter instead of reporting a single configuration",
    )
    args = parser.parse_args()

    drag = DragParameters(exponent=args.drag_exponent, yield_force=args.drag_yield)
    common = {
        "drag": drag,
        "gait_hz": args.gait_hz,
        "seconds": args.seconds,
        "rate_hz": args.physics_hz,
    }
    label = (
        f"drag exponent {drag.exponent:g}, yield {drag.yield_force:g} N"
        if not drag.is_linear
        else "linear drag (the committed law)"
    )

    if args.sweep == "none":
        r = measure(taper=args.taper, wavelength=args.wavelength, **common)
        print(f"{label}, wavelength {args.wavelength:g} BL, taper {args.taper:g}")
        print()
        print(f"  bend                       {r['bend']:8.1f} deg")
        print(f"  tail off the head's track  {r['track']:8.2f} mm   (95th {r['track95']:.1f})")
        print(f"  head excursion             {r['head']:8.2f} mm")
        print(f"  centroid lateral wobble    {r['wobble']:8.2f} mm")
        print(f"  speed                      {r['speed']:8.3f} BL/s")
        print()
        print("  rotation, decomposed per cycle:")
        print(f"    from the shape alone     {r['shape']:8.2f} deg")
        print(f"    rigid counter-rotation   {r['rigid']:8.2f} deg")
        print(f"    net chord swing          {r['net']:8.2f} deg")
        print(f"    mean of all directions   {r['mean']:8.2f} deg")
        removed = 1.0 - r["net"] / r["shape"] if r["shape"] else 0.0
        print(f"  -> the counter-rotation removes {removed * 100:.0f}% of the shape's swing")
        print()
        print("  The mean direction is the honest answer to whether the whole")
        print("  animal turns; the chord alone cannot tell rotation from a wave.")
        return 0

    sweeps = {
        "wavelength": (0.9, 0.8, 0.65, 0.55, 0.45, 0.35),
        "taper": (0.0, 0.08, 0.17, 0.25, 0.5),
    }
    name = "waves" if args.sweep == "wavelength" else "taper"
    print(f"{label}, sweeping {args.sweep}")
    print()
    header = (
        f"{name:>7} {'bend':>7} {'net rot':>9} {'mean rot':>9} "
        f"{'track':>8} {'wobble':>8} {'speed':>10}"
    )
    print(header)
    print("-" * len(header))
    for value in sweeps[args.sweep]:
        if args.sweep == "wavelength":
            r = measure(taper=args.taper, wavelength=value, **common)
            shown = 1.0 / value
        else:
            r = measure(taper=value, wavelength=args.wavelength, **common)
            shown = value
        print(
            f"{shown:7.2f} {r['bend']:6.1f}d {r['net']:8.1f}d {r['mean']:8.1f}d "
            f"{r['track']:7.2f}mm {r['wobble']:7.2f}mm {r['speed']:8.3f}L/s"
        )
    print()
    print("A mechanism that reduces the rotation only by reducing the bend has")
    print("not fixed anything: 5AG shows slip rises steeply as the bend drops.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
