"""How much does the body slide instead of crawling, and is the drag per-segment?

    .venv/Scripts/python tools/diagnose_slip.py

Watching the scripted gait, the body appeared to slide across the ground "as if
only calculating friction in the center of the animal". The first half of this
script checks that directly; the second puts a number on the sliding.

Slip is the standard measure for undulatory locomotion. A bending wave travels
backwards along the body at `frequency x wavelength`, and the animal advances at
some fraction of that. What it fails to advance, it slipped:

    slip = 1 - body speed / wave speed

A real *C. elegans* crawling on agar advances 70 to 90 per cent of its wave
speed, so slip is 0.1 to 0.3. Swimming in water it is far higher, which is why
the two gaits look different.

The per-segment check is worth keeping because it is cheap and it rules out a
whole class of error. Drag computed at one point could not distinguish a straight
body from a bent one, so a straight body and an L-shaped body are pulled sideways
and the net force compared against what resistive force theory says each should
give. Both agree to the digit, so the drag tensor varies along the body as it
must.

The sweep is over muscle torque, because slip is set by how deeply the body
undulates and nothing else here changes that. It is **not** monotonic: past about
23 degrees of bend the body curls rather than undulating, and slip rises again.
See model_assumptions 5AG.
"""

from __future__ import annotations

import argparse

import numpy as np

from worm.body.drag import DragParameters
from worm.body.geometry import BodyPlan
from worm.body.quasistatic import QuasiStaticBody

#: Slip of a real animal crawling on agar, for the sweep to be read against.
BIOLOGICAL_SLIP = (0.1, 0.3)


def _net_drag(body: QuasiStaticBody, direction: tuple[float, float]) -> np.ndarray:
    """Summed drag force for a rigid unit translation of the whole body."""
    velocity = np.tile(np.asarray(direction, dtype=np.float64), (body.plan.n_segments, 1))
    return -np.einsum("nab,nb->na", body._drag_tensors(), velocity).sum(axis=0)


def _per_segment_check(plan: BodyPlan, drag: DragParameters) -> None:
    """Drag at one point could not tell a straight body from a bent one."""
    n = plan.n_segments
    body = QuasiStaticBody(plan, drag=drag)
    along = float(np.linalg.norm(_net_drag(body, (1.0, 0.0))))
    across = float(np.linalg.norm(_net_drag(body, (0.0, 1.0))))

    # An L: half the segments turned 90 degrees, so a sideways pull meets half
    # the body side-on and half end-on.
    body.joint_angles = np.zeros(plan.n_joints)
    body.joint_angles[plan.n_joints // 2] = np.pi / 2.0
    bent = float(np.linalg.norm(_net_drag(body, (0.0, 1.0))))

    print("per-segment check -- net drag for a rigid 1 m/s translation:")
    print(f"  straight, along the body : {along:10.4g}   theory {n * drag.tangential:10.4g}")
    print(f"  straight, across it      : {across:10.4g}   theory {n * drag.perpendicular:10.4g}")
    print(
        f"  L-shaped, across it      : {bent:10.4g}   "
        f"theory {n / 2 * (drag.perpendicular + drag.tangential):10.4g}"
    )
    verdict = "ONE POINT ONLY -- bug" if np.isclose(bent, across) else "varies along the body, good"
    print(f"  straight vs L-shaped: {verdict}")


def _slip(
    plan: BodyPlan,
    drag: DragParameters,
    *,
    amplitude: float,
    gait_hz: float,
    waves: float,
    seconds: float,
    settle: float,
    rate_hz: float,
) -> tuple[float, float, float]:
    """Mean bend, body speed and slip, measured after the settling window."""
    body = QuasiStaticBody(plan, drag=drag)
    dt = 1.0 / rate_hz
    joints = np.arange(plan.n_joints)
    first = int(settle * rate_hz)
    bends: list[float] = []
    mark = body.segment_centres().mean(axis=0)
    for step in range(int(seconds * rate_hz)):
        phase = 2.0 * np.pi * (gait_hz * step * dt - waves * joints / plan.n_segments)
        body.step(amplitude * np.sin(phase), dt_s=dt)
        if step == first:
            mark = body.segment_centres().mean(axis=0)
        if step >= first:
            bends.append(float(np.abs(body.joint_angles).max()))
    moved = float(np.linalg.norm(body.segment_centres().mean(axis=0) - mark))
    speed = moved / plan.total_length_m / (seconds - settle)
    wave_speed = gait_hz / waves
    return float(np.degrees(np.mean(bends))), speed, 1.0 - speed / wave_speed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--seconds", type=float, default=24.0)
    parser.add_argument("--settle", type=float, default=12.0)
    parser.add_argument("--physics-hz", type=float, default=240.0)
    parser.add_argument("--gait-hz", type=float, default=0.5)
    parser.add_argument(
        "--waves",
        type=float,
        default=1.0 / 0.65,
        help="wavelengths carried on the body. The default is the runner's own "
        "wavelength_fraction of 0.65, so the sweep describes the committed gait; "
        "reading a sweep taken at one wavelength as though it did is how an "
        "earlier version of model_assumptions 5AG got it wrong.",
    )
    parser.add_argument("--drag-ratio", type=float, default=None)
    parser.add_argument(
        "--torque",
        type=float,
        nargs="+",
        default=[1e-5, 2e-5, 5e-5, 8e-5, 1e-4, 2e-4],
        help="joint torque amplitudes, N m. The runner's scripted gait works out "
        "at roughly 5e-5 equivalent; 2e-5 is the amplitude used by the tests.",
    )
    args = parser.parse_args()

    plan = BodyPlan()
    drag = DragParameters(**({} if args.drag_ratio is None else {"ratio": args.drag_ratio}))
    print(
        f"c_par {drag.tangential:.4g}  c_perp {drag.perpendicular:.4g}  "
        f"ratio {drag.perpendicular / drag.tangential:.1f}"
    )
    print()
    _per_segment_check(plan, drag)
    print()
    print(
        f"slip -- a real C. elegans on agar sits at "
        f"{BIOLOGICAL_SLIP[0]:.1f} to {BIOLOGICAL_SLIP[1]:.1f}"
    )
    print()
    header = f"{'torque':>9} {'max bend':>9} {'body spd':>11} {'wave spd':>11} {'slip':>7}"
    print(header)
    print("-" * len(header))
    for amplitude in args.torque:
        bend, speed, slip = _slip(
            plan,
            drag,
            amplitude=amplitude,
            gait_hz=args.gait_hz,
            waves=args.waves,
            seconds=args.seconds,
            settle=args.settle,
            rate_hz=args.physics_hz,
        )
        flag = "  <- biological" if BIOLOGICAL_SLIP[0] <= slip <= BIOLOGICAL_SLIP[1] else ""
        print(
            f"{amplitude:9.0e} {bend:8.1f}d {speed:10.3f}L/s "
            f"{args.gait_hz / args.waves:10.3f}L/s {slip:7.3f}{flag}"
        )
    print()
    print("Slip above the biological range means the body is sliding rather than")
    print("gripping. Linear drag has no static friction and no groove, which is")
    print("the known cause -- see worm/body/drag.py and model_assumptions 5AG.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
