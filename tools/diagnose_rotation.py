"""Why does a frozen body rotate? Measure the forces and momenta directly.

    C:\\isaacsim\\python.bat tools/diagnose_rotation.py

Holding the muscle drive constant produces a static shape, which generates no
thrust: anisotropic drag converts *undulation* into force, and there is no
undulation. The body should settle and stay put.

It does not. Measured over sixty seconds it holds its shape to two decimals
(`bend 5.1 deg`, `extent 0.94`) while its heading advances at a dead-linear
2.0 deg/s and it travels 0.8 body lengths. Three things say that cannot be the
drag model doing legitimate work:

* drag stops a spin in milliseconds -- the rotational time constant is well under
  a second at every body shape and every anisotropy ratio;
* the rate is constant while `amp` halves, and thrust from undulation scales as
  amplitude squared, so the rate should fall fourfold;
* internal joint torques cancel by action and reaction, so the muscles cannot
  rotate the body as a whole.

So something applies external torque. This script reports, every simulated
second: the net drag force, the net drag torque about the centroid, the body's
angular velocity, and its angular momentum. A conservative setup shows net force
and torque fluctuating about zero with angular momentum decaying. A leak shows a
persistent sign.

Written as a separate script rather than a flag on the runner because it needs to
read quantities the runner never computes, and because a diagnostic that changes
the thing it measures is worthless.
"""

from __future__ import annotations

import argparse

parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
parser.add_argument("--seconds", type=float, default=30.0)
parser.add_argument("--physics-hz", type=float, default=240.0)
parser.add_argument("--torque-scale", type=float, default=3.0e-3)
parser.add_argument(
    "--drive",
    type=float,
    default=0.05,
    help="Constant dorsal muscle activation. Uniform along the body, so a static "
    "arc is the only equilibrium it can have -- anything else is an instability.",
)
parser.add_argument(
    "--drag-ratio",
    type=float,
    default=None,
    help="Perpendicular/parallel drag anisotropy. 1 is isotropic. If the "
    "precession needs anisotropy it is a resistive-force-theory instability.",
)
parser.add_argument("--quiet", action="store_true", help="only the final precession rate")
parser.add_argument(
    "--headless",
    action="store_true",
    help="no window. Default is windowed, so the rotation can be seen as well as "
    "measured -- the two together are what identified it.",
)
args = parser.parse_args()

from isaacsim import SimulationApp  # noqa: E402

simulation_app = SimulationApp({"headless": args.headless})

import numpy as np  # noqa: E402
import omni.timeline  # noqa: E402
from isaacsim.core.experimental.prims import Articulation, RigidPrim  # noqa: E402
from isaacsim.core.simulation_manager import SimulationManager  # noqa: E402

from worm.body.drag import DragParameters, GroundDrag  # noqa: E402
from worm.body.geometry import BodyPlan  # noqa: E402
from worm.body.muscles import MuscleModel, MuscleParameters  # noqa: E402
from worm.isaac.stage import build_scene, dof_order  # noqa: E402


def main() -> int:
    dt = 1.0 / args.physics_hz
    plan = BodyPlan()
    root_path, joint_paths = build_scene(plan)
    SimulationManager.set_physics_dt(dt)
    masses = plan.masses_kg()
    articulation = Articulation(root_path)
    links = RigidPrim([f"{root_path}/segment_{i:02d}" for i in range(plan.n_segments)])
    drag = GroundDrag(
        plan,
        DragParameters(**({"ratio": args.drag_ratio} if args.drag_ratio is not None else {})),
    )
    muscle = MuscleModel(plan, MuscleParameters(peak_torque_scale=args.torque_scale))

    omni.timeline.get_timeline_interface().play()
    simulation_app.update()
    dofs = dof_order(articulation, joint_paths)
    articulation.set_dof_gains(
        stiffnesses=muscle.params.joint_stiffness,
        dampings=muscle.params.joint_damping,
    )
    # The same setup the runner does, and omitting it collapsed the body: every
    # segment piled within 4.4 mm of the centroid for a 100 mm animal, which made
    # the first reading of this diagnostic a measurement of a degenerate scene.
    articulation.set_dof_position_targets(np.zeros((1, plan.n_joints)), dof_indices=dofs)
    articulation.set_dof_armatures(2.0e-8)
    # Let the body settle into its built shape before anything is recorded.
    for _ in range(int(0.5 / dt)):
        simulation_app.update()
    spread0 = float(
        np.linalg.norm(
            np.asarray(links.get_world_poses()[0])[:, :2]
            - np.asarray(links.get_world_poses()[0])[:, :2].mean(axis=0),
            axis=1,
        ).max()
    )
    expected = 0.4 * plan.total_length_m
    if spread0 < expected:
        raise SystemExit(
            f"body is collapsed: furthest segment {spread0 * 1000:.1f} mm from the "
            f"centroid, expected at least {expected * 1000:.0f} mm for a "
            f"{plan.total_length_m * 1000:.0f} mm animal. Refusing to report "
            f"mechanics measured on a degenerate scene."
        )
    print(f"  body extends {spread0 * 1000:.1f} mm from its centroid, as built")

    # A constant, uniform dorsal bias. No neurons, no feedback, no rhythm: the
    # simplest possible drive, so anything that moves is mechanics.
    drive = np.zeros((4, plan.n_segments), dtype=np.float64)
    drive[0] = args.drive

    if not args.quiet:
        print(
            f"{'t':>6} {'heading':>9} {'omega':>10} {'net |F|':>11} {'net torque':>12} "
            f"{'L about cen':>12} {'moved':>8} {'inertia':>10} {'spread':>8}"
        )
    headings: list[tuple[float, float]] = []
    start = None
    t = 0.0
    while t < args.seconds:
        muscle.step(drive, dt_ms=dt * 1000.0)
        articulation.set_dof_efforts(muscle.joint_torques().reshape(1, -1), dof_indices=dofs)
        positions = np.asarray(links.get_world_poses()[0])
        xy = positions[:, :2]
        velocities = np.asarray(links.get_velocities()[0])
        vxy = velocities[:, :2]
        forces = np.zeros((plan.n_segments, 3), dtype=np.float32)
        planar = drag.forces(xy, vxy, dt_s=dt, masses_kg=masses)
        forces[:, :2] = planar
        links.apply_forces(forces)
        simulation_app.update()
        t += dt

        if start is None:
            start = xy.mean(axis=0)
        if round(t / dt) % round(1.0 / dt) == 0:
            centre = xy.mean(axis=0)
            r = xy - centre
            axis = positions[-1, :2] - positions[0, :2]
            heading = float(np.degrees(np.arctan2(axis[1], axis[0])))
            # Angular momentum and velocity of the body about its own centroid.
            v_rel = vxy - vxy.mean(axis=0)
            angular_momentum = float(
                np.sum(masses * (r[:, 0] * v_rel[:, 1] - r[:, 1] * v_rel[:, 0]))
            )
            inertia = float(np.sum(masses * np.sum(r**2, axis=1)))
            omega = np.degrees(angular_momentum / inertia) if inertia > 0 else 0.0
            net_force = float(np.linalg.norm(planar.sum(axis=0)))
            net_torque = float(np.sum(r[:, 0] * planar[:, 1] - r[:, 1] * planar[:, 0]))
            moved = float(np.linalg.norm(centre - start)) / plan.total_length_m
            # The body's own extent and inertia, reported so that `omega` can be
            # checked against the heading rate rather than compared with an
            # inertia computed elsewhere under different assumptions -- which is
            # how a factor of 43 crept into an earlier reading of this output.
            spread = float(np.linalg.norm(r, axis=1).max())
            # Checked every second, not only at startup. The body begins 47.9 mm
            # from its centroid and a uniform drive of 0.05 curls it to 5 mm
            # within one second -- a ball, not an animal -- and the first reading
            # of this diagnostic reported that ball's spin as a mechanical limit
            # cycle. A guard that only runs before the loop cannot catch it.
            if spread < 0.25 * plan.total_length_m:
                print(
                    f"{chr(10)}COLLAPSED at t={t:.1f}s: furthest segment "
                    f"{spread * 1000:.1f} mm from the centroid, was "
                    f"{spread0 * 1000:.1f} mm at rest. The drive has curled the "
                    f"body into a ball; nothing measured past here is mechanics "
                    f"of a worm."
                )
                simulation_app.close()
                return 1
            headings.append((t, heading))
            if not args.quiet:
                print(
                    f"{t:6.1f} {heading:+8.2f} {omega:+9.3f} {net_force:11.3e} "
                    f"{net_torque:+12.3e} {angular_momentum:+12.3e} {moved:8.3f} "
                    f"{inertia:10.3e} {spread:8.4f}"
                )

    # Unwrapped heading rate: the precession frequency, if there is one. The
    # residual says whether it is a steady rate or a decaying transient -- a limit
    # cycle fits a straight line, a transient does not.
    if len(headings) > 2:
        times = np.array([x for x, _ in headings])
        angles = np.unwrap(np.radians([y for _, y in headings]))
        fit = np.polyfit(times, angles, 1)
        rate = float(np.degrees(fit[0]))
        residual = float(np.degrees(np.std(angles - np.polyval(fit, times))))
        period = abs(360.0 / rate) if abs(rate) > 1e-9 else float("inf")
        print(
            f"\nprecession {rate:+.3f} deg/s   period {period:.0f} s   "
            f"linearity residual {residual:.2f} deg"
        )
        print(
            f"  heading rate {rate:+.3f} deg/s vs body rotation from angular "
            f"momentum: compare the omega column. If they differ, the shape is "
            f"changing as well as (or instead of) the body turning."
        )
        print(
            f"  drive={args.drive:g}  torque={args.torque_scale:g}  "
            f"ratio={args.drag_ratio if args.drag_ratio is not None else 'default'}"
        )
    print(
        "net torque holding one sign is a leak; fluctuating about zero with "
        "angular momentum conserved is not."
    )
    simulation_app.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
