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
    drag = GroundDrag(plan, DragParameters())
    muscle = MuscleModel(plan, MuscleParameters(peak_torque_scale=args.torque_scale))

    omni.timeline.get_timeline_interface().play()
    simulation_app.update()
    dofs = dof_order(articulation, joint_paths)
    articulation.set_dof_gains(
        stiffnesses=muscle.params.joint_stiffness,
        dampings=muscle.params.joint_damping,
    )

    # A constant, uniform dorsal bias. No neurons, no feedback, no rhythm: the
    # simplest possible drive, so anything that moves is mechanics.
    drive = np.zeros((4, plan.n_segments), dtype=np.float64)
    drive[0] = 0.05

    print(
        f"{'t':>6} {'heading':>9} {'omega':>10} {'net |F|':>11} {'net torque':>12} "
        f"{'L about cen':>12} {'moved':>8}"
    )
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
            print(
                f"{t:6.1f} {heading:+8.2f} {omega:+9.3f} {net_force:11.3e} "
                f"{net_torque:+12.3e} {angular_momentum:+12.3e} {moved:8.3f}"
            )

    print(
        "\nnet torque holding one sign is a leak; fluctuating about zero with "
        "angular momentum decaying is conservative."
    )
    simulation_app.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
