"""Does the articulation end up where the quasi-static solver put it?

    C:\\isaacsim\\python.bat tools/check_kinematic_drive.py

Under ``--quasistatic`` the solver owns the configuration and the articulation is
meant to be a display of it: the runner writes the root pose and every joint
angle, and PhysX integrates nothing. If that write does not reproduce the solved
body, every number the runner reports is measured on a different animal from the
one the solver computed -- and the viewport shows that different animal.

There is reason to think it does not. The solver's segment-centre chord swings
12.6 degrees per gait cycle at 1.82 waves; the same column read out of the
articulation swings 92, which is close to the 104 the body's shape would give
with no counter-rotation at all. Two explanations have already been ruled out by
measurement: PhysX's joint drive (disabling it changed the output not at all) and
the quaternion convention (``tools/check_root_pose.py`` confirms w-x-y-z is
right).

So this compares the two bodies directly, segment by segment, instead of
comparing summary statistics of each. It reports the worst position error, where
along the body it falls, and the two chord angles side by side.
"""

from __future__ import annotations

import argparse

parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
parser.add_argument("--seconds", type=float, default=6.0)
parser.add_argument("--physics-hz", type=float, default=240.0)
parser.add_argument("--wavelength", type=float, default=0.55)
parser.add_argument("--drag-ratio", type=float, default=40.0)
parser.add_argument(
    "--headless",
    action="store_true",
    help="no window. Default is windowed; this prints numbers rather than showing "
    "motion, but the standing rule is that Isaac runs are watched.",
)
args = parser.parse_args()

from isaacsim import SimulationApp  # noqa: E402

simulation_app = SimulationApp({"headless": args.headless})

import numpy as np  # noqa: E402
import omni.timeline  # noqa: E402
from isaacsim.core.experimental.prims import Articulation, RigidPrim  # noqa: E402
from isaacsim.core.simulation_manager import SimulationManager  # noqa: E402

from worm.body.drag import DragParameters  # noqa: E402
from worm.body.geometry import BodyPlan  # noqa: E402
from worm.body.muscles import (  # noqa: E402
    MuscleModel,
    MuscleParameters,
    sine_wave_drive,
)
from worm.body.quasistatic import QuasiStaticBody  # noqa: E402
from worm.isaac.stage import build_scene, dof_order  # noqa: E402

SEED_GAIN = 5.0e-5 / 3.0e-3


def main() -> int:
    dt = 1.0 / args.physics_hz
    plan = BodyPlan()
    root_path, joint_paths = build_scene(plan)
    SimulationManager.set_physics_dt(dt)
    articulation = Articulation(root_path)
    links = RigidPrim([f"{root_path}/segment_{i:02d}" for i in range(plan.n_segments)])

    drag = DragParameters(ratio=args.drag_ratio, exponent=0.6, yield_force=2.0e-5)
    body = QuasiStaticBody(plan, drag=drag, muscle=MuscleParameters(peak_torque_scale=3.0e-3))
    muscle = MuscleModel(plan, MuscleParameters(peak_torque_scale=3.0e-3))

    omni.timeline.get_timeline_interface().play()
    simulation_app.update()
    dofs = dof_order(articulation, joint_paths)
    articulation.set_dof_gains(
        stiffnesses=np.zeros((1, plan.n_joints)), dampings=np.zeros((1, plan.n_joints))
    )

    # Which segment PhysX treats as the root: straighten the joints, write a
    # known pose, see which one lands on it.
    articulation.set_dof_positions(np.zeros((1, plan.n_joints), dtype=np.float32), dof_indices=dofs)
    articulation.set_world_poses(
        positions=np.zeros((1, 3), dtype=np.float32),
        orientations=np.array([[1.0, 0.0, 0.0, 0.0]], dtype=np.float32),
    )
    simulation_app.update()
    probe_xy = np.asarray(links.get_world_poses()[0])[:, :2]
    root_link = int(np.argmin(np.linalg.norm(probe_xy, axis=1)))
    print(f"articulation root link: segment {root_link} of {plan.n_segments}")
    print()
    print(f"{'t':>6} {'worst err':>11} {'at seg':>7} {'solver chord':>13} {'isaac chord':>12}")
    print("-" * 54)
    worst_overall = 0.0
    worst_seg = -1
    for step in range(int(args.seconds * args.physics_hz)):
        drive = sine_wave_drive(
            plan,
            step * dt * 1000.0,
            frequency_hz=0.5,
            wavelength_fraction=args.wavelength,
        )
        muscle.step(SEED_GAIN * drive, dt_ms=dt * 1000.0)
        body.step(muscle.joint_torques(), dt_s=dt)

        root = body.segment_centres()[root_link]
        angle = body.segment_angles()[root_link]
        articulation.set_world_poses(
            positions=np.array([[*root, 0.0]], dtype=np.float32),
            orientations=np.array(
                [[np.cos(angle / 2.0), 0.0, 0.0, np.sin(angle / 2.0)]], dtype=np.float32
            ),
        )
        articulation.set_dof_positions(
            body.joint_angles.reshape(1, -1).astype(np.float32), dof_indices=dofs
        )
        simulation_app.update()

        if step % int(0.25 * args.physics_hz):
            continue
        solved = body.segment_centres()
        shown = np.asarray(links.get_world_poses()[0])[:, :2]
        error = np.linalg.norm(shown - solved, axis=1)
        if error.max() > worst_overall:
            worst_overall, worst_seg = float(error.max()), int(error.argmax())
        s_axis = solved[-1] - solved[0]
        i_axis = shown[-1] - shown[0]
        print(
            f"{step * dt:6.2f} {error.max() * 1e3:8.3f}mm {int(error.argmax()):7d} "
            f"{np.degrees(np.arctan2(s_axis[1], s_axis[0])):12.2f}d "
            f"{np.degrees(np.arctan2(i_axis[1], i_axis[0])):11.2f}d"
        )

    print()
    print(f"worst error over the run: {worst_overall * 1e3:.3f} mm at segment {worst_seg}")
    print(f"  (a segment is {plan.segment_length_m * 1e3:.2f} mm long)")
    if worst_overall > 0.5 * plan.segment_length_m:
        print("  -> the articulation is NOT showing the solved body")
    else:
        print("  -> the articulation matches the solver; the chords must agree too")
    simulation_app.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
