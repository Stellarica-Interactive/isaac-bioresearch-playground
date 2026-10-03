"""Does the articulation's root orientation take effect as written?

    C:\\isaacsim\\python.bat tools/check_root_pose.py

Under ``--quasistatic`` the runner owns the configuration and writes it into the
articulation: the root's world pose, then every joint angle. If the root
orientation does not take effect as intended, the body still has the right
*shape* but is not turned to face the right way -- and because the solver's
heading swings about 100 degrees per gait cycle, counter-rotating against the
shape's own swing to leave only 10 or 20 degrees net, losing it means the full
shape swing shows up on screen instead.

That is the symptom this exists to rule in or out. Measured standalone, the
segment-centre chord of the committed gait swings 12.6 degrees per cycle at 1.82
waves; the same column read out of Isaac showed 92, which is close to the
shape-only term of 104. Something is eating the root rotation.

No simulation: build the scene, write a known pose with the joints held
straight, read the segment positions back, and compare against arithmetic. A
straight body at heading ``h`` must lie along ``(cos h, sin h)``.
"""

from __future__ import annotations

import argparse

parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
parser.add_argument("--heading-deg", type=float, default=40.0)
parser.add_argument("--headless", action="store_true")
args = parser.parse_args()

from isaacsim import SimulationApp  # noqa: E402

simulation_app = SimulationApp({"headless": args.headless})

import numpy as np  # noqa: E402
import omni.timeline  # noqa: E402
from isaacsim.core.experimental.prims import Articulation, RigidPrim  # noqa: E402
from isaacsim.core.simulation_manager import SimulationManager  # noqa: E402

from worm.body.geometry import BodyPlan  # noqa: E402
from worm.isaac.stage import build_scene, dof_order  # noqa: E402


def main() -> int:
    plan = BodyPlan()
    root_path, joint_paths = build_scene(plan)
    SimulationManager.set_physics_dt(1.0 / 240.0)
    articulation = Articulation(root_path)
    links = RigidPrim([f"{root_path}/segment_{i:02d}" for i in range(plan.n_segments)])

    omni.timeline.get_timeline_interface().play()
    simulation_app.update()
    dofs = dof_order(articulation, joint_paths)
    # Nothing driving the joints, so the body cannot bend on its own.
    articulation.set_dof_gains(
        stiffnesses=np.zeros((1, plan.n_joints)), dampings=np.zeros((1, plan.n_joints))
    )

    heading = np.radians(args.heading_deg)
    half = heading / 2.0
    conventions = {
        "wxyz": [np.cos(half), 0.0, 0.0, np.sin(half)],
        "xyzw": [0.0, 0.0, np.sin(half), np.cos(half)],
    }

    print(f"asking for heading {args.heading_deg:+.1f} deg with the joints held straight")
    print(f"a straight body must then lie along ({np.cos(heading):+.3f}, {np.sin(heading):+.3f})")
    print()
    for name, quat in conventions.items():
        articulation.set_world_poses(
            positions=np.array([[0.0, 0.0, 0.0]], dtype=np.float32),
            orientations=np.array([quat], dtype=np.float32),
        )
        articulation.set_dof_positions(
            np.zeros((1, plan.n_joints), dtype=np.float32), dof_indices=dofs
        )
        for _ in range(4):
            simulation_app.update()

        xyz = np.asarray(links.get_world_poses()[0])
        axis = xyz[-1, :2] - xyz[0, :2]
        got = float(np.degrees(np.arctan2(axis[1], axis[0])))
        out_of_plane = float(np.abs(xyz[:, 2]).max())
        error = abs(((got - args.heading_deg + 180.0) % 360.0) - 180.0)
        verdict = "CORRECT" if error < 2.0 else f"WRONG by {error:.1f} deg"
        print(
            f"  passing {name:5s} -> body lies at {got:+8.2f} deg, "
            f"max |z| {out_of_plane * 1e3:6.2f} mm   {verdict}"
        )

    print()
    print("Isaac's own convention decides which of these the runner must pass.")
    print("A large max |z| means the rotation went out of plane, which the planar")
    print("constraint then fights -- leaving the body unturned.")
    simulation_app.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
