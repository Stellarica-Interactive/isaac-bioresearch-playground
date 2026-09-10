"""Drive the worm body with a scripted wave. **No neurons involved.**

    C:\\isaacsim\\python.bat worm\\isaac\\run_body.py --headless --seconds 20

This is the control experiment, and it exists to answer one question before the
connectome is ever connected: *can this body crawl at all?*

If it can, and the connectome-driven version later cannot, the fault is in the
neural model. If it cannot, no neural model would have saved it. Wiring the
connectome in first would confound the two, and the failure would look
biological when it was mechanical.

The scripted sine wave here is exactly the kind of code the project forbids in a
biological result. It is legitimate only because it is a control, and it gets
deleted at the point a real claim is made.

What is being tested
--------------------

1. The articulation holds together and does not explode.
2. Anisotropic drag converts undulation into net forward displacement. Run with
   ``--drag-ratio 1`` to see the same wave produce almost no progress, which is
   the check that the anisotropy is doing the work rather than something else.
"""

from __future__ import annotations

import argparse
import sys

# --- 1. Argument parsing and app launch must happen before any Isaac import ---

parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
parser.add_argument("--headless", action="store_true", help="no viewport")
parser.add_argument("--seconds", type=float, default=20.0, help="simulated duration")
parser.add_argument("--segments", type=int, default=24)
parser.add_argument("--physics-hz", type=float, default=240.0)
parser.add_argument("--frequency-hz", type=float, default=0.5, help="undulation frequency")
parser.add_argument("--wavelength", type=float, default=0.65, help="wavelengths per body")
parser.add_argument("--torque-scale", type=float, default=None)
parser.add_argument(
    "--stiffness",
    type=float,
    default=None,
    help="passive cuticle stiffness, N m/rad. Bend amplitude is set by the "
    "ratio torque/stiffness; the absolute value sets how well the body resists "
    "buckling into a coil.",
)
parser.add_argument("--damping", type=float, default=None)
parser.add_argument(
    "--drag-ratio",
    type=float,
    default=None,
    help="perpendicular/tangential drag. Pass 1 for the isotropic control, which "
    "should barely move.",
)
parser.add_argument(
    "--armature",
    type=float,
    default=2.0e-8,
    help="extra rotational inertia per joint, kg m^2. Numerical stability aid; "
    "pass 0 to disable and watch it diverge.",
)
parser.add_argument(
    "--gravity",
    action="store_true",
    help="enable gravity and a ground plane. Off by default: the drag model IS the "
    "substrate, and adding contact friction on top resists the body twice with an "
    "isotropic force that masks the anisotropy the gait depends on.",
)
parser.add_argument(
    "--no-drag",
    action="store_true",
    help="disable ground drag entirely. Diagnostic: if the body still coils with no "
    "drag, the instability is in the muscle drive rather than drag compression.",
)
parser.add_argument(
    "--no-grid",
    action="store_true",
    help="hide the checkerboard. It is scenery with no collider, so this changes "
    "nothing a run measures.",
)
parser.add_argument(
    "--self-collision",
    action="store_true",
    help="Stop the body passing through itself when it folds. Off by default "
    "because it adds contact forces beside the drag model that stands in for "
    "the substrate, and every W2 number was measured without it.",
)
parser.add_argument("--report-every", type=float, default=2.0, help="seconds between reports")
args, _ = parser.parse_known_args()

from isaacsim import SimulationApp  # noqa: E402

simulation_app = SimulationApp({"headless": args.headless})

# --- 2. Everything else, now that the app exists ------------------------------

import numpy as np  # noqa: E402
import omni.timeline  # noqa: E402
from isaacsim.core.experimental.prims import Articulation, RigidPrim  # noqa: E402
from isaacsim.core.simulation_manager import SimulationManager  # noqa: E402

from worm.body.drag import DragParameters, GroundDrag  # noqa: E402
from worm.body.geometry import BodyPlan  # noqa: E402
from worm.body.muscles import MuscleModel, MuscleParameters, sine_wave_drive  # noqa: E402
from worm.isaac.stage import add_camera, add_ground_grid, build_scene, dof_order  # noqa: E402


def main() -> int:
    plan = BodyPlan(n_segments=args.segments)
    overrides: dict[str, float] = {}
    if args.torque_scale:
        overrides["peak_torque_scale"] = args.torque_scale
    if args.stiffness:
        overrides["joint_stiffness"] = args.stiffness
    if args.damping:
        overrides["joint_damping"] = args.damping
    muscle_params = MuscleParameters(**overrides)
    drag_params = DragParameters(
        **({"ratio": args.drag_ratio} if args.drag_ratio is not None else {})
    )

    print(plan.describe())
    muscles = MuscleModel(plan, muscle_params)
    drag = GroundDrag(plan, drag_params)
    print(f"  {drag.describe()}")
    print(
        f"  torque {muscles.params.peak_torque_scale:.3g} N m, "
        f"stiffness {muscles.params.joint_stiffness:.3g}, "
        f"target bend {
            np.degrees(muscles.params.peak_torque_scale / muscles.params.joint_stiffness):.1f} deg"
    )

    root_path, joint_paths = build_scene(
        plan,
        gravity=args.gravity,
        ground=args.gravity,
        self_collision=args.self_collision,
    )
    camera = add_camera(plan)
    if not args.no_grid:
        add_ground_grid(plan)
    print(f"  built {len(joint_paths)} joints under {root_path}\n")

    dt = 1.0 / args.physics_hz
    SimulationManager.set_physics_dt(dt)

    masses = plan.masses_kg()
    articulation = Articulation(root_path)
    segment_paths = [f"{root_path}/segment_{i:02d}" for i in range(plan.n_segments)]
    links = RigidPrim(segment_paths)

    # Physics tensors only exist once the timeline is playing; reading DOF state
    # before this point asserts rather than returning something plausible.
    omni.timeline.get_timeline_interface().play()
    simulation_app.update()

    # Passive cuticle elasticity goes into the JOINT DRIVE, not into our torque.
    #
    # This is the single most important numerical decision in the body model. The
    # effective inertia of one joint is around 3e-11 kg m^2, so a passive damping
    # of 2e-5 N m s/rad has a time constant near 1e-6 s -- roughly three thousand
    # times shorter than a 240 Hz step. Applying that damping explicitly, as part
    # of the torque we compute, is that far past its stability limit: it does not
    # damp, it oscillates and slams every joint against its limit, which is
    # exactly what it did. The body coiled into a knot and tumbled, and no torque
    # scale over three orders of magnitude changed the symptom.
    #
    # PhysX solves the joint drive implicitly, so the same stiffness and damping
    # are unconditionally stable there. It is also the more honest model: the
    # cuticle really is a passive spring-damper acting at every joint, not part of
    # the muscle's active command.
    articulation.set_dof_gains(
        stiffnesses=muscles.params.joint_stiffness,
        dampings=muscles.params.joint_damping,
    )
    # PhysX numbers DOFs from its own chosen root outward, not head to tail.
    # Every read and write below goes through this permutation; see dof_order().
    dofs = dof_order(articulation, joint_paths)
    articulation.set_dof_position_targets(np.zeros((1, plan.n_joints)), dof_indices=dofs)
    print(
        f"  passive cuticle via joint drive: stiffness "
        f"{muscles.params.joint_stiffness:.3g}, damping {muscles.params.joint_damping:.3g}"
    )

    # Joint armature: extra rotational inertia on each degree of freedom.
    #
    # These segments are tiny -- about 1e-4 kg with a transverse inertia near
    # 4e-10 kg m^2 -- so any torque large enough to overcome ground drag produces
    # an enormous angular acceleration, and an explicit solver at a few hundred Hz
    # cannot follow it. Armature is the standard robotics remedy: in a real robot
    # it represents motor rotor inertia, here it is purely numerical and is
    # recorded as an engineering addition rather than biology.
    if args.armature:
        articulation.set_dof_armatures(args.armature)
        print(f"  joint armature {args.armature:.3g} kg m^2 (numerical, not biological)")

    start = _centroid(links)
    reported = 0.0
    t_ms = 0.0
    steps = int(args.seconds * args.physics_hz)

    for step in range(steps):
        drive = sine_wave_drive(
            plan,
            t_ms,
            frequency_hz=args.frequency_hz,
            wavelength_fraction=args.wavelength,
        )
        muscles.step(drive, dt_ms=dt * 1000.0)

        # Active muscle torque only. The passive spring-damper terms are NOT
        # passed here: they are handled implicitly by the joint drive configured
        # above. Passing them would apply them twice, and explicitly, which is the
        # instability this model was built wrong with the first time.
        torque = muscles.joint_torques()
        articulation.set_dof_efforts(torque.reshape(1, -1), dof_indices=dofs)

        if not args.no_drag:
            _apply_drag(links, drag, dt, masses)

        simulation_app.update()
        t_ms += dt * 1000.0
        # Keep the animal in frame; it travels several body lengths.
        if step % 4 == 0:
            camera.follow(_centroid(links))

        if t_ms / 1000.0 - reported >= args.report_every:
            reported = t_ms / 1000.0
            _report(reported, start, links, plan, articulation=articulation, dof_indices=dofs)

    _report(
        args.seconds,
        start,
        links,
        plan,
        final=True,
        articulation=articulation,
        dof_indices=dofs,
    )
    simulation_app.close()
    return 0


def _centroid(links: RigidPrim) -> np.ndarray:
    positions, _ = links.get_world_poses()
    return np.asarray(positions)[:, :2].mean(axis=0)


def _apply_drag(links: RigidPrim, drag: GroundDrag, dt_s: float, masses_kg: np.ndarray) -> None:
    """Anisotropic ground drag, as an explicit force on each segment.

    Not left to contact friction: PhysX material friction is isotropic and cannot
    express the distinction between sliding along the body and sliding across it,
    which is the distinction that makes crawling work at all.

    ``dt_s`` and the masses are passed so the drag is computed as an exact
    exponential decay over the step rather than a raw ``-c v``. Without that the
    perpendicular coefficient is large enough relative to the segment masses to
    make explicit damping unstable at 240 Hz, and the body diverges to NaN within
    a few seconds.
    """
    positions, _ = links.get_world_poses()
    velocities = links.get_velocities()[0]
    xy = np.asarray(positions)[:, :2]
    vxy = np.asarray(velocities)[:, :2]

    planar = drag.forces(xy, vxy, dt_s=dt_s, masses_kg=masses_kg)
    forces = np.zeros((len(xy), 3), dtype=np.float32)
    forces[:, :2] = planar
    if not np.all(np.isfinite(forces)):
        raise FloatingPointError(
            "drag force is not finite; the body has already diverged. Reduce "
            "--torque-scale, or raise --physics-hz."
        )
    # Applied at each link's own transform, in the world frame: drag acts on the
    # segment as a whole, so there is no offset to specify.
    links.apply_forces(forces)


def _shape(links: RigidPrim) -> tuple[np.ndarray, float]:
    """Segment centres in the plane, and the body's heading in degrees.

    Heading is the head-to-tail direction. A crawling worm's heading drifts
    slowly; a *spinning* one's does not, and telling those apart is the
    difference between locomotion and being stirred.
    """
    positions, _ = links.get_world_poses()
    xy = np.asarray(positions)[:, :2]
    axis = xy[-1] - xy[0]
    return xy, float(np.degrees(np.arctan2(axis[1], axis[0])))


def _report(
    seconds: float,
    start: np.ndarray,
    links: RigidPrim,
    plan: BodyPlan,
    *,
    final: bool = False,
    articulation: Articulation | None = None,
    dof_indices: list[int] | None = None,
) -> None:
    now = _centroid(links)
    if not np.all(np.isfinite(now)):
        print(f"  t={seconds:5.1f}s   DIVERGED: body position is not finite")
        return
    delta = now - start
    body_lengths = float(np.linalg.norm(delta)) / plan.total_length_m
    label = "FINAL" if final else f"t={seconds:5.1f}s"

    xy, heading = _shape(links)
    # A body bent into a wave is shorter end-to-end than it is along its length.
    extent = float(np.linalg.norm(xy[-1] - xy[0])) / plan.total_length_m
    # The model is meant to be planar. If it is not, `extent` is measuring a
    # foreshortened projection rather than a bend, so report z as well.
    xyz = np.asarray(links.get_world_poses()[0])
    z_span = float(xyz[:, 2].max() - xyz[:, 2].min()) / plan.total_length_m
    extent3 = float(np.linalg.norm(xyz[-1] - xyz[0])) / plan.total_length_m

    detail = ""
    if articulation is not None:
        angles = np.degrees(
            np.asarray(articulation.get_dof_positions(dof_indices=dof_indices)).reshape(-1)
        )
        # Net turn distinguishes a travelling wave from a body winding into a
        # ring: a wave's joint angles alternate and largely cancel, a ring's do
        # not and sum toward 360 degrees.
        detail = (
            f"  bend max {np.abs(angles).max():5.1f} rms {angles.std():5.1f}"
            f"  net turn {angles.sum():+7.1f} |turn| {np.abs(angles).sum():6.1f}"
        )

    print(
        f"  {label}  moved {np.linalg.norm(delta) * 1e3:7.2f} mm "
        f"({body_lengths:5.3f} BL)  heading {heading:+7.1f}deg  "
        f"extent {extent:4.2f} (3d {extent3:4.2f}, z span {z_span:4.2f}){detail}"
    )
    if final:
        print(
            "\n  A worm crawling forward should cover a meaningful fraction of a "
            "body length.\n  Re-run with --drag-ratio 1 as a control: the same wave "
            "should barely move it,\n  which is what shows the anisotropy is doing "
            "the work."
        )


if __name__ == "__main__":
    sys.exit(main())
