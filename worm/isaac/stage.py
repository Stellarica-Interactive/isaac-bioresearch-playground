"""Author the worm body as a USD articulation.

Split from :mod:`worm.body` deliberately: everything here needs ``pxr`` and a
running Isaac Sim, everything there does not. The body plan is the contract
between them, so the biomechanics stay testable in milliseconds while only this
file requires a simulator.

The articulation is a serial chain: ``n`` capsule links joined by ``n-1``
revolute joints, all rotating about the world Z axis so the body bends in the
XY plane. See :mod:`worm.body.geometry` for why the model is planar and why
there are 24 segments.

Drive mode
----------

Joints are created with **zero stiffness and damping** so PhysX applies no
position control of its own. All torque comes from the muscle model via
``set_dof_efforts``. Leaving a default drive in place would mean an invisible
PD controller quietly fighting the muscles, which is exactly the kind of thing
that makes a locomotion result meaningless without ever looking wrong.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from worm.body.geometry import BodyPlan

if TYPE_CHECKING:  # pragma: no cover - only available inside Isaac Sim
    from pxr import Usd

DEFAULT_ROOT = "/World/Worm"

#: Joint limit in degrees. Real worms bend hard, but a limit stops the chain
#: folding through itself when a torque is misconfigured.
DEFAULT_JOINT_LIMIT_DEG = 60.0


def build_worm(
    stage: Usd.Stage,
    plan: BodyPlan | None = None,
    *,
    root_path: str = DEFAULT_ROOT,
    joint_limit_deg: float = DEFAULT_JOINT_LIMIT_DEG,
    z_offset_m: float | None = None,
    planar: bool = True,
    self_collision: bool = False,
) -> tuple[str, list[str]]:
    """Create the articulation on ``stage``. Returns ``(root_path, joint_paths)``.

    ``joint_paths`` is ordered head to tail, which is the order every torque and
    angle array in this package uses.
    """
    from pxr import Gf, PhysxSchema, Sdf, UsdGeom, UsdPhysics

    plan = plan or BodyPlan()
    segments = plan.segments()
    masses = plan.masses_kg()
    # Rest on the ground plane rather than intersecting it.
    z = plan.max_radius_m if z_offset_m is None else z_offset_m

    root = UsdGeom.Xform.Define(stage, root_path)
    UsdPhysics.ArticulationRootAPI.Apply(root.GetPrim())
    # A worm is not bolted down: it is a floating-base articulation.
    #
    # Self collision defaults off, and the consequence is visible: a strongly bent
    # body passes straight through itself, which a real animal obviously cannot do.
    # Turning it on is safe here, contrary to what the geometry first suggests.
    # Adjacent capsules do overlap by design -- a segment is 4.17 mm long and up to
    # 3.25 mm in radius -- but PhysX filters collisions between links joined by a
    # joint, which is exactly that set. The nearest unfiltered pair is i and i+2, at
    # 8.33 mm centre to centre against a radius sum of at most 6.47 mm, so they are
    # clear at rest and only meet when the body genuinely folds onto itself.
    #
    # It is off by default because it adds contact forces alongside the drag model
    # that represents the substrate (see build_scene), and every W2 locomotion
    # number was measured without it.
    physx_root = PhysxSchema.PhysxArticulationAPI.Apply(root.GetPrim())
    physx_root.CreateEnabledSelfCollisionsAttr().Set(self_collision)

    link_paths: list[str] = []
    for segment in segments:
        path = f"{root_path}/segment_{segment.index:02d}"
        capsule = UsdGeom.Capsule.Define(stage, path)
        # The capsule's own axis is X, matching the body axis.
        capsule.CreateAxisAttr().Set("X")
        capsule.CreateRadiusAttr().Set(float(segment.radius_m))
        # Height excludes the hemispherical caps, so shorten it to keep the
        # segment's total extent equal to its nominal length.
        capsule.CreateHeightAttr().Set(float(max(segment.length_m - 2 * segment.radius_m, 1e-6)))
        capsule.AddTranslateOp().Set(Gf.Vec3d(float(segment.centre_x_m), 0.0, float(z)))

        prim = capsule.GetPrim()
        UsdPhysics.CollisionAPI.Apply(prim)
        UsdPhysics.RigidBodyAPI.Apply(prim)
        mass_api = UsdPhysics.MassAPI.Apply(prim)
        mass_api.CreateMassAttr().Set(float(masses[segment.index]))
        link_paths.append(path)

    joint_paths: list[str] = []
    for i in range(plan.n_joints):
        parent, child = segments[i], segments[i + 1]
        path = f"{root_path}/joint_{i:02d}"
        joint = UsdPhysics.RevoluteJoint.Define(stage, path)
        joint.GetBody0Rel().SetTargets([Sdf.Path(link_paths[i])])
        joint.GetBody1Rel().SetTargets([Sdf.Path(link_paths[i + 1])])
        # Rotating about Z bends the body within the XY plane.
        joint.GetAxisAttr().Set("Z")
        # The joint sits on the boundary between the two capsules, expressed in
        # each body's own frame.
        half = parent.length_m / 2.0
        joint.GetLocalPos0Attr().Set(Gf.Vec3f(float(half), 0.0, 0.0))
        joint.GetLocalPos1Attr().Set(Gf.Vec3f(float(-child.length_m / 2.0), 0.0, 0.0))
        joint.GetLowerLimitAttr().Set(-float(joint_limit_deg))
        joint.GetUpperLimitAttr().Set(float(joint_limit_deg))

        # An explicit drive with zero gains: torque comes from the muscle model,
        # and PhysX must not add a position controller of its own.
        drive = UsdPhysics.DriveAPI.Apply(joint.GetPrim(), "angular")
        drive.CreateTypeAttr().Set("force")
        drive.CreateStiffnessAttr().Set(0.0)
        drive.CreateDampingAttr().Set(0.0)
        drive.CreateTargetPositionAttr().Set(0.0)
        joint_paths.append(path)

    if planar:
        _lock_to_plane(stage, root_path, link_paths[0])

    return root_path, joint_paths


def _lock_to_plane(stage: Usd.Stage, root_path: str, first_link: str) -> None:
    """Constrain the body to the XY plane with a world-anchored D6 joint.

    **Without this the model silently stops being the model.** The body plan is
    planar -- every revolute joint turns about Z, every drag force is in-plane --
    but the articulation has a free 6-DOF floating base, and nothing in the scene
    opposes rotation out of the plane once gravity is off. The explicit,
    per-segment drag forces are not perfectly symmetric, and that asymmetry is
    enough: the body slowly tips and then stands up vertically.

    It is a nasty failure because the body itself stays perfectly healthy. Its
    joint angles remain a clean travelling wave, its true 3-D end-to-end extent
    holds at 0.90 the whole time. Only the *projection* into XY collapses, which
    looks exactly like the body coiling into a ball -- and because
    :class:`~worm.body.drag.GroundDrag` builds its tangents from XY positions
    alone, the drag then becomes meaningless and drives the tip further.

    Locking translation in Z and rotation about X and Y leaves exactly the three
    degrees of freedom a worm crawling on a surface has: slide in X, slide in Y,
    turn about Z.
    """
    from pxr import Sdf, UsdPhysics

    joint = UsdPhysics.Joint.Define(stage, f"{root_path}/planar_constraint")
    joint.GetBody1Rel().SetTargets([Sdf.Path(first_link)])
    joint.CreateExcludeFromArticulationAttr().Set(True)

    # In UsdPhysics a limit whose low exceeds its high locks that axis, and an
    # axis with no LimitAPI at all is free. So naming only these three leaves
    # transX, transY and rotZ untouched.
    for axis in ("transZ", "rotX", "rotY"):
        limit = UsdPhysics.LimitAPI.Apply(joint.GetPrim(), axis)
        limit.CreateLowAttr().Set(1.0)
        limit.CreateHighAttr().Set(-1.0)


def build_scene(
    plan: BodyPlan | None = None,
    *,
    root_path: str = DEFAULT_ROOT,
    gravity: bool = False,
    ground: bool = False,
    planar: bool = True,
    self_collision: bool = False,
) -> tuple[str, list[str]]:
    """Create a fresh stage containing a physics scene, a light and the worm.

    Must be called after ``SimulationApp`` has been constructed.

    **Gravity and ground contact are off by default, and that is a modelling
    decision rather than a shortcut.** In this model the substrate is represented
    entirely by :class:`~worm.body.drag.GroundDrag`: an anisotropic resistive
    force standing in for a worm lying in a thin film of water on agar. Adding a
    real ground plane on top of that means the body is resisted twice, once by a
    drag model tuned to represent the surface and once by isotropic PhysX contact
    friction with a surface, and the two fight.

    Worse, the isotropic half wins in the sense that matters: it cannot tell
    sliding along the body from sliding across it, which is the exact distinction
    that makes undulation produce forward motion. With both enabled the body
    slides and bounces roughly as far with anisotropic drag as without it, so the
    anisotropy -- the whole mechanism -- stops being measurable.

    A planar model with no gravity is also what the standard neuromechanical
    models do (Boyle, Berri and Cohen 2012). The body stays in the XY plane on its
    own: every joint rotates about Z, every drag force is in-plane, and with no
    gravity there is nothing to push it out.

    Set ``gravity=True, ground=True`` to reproduce the contact-based version and
    see the anisotropy disappear.
    """
    import isaacsim.core.experimental.utils.stage as stage_utils
    from isaacsim.core.experimental.objects import DistantLight, GroundPlane
    from pxr import Gf, UsdPhysics

    stage_utils.create_new_stage()
    stage = stage_utils.get_current_stage()

    scene = UsdPhysics.Scene.Define(stage, "/World/PhysicsScene")
    scene.CreateGravityDirectionAttr().Set(Gf.Vec3f(0.0, 0.0, -1.0))
    scene.CreateGravityMagnitudeAttr().Set(9.81 if gravity else 0.0)

    if ground:
        GroundPlane("/World/GroundPlane", positions=[0, 0, 0])
    light = DistantLight("/World/DistantLight")
    light.set_intensities(2000)

    # Without gravity there is nothing to rest on, so the body sits in the z = 0
    # plane rather than being lifted clear of a ground plane.
    z_offset = None if ground else 0.0
    return build_worm(
        stage,
        plan,
        root_path=root_path,
        z_offset_m=z_offset,
        planar=planar,
        self_collision=self_collision,
    )


#: How many body lengths the camera frames across the view.
DEFAULT_CAMERA_FRAMING = 2.5


def add_ground_grid(
    plan: BodyPlan | None = None,
    *,
    path: str = "/World/GroundGrid",
    extent_body_lengths: float = 24.0,
    checker_body_lengths: float = 0.25,
) -> str:
    """A checkerboard under the worm, purely so movement is visible.

    **This has no collider and no rigid body, and that is deliberate rather than
    lazy.** In this model the substrate is represented entirely by
    :class:`~worm.body.drag.GroundDrag`, an anisotropic resistive force standing in
    for a worm lying in a film of water on agar. A real PhysX ground plane would
    resist the body a second time with *isotropic* contact friction, which cannot
    tell sliding along the body from sliding across it -- and that distinction is
    the entire mechanism by which undulation becomes forward motion. With contact
    friction added the anisotropy stops being measurable (see :func:`build_scene`).

    So this is scenery. It is drawn below the body, it is never touched by the
    solver, and switching it on or off cannot change a single number in a run.

    The checkers give the eye a fixed reference. Without one, a worm undulating in
    an empty void looks much the same whether it is travelling or thrashing in
    place -- which is exactly the distinction the headline results turn on.
    """
    import isaacsim.core.experimental.utils.stage as stage_utils
    from pxr import Gf, UsdGeom, Vt

    plan = plan or BodyPlan()
    stage = stage_utils.get_current_stage()

    half = 0.5 * extent_body_lengths * plan.total_length_m
    cell = checker_body_lengths * plan.total_length_m
    n = max(2, int(round(2.0 * half / cell)))
    # Sit just below the deepest part of the body so nothing z-fights with it.
    z = -1.05 * plan.max_radius_m

    step = 2.0 * half / n
    points = [
        Gf.Vec3f(float(-half + i * step), float(-half + j * step), float(z))
        for j in range(n + 1)
        for i in range(n + 1)
    ]

    counts: list[int] = []
    indices: list[int] = []
    colours: list[Gf.Vec3f] = []
    light, dark = Gf.Vec3f(0.32, 0.34, 0.38), Gf.Vec3f(0.20, 0.21, 0.24)
    for j in range(n):
        for i in range(n):
            base = j * (n + 1) + i
            # Counter-clockwise seen from +Z, which is where the camera is.
            indices += [base, base + 1, base + n + 2, base + n + 1]
            counts.append(4)
            colours.append(light if (i + j) % 2 == 0 else dark)

    mesh = UsdGeom.Mesh.Define(stage, path)
    mesh.CreatePointsAttr().Set(Vt.Vec3fArray(points))
    mesh.CreateFaceVertexCountsAttr().Set(Vt.IntArray(counts))
    mesh.CreateFaceVertexIndicesAttr().Set(Vt.IntArray(indices))
    # One colour per face, so the checkers need no texture file and no material --
    # nothing binary enters the repository for a piece of scenery.
    colour_attr = mesh.CreateDisplayColorAttr()
    colour_attr.Set(Vt.Vec3fArray(colours))
    colour_attr.SetMetadata("interpolation", UsdGeom.Tokens.uniform)
    mesh.CreateDoubleSidedAttr().Set(True)
    mesh.CreateExtentAttr().Set(Vt.Vec3fArray([Gf.Vec3f(-half, -half, z), Gf.Vec3f(half, half, z)]))
    # No CollisionAPI, no RigidBodyAPI, no physics schema of any kind. Anything
    # added here stops this being scenery and starts it being an experiment.
    return path


@dataclass
class WormCamera:
    """A top-down camera that keeps the worm in frame.

    Isaac's default perspective camera is placed for a scene measured in metres --
    several metres out, with its near clipping plane at 1 m. This body is 0.1 m
    long, so the default view is roughly fifty body lengths away and the animal is
    a speck near the origin. That is a property of the default viewport, not of the
    model, but it makes every visual check unnecessarily hard.

    The camera also *follows*, because a worm that is working travels several body
    lengths and would otherwise leave a fixed frame within seconds.
    """

    path: str
    height_m: float
    _translate: Any

    def follow(self, xy: Any) -> None:
        """Centre the view on ``xy``, typically the body centroid."""
        from pxr import Gf

        self._translate.Set(Gf.Vec3d(float(xy[0]), float(xy[1]), self.height_m))


def add_camera(
    plan: BodyPlan | None = None,
    *,
    path: str = "/World/WormCamera",
    framing: float = DEFAULT_CAMERA_FRAMING,
    activate: bool = True,
) -> WormCamera:
    """Add a top-down camera framed on the worm and make it the active view.

    Call after :func:`build_scene`. Safe to call headless, where there is no
    viewport to activate and the camera is simply authored into the stage.
    """
    import isaacsim.core.experimental.utils.stage as stage_utils
    from pxr import Gf, Sdf, UsdGeom

    plan = plan or BodyPlan()
    stage = stage_utils.get_current_stage()

    camera = UsdGeom.Camera.Define(stage, path)
    # A USD camera looks down its own -Z, so an unrotated camera above the plane
    # looks straight down at it. The body is planar, so that is the view that
    # shows the actual shape rather than a foreshortening of it.
    focal_length, aperture = 18.0, 20.955
    width_m = framing * plan.total_length_m
    height_m = width_m * focal_length / aperture

    camera.CreateFocalLengthAttr().Set(focal_length)
    camera.CreateHorizontalApertureAttr().Set(aperture)
    # The default near plane of 1 m sits far behind a body 0.1 m across.
    camera.CreateClippingRangeAttr().Set(Gf.Vec2f(0.001, 100.0))
    translate = camera.AddTranslateOp()
    translate.Set(Gf.Vec3d(0.0, 0.0, height_m))

    if activate:
        try:
            from omni.kit.viewport.utility import get_active_viewport

            viewport = get_active_viewport()
        except (ImportError, AttributeError):  # headless: no viewport to point
            viewport = None
        if viewport is not None:
            viewport.set_active_camera(Sdf.Path(path))

    return WormCamera(path=path, height_m=height_m, _translate=translate)


def dof_order(articulation: Any, joint_paths: list[str]) -> list[int]:
    """Map our head-to-tail joint order onto PhysX's internal DOF order.

    **This is not optional, and getting it wrong is nearly invisible.**

    PhysX does not index an articulation's degrees of freedom in the order the
    joints were authored. It picks its own root link -- for a serial chain it
    takes the middle, to keep the tree balanced -- and numbers outward from
    there. For a 24-segment worm the DOF order comes out as::

        joint_10, joint_11, joint_09, joint_12, joint_08, joint_13, ...

    so writing a head-to-tail torque array straight into ``set_dof_efforts``
    applies the head's command to a mid-body joint. The result still has exactly
    the right bend amplitude and the right distribution of angles -- every summary
    statistic matches what the model intended -- while the body's actual shape is
    scrambled. It coils instead of undulating, and no amount of retuning the
    torque fixes it, because nothing about the magnitudes is wrong.

    Returns indices to pass as ``dof_indices`` when reading or writing DOF state,
    so that array position 0 always means the head-most joint.
    """
    names = [path.rsplit("/", 1)[-1] for path in joint_paths]
    import numpy as np

    return np.asarray(articulation.get_dof_indices(names)).reshape(-1).tolist()
