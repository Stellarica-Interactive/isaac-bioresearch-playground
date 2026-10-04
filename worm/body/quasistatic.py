"""Quasi-static body mechanics: solve the force balance, do not integrate momentum.

    body = QuasiStaticBody(plan)
    body.step(joint_torques, dt_s=1/60)
    positions = body.segment_centres()      # for rendering and for touch
    angles = body.joint_angles              # for proprioception

Why this exists
---------------

A worm on agar is deeply overdamped. Measured on this body, a segment's velocity
relaxes in ``m/c = 2.35 ms`` while the committed physics step is ``4.17 ms``, so
drag erases velocity faster than the simulator advances. Integrating momentum in
that regime makes the result depend on the step rather than on the mechanics:
the scripted travelling wave covers **7.216 BL at 240 Hz and 1.442 BL at 60 Hz**
with identical amplitude and phase. Same gait, five times less travel. See
``docs/model_assumptions.md`` §5Y and §5AE.

Making that solver correct means ``c dt/m`` well below one, which needs about
4250 Hz -- roughly forty-five minutes of wall clock per ten simulated seconds.
Technically valid, practically unusable.

So inertia is dropped, which is what the medium already does physically and what
the standard neuromechanical models do (Boyle, Berri & Cohen 2012; Sibernetic).
At every instant the body moves at whatever velocity makes the external forces
balance. There is no stability constraint from inertia at all, so the timestep is
limited only by how fast the *shape* changes.

The mathematics
---------------

The configuration is the first node's position, the body's heading, and the joint
angles: ``(x, y, theta, q)``. Every segment's velocity is linear in the
configuration rates, so with the shape rates ``qdot`` known, the three rigid
degrees of freedom follow from requiring no net external force or torque --
nothing pushes the animal but the ground it is pushing against::

    sum_i  F_i = 0           F_i = -D_i v_i
    sum_i  r_i x F_i = 0     D_i = c_par t_i t_i^T + c_perp n_i n_i^T

which is a 3x3 linear system. ``D_i`` is the anisotropic drag of
:mod:`worm.body.drag` written as a tensor; the anisotropy is what converts a
travelling shape change into forward motion, and with ``c_par == c_perp`` the
solution is exactly zero net displacement, which the tests assert.

What is simplified, and what that costs
---------------------------------------

**The joints were closed separately in the first version, and should not have
been.** Their rates came from their own damping alone, ``qdot = (torque - k q)/b``,
leaving the medium's resistance to *bending* out. The shape then changed too
freely, and the zero-net-torque condition answered the large shape rates with a
large counter-rotation of the whole body -- visible in the viewport as a worm
yawing back and forth rather than travelling. Every configuration rate is now an
unknown of the same balance, so the medium resists bending through the same
``J^T D J`` that resists translation.

**No contact, no self-collision.** This solver knows only drag. The existing
touch model is already proximity-based rather than contact-based -- the probe
deliberately carries no collider (§5D) -- so touch is unaffected, but a body that
should not pass through itself needs a separate check.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from worm.body.drag import DragParameters
from worm.body.geometry import BodyPlan
from worm.body.muscles import MuscleParameters
from worm.isaac_limits import JOINT_LIMIT_RAD

#: How close to a joint limit counts as being at it, radians.
#:
#: ``JOINT_LIMIT_RAD`` is the same 60 degree ceiling the Isaac articulation
#: enforces, imported rather than restated so the two solvers constrain the body
#: identically -- one free to fold further would not be comparable with the
#: committed runs. This tolerance is one part in ten thousand of it: far below
#: any angle that matters biologically, and far above the rounding of a single
#: integration step.
_LIMIT_TOL = 1.0e-4 * JOINT_LIMIT_RAD

#: Fixed-point passes allowed when the drag law is nonlinear, and the relative
#: change below which the iteration stops.
#:
#: Measured over sixty steps of the committed gait with both nonlinearities on,
#: as net unbalanced force over total drag force:
#:
#: ===========  ======  ==============  ======
#: tolerance    passes  residual        cost
#: ===========  ======  ==============  ======
#: 1e-6         24      5.0e-5          0.70 s
#: **1e-10**    **80**  **5.3e-9**      1.27 s
#: 1e-14        250     4.8e-14         27.3 s
#: ===========  ======  ==============  ======
#:
#: 1e-10 costs 1.8x the loosest setting and buys four orders of magnitude, so it
#: is the default. 1e-14 is below what the iteration reaches quickly and burns
#: the whole budget for no physical gain. See model_assumptions 5AH.4.
_NONLINEAR_PASSES = 80
_NONLINEAR_TOL = 1.0e-10

#: Band around touching within which self-contact resists approach, metres.
#:
#: ARBITRARY ENGINEERING, but bounded above by a measurement rather than chosen
#: freely: **the scripted gait's closest self-approach is 1.57 mm**, so a band
#: wider than that would engage during ordinary crawling and make every
#: locomotion number contingent on this contact model. 1 mm leaves a margin and,
#: measured, changes the gait's speed by 0.00 per cent -- 0.2176 BL/s with
#: contact on and off alike.
#:
#: Bounded below by the timestep: the band has to be wide enough that an
#: approaching segment is seen inside it before it passes through.
CONTACT_BAND_M = 1.0e-3

#: Drag coefficient along a contact normal, as a multiple of the perpendicular
#: body drag. Large enough that approach stops, small enough to leave the solve
#: well conditioned.
#:
#: ARBITRARY ENGINEERING, and the reason this is a penalty rather than a
#: constraint: a true non-penetration constraint would add a row per contact to
#: the active set of ``_solve_with_limits`` and iterate over contacts and joint
#: limits together. This is the smaller change and it behaves like the groove of
#: model_assumptions 5AH -- very stiff resistance to one direction of motion.
#:
#: Measured together with the band above, on a curl of the depth the connectome
#: actually produces: worst self-clearance goes from **-5.80 mm to +0.79 mm**,
#: so the body stops folding through itself. It cannot undo every overlap -- at
#: the 60 degree joint limit the body wraps 3.8 times and self-intersection is
#: geometrically forced, since a closed circle of this body needs only 15.7
#: degrees per joint. That is a fact about the joint limit, not about contact.
CONTACT_STIFFNESS = 2.0e4

#: How much of the body's length of track to remember, in body lengths. The
#: groove only matters where the body still lies in it, so remembering much more
#: than one body length is wasted work; remembering less leaves the tail without
#: a groove to sit in.
#:
#: ARBITRARY ENGINEERING, bounded by geometry on both sides.
GROOVE_MEMORY_BODY_LENGTHS = 1.5

#: Extra drag across the groove, as a multiple of the perpendicular body drag.
#:
#: ASSUMED. Rabets et al. 2014 established that groove formation dominates agar's
#: response; they published no stiffness for it, and this project has no
#: measurement to set one from. It is swept rather than fitted and the sweep is
#: the result -- see model_assumptions 5AL.
GROOVE_STIFFNESS = 4.0

#: Spacing of remembered track points, in segment lengths. Finer costs more
#: distance computations for a track that is already smooth at this scale.
GROOVE_SAMPLE_SEGMENTS = 0.5

#: How near a remembered track point a segment must be to be *in* the groove, as
#: a multiple of its own radius.
#:
#: Not a tuning parameter but a missing condition. The first version applied the
#: groove to every segment regardless of distance, and measured, the nearest
#: track point was **25 to 45 mm away** on a 100 mm body, with the groove's
#: direction 25 to 34 degrees off the segment's own (worst case 69). So every
#: segment was being resisted along a direction taken from somewhere it had never
#: been. The body stopped -- 0.002 BL/s at the stiffest setting -- and that
#: deadlocked the mechanism: the groove forms where the body has been, the body
#: only gets anywhere if the groove does not resist it, and a groove applied
#: before a track exists resists everything forever.
#:
#: A groove only confines a body that is lying in it. Two radii is "touching the
#: wall of its own channel".
GROOVE_REACH_RADII = 2.0

#: Half-width of the channel, as a multiple of the segment's radius. Inside it
#: the segment is free; beyond it, and only while moving further out, the groove
#: constrains.
#:
#: A knife-edge groove does not work and the reason is instructive. Requiring
#: ``v . n = 0`` exactly, at every segment lying near the track, pins the body
#: completely -- measured 0.000 BL/s with only 10 or 11 active rows and 15 DOF
#: still free, because each row carries the translation columns and ten rows with
#: differing normals leave translation nowhere to go. It demands the body conform
#: exactly to the remembered track, and any mismatch between the body's shape and
#: the channel's kills all motion.
#:
#: A real groove is a shallow depression of finite width that the animal rides
#: within. So the constraint is unilateral with a dead zone -- the same structure
#: as the joint limits, which work (§5AF.6).
GROOVE_CHANNEL_RADII = 1.0


def _perp(v: np.ndarray) -> np.ndarray:
    """Rotate planar vectors 90 degrees: ``(x, y) -> (-y, x)``."""
    return np.stack([-v[..., 1], v[..., 0]], axis=-1)


@dataclass
class QuasiStaticBody:
    """A planar chain whose velocity solves the force balance each step."""

    plan: BodyPlan
    drag: DragParameters = field(default_factory=DragParameters)
    muscle: MuscleParameters = field(default_factory=MuscleParameters)

    #: How the groove is enforced, when :attr:`groove` is on. Three forms, and
    #: the first two are rejected in model_assumptions 5AL:
    #:
    #: ``"penalty"`` -- stiff drag across the *track*. Immobilises the animal,
    #: because the track's normal sits 25-34 degrees off the body's, so the term
    #: leaks 7 to 12 times the body's own tangential drag into the direction it
    #: is travelling.
    #:
    #: ``"constraint"`` -- forbid motion across the track exactly, by Lagrange
    #: multiplier. Also immobilises it, for the opposite reason: a hard condition
    #: at a dozen points demands the body conform exactly to its remembered
    #: track, and the 0.74 mm it never conforms by is enough for the rows to
    #: conflict and the solution to collapse.
    #:
    #: ``"restoring"`` -- a force toward the channel centreline along the
    #: **body's own normal**, which resistive force theory already resists, so it
    #: cannot leak into the tangent. Zero inside the channel. This is the form the
    #: other two diagnose rather than a new guess.
    groove_mode: str = "restoring"

    #: Whether the medium remembers the track the body has cut, and resists
    #: motion across it. Off by default. See :data:`GROOVE_STIFFNESS`.
    groove: bool = False

    #: Whether the body resists passing through itself. Off by default: every
    #: number in model_assumptions was measured without it, and the connectome
    #: run is the only configuration that ever reaches contact (`clear` -1.39 mm,
    #: against +1.57 and better for the scripted gait).
    self_contact: bool = False

    #: First node position, body heading, and joint angles.
    origin: np.ndarray = field(init=False)
    heading: float = field(init=False, default=0.0)
    joint_angles: np.ndarray = field(init=False)

    #: Segment radii, cached: the contact test needs them every step.
    _radii: np.ndarray = field(init=False, repr=False)

    #: Remembered track: positions the head has occupied, oldest first, and the
    #: unit direction it was travelling at each. The groove's geometry.
    _track: np.ndarray = field(init=False, repr=False)
    _track_dirs: np.ndarray = field(init=False, repr=False)

    #: Groove drag tensors, recomputed once per step rather than once per pass of
    #: the nonlinear fixed point: they depend on positions, which do not change
    #: within a step.
    _groove_tensors: np.ndarray | None = field(init=False, default=None, repr=False)

    #: Last solved segment velocities, used to warm-start the nonlinear fixed
    #: point. Not part of the configuration -- a quasi-static body has no
    #: velocity state, and this is purely the previous iteration's answer kept as
    #: the next one's starting guess. The converged solution does not depend on
    #: it, only the number of passes needed to reach it.
    _last_velocities: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.origin = np.zeros(2, dtype=np.float64)
        self.joint_angles = np.zeros(self.plan.n_joints, dtype=np.float64)
        self._last_velocities = np.zeros((self.plan.n_segments, 2), dtype=np.float64)
        self._radii = np.array(
            [self.plan.radius_at(i) for i in range(self.plan.n_segments)],
            dtype=np.float64,
        )
        self._track = np.zeros((0, 2), dtype=np.float64)
        self._track_dirs = np.zeros((0, 2), dtype=np.float64)

    # -- kinematics --------------------------------------------------------

    @property
    def segment_length_m(self) -> float:
        return self.plan.total_length_m / self.plan.n_segments

    def segment_angles(self) -> np.ndarray:
        """Direction of each segment as an angle, radians.

        Segment 0 carries the body heading and each joint adds its own, so this
        is ``heading`` only for the head. Anything that needs a particular
        segment's orientation -- such as driving an articulation whose root link
        is not the head -- has to read it from here rather than using
        ``heading``.
        """
        return self.heading + np.concatenate([[0.0], np.cumsum(self.joint_angles)])

    def _directions(self) -> np.ndarray:
        """Unit heading of each segment, from the chain of joint angles."""
        angles = self.segment_angles()
        return np.stack([np.cos(angles), np.sin(angles)], axis=1)

    def nodes(self) -> np.ndarray:
        """The ``n_segments + 1`` joint positions along the body."""
        step = self.segment_length_m * self._directions()
        return np.vstack([self.origin, self.origin + np.cumsum(step, axis=0)])

    def segment_centres(self) -> np.ndarray:
        node = self.nodes()
        return 0.5 * (node[:-1] + node[1:])

    # -- dynamics ----------------------------------------------------------

    def _jacobian(self) -> np.ndarray:
        """Segment velocities per configuration rate, shape ``(n, 2, 3 + joints)``.

        Columns are: translation in x, translation in y, rotation about the first
        node, then one per joint. A joint swings everything downstream of it
        about its own node, which is why its column is ``perp(centre - node)``.
        """
        length = self.segment_length_m
        directions = self._directions()
        nodes = self.nodes()
        centres = self.segment_centres()
        n, joints = self.plan.n_segments, self.plan.n_joints

        jac = np.zeros((n, 2, 3 + joints), dtype=np.float64)
        jac[:, 0, 0] = 1.0
        jac[:, 1, 1] = 1.0
        jac[:, :, 2] = _perp(centres - self.origin)
        # Joint j moves segment i only if it lies upstream of it.
        for j in range(joints):
            downstream = np.arange(n) > j
            jac[downstream, :, 3 + j] = _perp(centres[downstream] - nodes[j + 1])
        del length, directions
        return jac

    def _drag_tensors(self, velocities: np.ndarray | None = None) -> np.ndarray:
        """Anisotropic drag per segment as a 2x2 tensor.

        ``c_par t t^T + c_perp n n^T`` -- the same coefficients as
        :mod:`worm.body.drag`, written so they can enter a linear solve. The
        difference between the two coefficients is the entire reason undulation
        becomes thrust; with them equal the solution is no net travel.

        With a sublinear exponent the coefficients depend on the component speeds
        and ``velocities`` must be supplied; each direction is scaled by its own
        speed, because the groove the exponent stands in for confines the body
        laterally without resisting it sliding along its own length.
        """
        directions = self._directions()
        normals = _perp(directions)
        along = np.einsum("ij,ik->ijk", directions, directions)
        across = np.einsum("ij,ik->ijk", normals, normals)

        c_par, c_perp = self.drag.tangential, self.drag.perpendicular
        if not self.drag.is_linear:
            if velocities is None:
                raise ValueError("a nonlinear drag law needs the segment velocities")
            speed_par = np.abs(np.einsum("ij,ij->i", velocities, directions))
            speed_perp = np.abs(np.einsum("ij,ij->i", velocities, normals))
            # The power law rescales the viscous coefficients; the yield term is
            # added on top, because a threshold is a separate physical claim
            # from a nonlinear slope and the two are swept independently.
            c_par = c_par * self.drag.speed_factor(speed_par) + self.drag.yield_coefficient(
                speed_par
            )
            c_perp = c_perp * self.drag.speed_factor(speed_perp) + self.drag.yield_coefficient(
                speed_perp
            )
            tensors = c_par[:, None, None] * along + c_perp[:, None, None] * across
        else:
            tensors = c_par * along + c_perp * across
        if self.self_contact:
            extra = self._contact_drag(velocities)
            if extra is not None:
                tensors = tensors + extra
        if self.groove and self._groove_tensors is not None:
            tensors = tensors + self._groove_tensors
        return tensors

    def _remember_track(self) -> None:
        """Record where the head is, if it has moved far enough to be new.

        The head cuts the groove, so the head's path *is* the track. Recorded by
        distance rather than by time, so a body that has stopped does not fill the
        buffer with one repeated point and does not forget the track it is in.
        """
        nose = self.nodes()[0]
        spacing = GROOVE_SAMPLE_SEGMENTS * self.segment_length_m
        if len(self._track) and np.linalg.norm(nose - self._track[-1]) < spacing:
            return

        # The direction the head *travelled*, not the direction it was pointing.
        # Those differ by up to 48 degrees, because the head's own tangent sweeps
        # that far every cycle while its path does not: a groove is cut along the
        # path. Using the heading instead gave the groove the head's undulation,
        # so every segment was resisted across a direction oscillating wildly
        # about the one it was actually moving along, and the body stopped.
        # Same confusion as 5AE -- the head is not the body.
        if len(self._track):
            step = nose - self._track[-1]
            norm = float(np.linalg.norm(step))
            direction = step / norm if norm > 1e-12 else self._directions()[0]
        else:
            direction = self._directions()[0]
        self._track = np.vstack([self._track, nose])
        self._track_dirs = np.vstack([self._track_dirs, direction])
        keep = int(GROOVE_MEMORY_BODY_LENGTHS * self.plan.total_length_m / spacing)
        if len(self._track) > keep:
            self._track = self._track[-keep:]
            self._track_dirs = self._track_dirs[-keep:]

    def _groove_drag(self) -> np.ndarray | None:
        """Extra drag across the track, per segment, or ``None`` if there is none.

        For each segment the nearest remembered track point gives the groove's
        local direction, and the drag is added along the normal to it. Where the
        body lies in its own track that normal is the body's own normal and this
        duplicates the existing anisotropy; where the body has swept sideways the
        two differ, and the difference is what resists.
        """
        if len(self._track) < 2:
            return None
        centres = self.segment_centres()
        # Nearest remembered point for each segment: one 24-by-N distance matrix.
        offsets = centres[:, None, :] - self._track[None, :, :]
        square = np.einsum("nmi,nmi->nm", offsets, offsets)
        nearest = np.argmin(square, axis=1)
        distance = np.sqrt(square[np.arange(len(centres)), nearest])

        # A groove only confines a body lying in it. Without this the groove's
        # direction is taken from wherever the track happens to be nearest, which
        # measured 25 to 45 mm away, and the body is resisted along a direction
        # it has never travelled. See GROOVE_REACH_RADII.
        inside = distance < GROOVE_REACH_RADII * self._radii
        if not inside.any():
            return None

        normals = _perp(self._track_dirs[nearest])
        scale = GROOVE_STIFFNESS * self.drag.perpendicular * inside
        return scale[:, None, None] * np.einsum("ij,ik->ijk", normals, normals)

    def _restoring_drag(self) -> np.ndarray | None:
        """Lateral drag along the body's own normal, for segments outside their
        channel.

        The two rejected forms fail in opposite directions: a penalty across the
        *track* leaks into the body's tangent because the two normals disagree by
        25 to 34 degrees, and a hard constraint demands a conformance the body
        never has. This keeps the penalty's compliance and the constraint's
        honesty about direction: the force acts along the body's own normal --
        a direction the medium already resists, so forward motion is untouched --
        and its size grows with how far outside the channel the segment is.
        """
        if len(self._track) < 2:
            return None
        centres = self.segment_centres()
        offsets = centres[:, None, :] - self._track[None, :, :]
        square = np.einsum("nmi,nmi->nm", offsets, offsets)
        nearest = np.argmin(square, axis=1)
        distance = np.sqrt(square[np.arange(len(centres)), nearest])

        track_normals = _perp(self._track_dirs[nearest])
        across = np.abs(np.einsum("na,na->n", centres - self._track[nearest], track_normals))
        width = GROOVE_CHANNEL_RADII * self._radii
        # How far outside the channel, in radii. Zero inside it.
        excess = np.maximum(0.0, across - width) / self._radii
        engaged = (distance < GROOVE_REACH_RADII * self._radii) & (excess > 0.0)
        if not engaged.any():
            return None

        body_normals = _perp(self._directions())
        scale = GROOVE_STIFFNESS * self.drag.perpendicular * excess * engaged
        return scale[:, None, None] * np.einsum("ij,ik->ijk", body_normals, body_normals)

    def _contact_drag(self, velocities: np.ndarray | None) -> np.ndarray | None:
        """Extra per-segment drag resisting self-approach, or ``None`` if clear.

        A contact is really a constraint -- two segments that have closed to
        touching must not approach further -- and the honest treatment is an
        extra row in the active set of :meth:`_solve_with_limits`. This is a
        penalty instead: along the line between two touching segments, and only
        when they are closing, the medium becomes very stiff. Same shape as the
        yield force of :mod:`worm.body.drag`, and the same physics as the groove
        it stands in for.

        Only non-neighbouring pairs count. Adjacent segments share a node and
        always overlap by construction, so including them would make a straight
        body self-collide.
        """
        centres = self.segment_centres()
        radii = self._radii
        delta = centres[:, None, :] - centres[None, :, :]
        distance = np.linalg.norm(delta, axis=-1)
        gap = distance - (radii[:, None] + radii[None, :])
        index = np.arange(self.plan.n_segments)
        far_enough = np.abs(index[:, None] - index[None, :]) > 1
        close = far_enough & (gap < CONTACT_BAND_M) & (distance > 1e-12)
        pairs = np.argwhere(np.triu(close))
        if pairs.size == 0:
            return None

        extra = np.zeros((self.plan.n_segments, 2, 2), dtype=np.float64)
        scale = CONTACT_STIFFNESS * self.drag.perpendicular
        for i, j in pairs:
            normal = delta[i, j] / distance[i, j]
            if velocities is not None:
                # Only on approach: a contact resists closing, never separation,
                # and a penalty that resisted both would glue the body together.
                closing = float((velocities[j] - velocities[i]) @ normal)
                if closing > 0.0:
                    continue
            outer = scale * np.outer(normal, normal)
            extra[i] += outer
            extra[j] += outer
        return extra

    def _groove_constraints(self, jac: np.ndarray) -> np.ndarray | None:
        """Rows forbidding motion across the track, for segments lying in it.

        One row per in-groove segment: ``n^T J_i u = 0``, where ``n`` is the
        normal to the track at that segment's nearest remembered point. The force
        that enforces it is a Lagrange multiplier and comes out of the solve, so
        there is no stiffness to choose -- which is the difference from the
        penalty of §5AL, where the force was a constant somebody picked and ended
        up resisting the direction the body was travelling.
        """
        if len(self._track) < 2:
            return None
        centres = self.segment_centres()
        offsets = centres[:, None, :] - self._track[None, :, :]
        square = np.einsum("nmi,nmi->nm", offsets, offsets)
        nearest = np.argmin(square, axis=1)
        distance = np.sqrt(square[np.arange(len(centres)), nearest])
        near = distance < GROOVE_REACH_RADII * self._radii
        normals = _perp(self._track_dirs[nearest])
        # Signed offset across the channel, and the channel's half-width.
        offset = np.einsum("na,na->n", centres - self._track[nearest], normals)
        beyond = np.abs(offset) > GROOVE_CHANNEL_RADII * self._radii
        which = np.flatnonzero(near & beyond)
        if which.size == 0:
            return None
        # Constrain only the direction that would take the segment further out,
        # so a segment riding back toward its channel is free -- unilateral, with
        # the channel as the dead zone, exactly like a joint limit.
        signed = normals[which] * np.sign(offset[which])[:, None]
        return np.einsum("ka,kaj->kj", signed, jac[which])

    def _solve_rates(self, joint_torques: np.ndarray) -> np.ndarray:
        """``(vx, vy, omega, qdot)`` from one generalised force balance.

        The three rigid rows are homogeneous: no external force or torque acts on
        the animal beyond the ground it is pushing against. The joint rows carry
        the muscle torque against the passive stiffness, and the medium resists
        bending through the same ``J^T D J`` that resists translation -- which is
        what the earlier two-stage version left out.

        Joints held at their limit are **locked out of the unknowns** and the
        balance re-solved over the rest, rather than solved freely and then
        clipped. Clipping is not a harmless safety net here: the rigid rates are
        solved jointly with ``qdot``, so discarding part of ``qdot`` leaves a
        velocity field that no longer satisfies zero net torque, and the leftover
        rotation accumulates every step without bound. A uniform drive of
        1e-3 N m saturates all 23 joints and spun the body at **-26814 deg/s,
        sustained**. The committed drive does not reach that, so this is not what
        anyone has seen on screen; it fires on any drive strong enough to pin a
        joint. A joint limit is an *internal* constraint, equal and
        opposite across the joint like the muscle torque, so it cannot torque the
        animal as a whole; locking the joint and re-solving is therefore exact,
        not an approximation. See model_assumptions 5AF.6.
        """
        jac = self._jacobian()
        joints = self.plan.n_joints
        forces = np.zeros(3 + joints, dtype=np.float64)
        forces[3:] = np.asarray(joint_torques, dtype=np.float64) - (
            self.muscle.joint_stiffness * self.joint_angles
        )

        rows = (
            self._groove_constraints(jac)
            if self.groove and self.groove_mode == "constraint"
            else None
        )

        if self.drag.is_linear:
            return self._solve_with_limits(jac, self._drag_tensors(), forces, rows)

        # A sublinear law makes the coefficients depend on the velocities they
        # determine, so the balance is nonlinear and is closed by a fixed point:
        # solve with the current coefficients, re-read the velocities, repeat.
        # A power-law coefficient is a contraction and converges in two or three
        # passes. The yield term varies as 1/|v| and is not, so the pass budget
        # is set for it; the tolerance ends the loop early in the common case.
        # The linear solution is the starting guess.
        # Warm start from the previous solve. The body moves slowly against a
        # physics step, so these coefficients are nearly right already and the
        # loop usually exits after one or two passes instead of running its
        # budget. The fixed point is the same either way; only its cost changes.
        rates = self._solve_with_limits(
            jac, self._drag_tensors(self._last_velocities), forces, rows
        )
        velocities = np.einsum("naj,j->na", jac, rates)
        for _ in range(_NONLINEAR_PASSES):
            previous = rates
            rates = self._solve_with_limits(jac, self._drag_tensors(velocities), forces, rows)
            velocities = np.einsum("naj,j->na", jac, rates)
            if np.max(np.abs(rates - previous)) <= _NONLINEAR_TOL * max(
                float(np.max(np.abs(rates))), 1e-30
            ):
                break
        self._last_velocities = velocities
        return rates

    def _solve_with_limits(
        self,
        jac: np.ndarray,
        drag: np.ndarray,
        forces: np.ndarray,
        constraints: np.ndarray | None = None,
    ) -> np.ndarray:
        """One force balance, with joints at their limits locked out of it.

        ``constraints`` are rows ``C`` requiring ``C u = 0`` -- a channel the body
        may move along and not out of. They enter as Lagrange multipliers rather
        than as forces, so the medium's reaction comes out of the solve instead of
        being a stiffness somebody chose. Solved by least squares, because up to
        24 such rows on 26 degrees of freedom can be redundant or momentarily
        inconsistent and a direct solve would raise rather than degrade.
        """
        joints = self.plan.n_joints
        # J^T D J, summed over segments.
        matrix = np.einsum("nai,nab,nbj->ij", jac, drag, jac)
        # The joints' own internal damping, which is not part of the medium.
        matrix[3:, 3:] += np.eye(joints) * self.muscle.joint_damping

        # Active-set iteration. Each pass solves over the free joints, then locks
        # any that the solution drives further into a limit they already sit on.
        # Locking only ever shrinks the free set, so this terminates in at most
        # one pass per joint; in practice it converges in one or two.
        locked = np.zeros(joints, dtype=bool)
        at_high = self.joint_angles >= JOINT_LIMIT_RAD - _LIMIT_TOL
        at_low = self.joint_angles <= -JOINT_LIMIT_RAD + _LIMIT_TOL
        rates = np.zeros(3 + joints, dtype=np.float64)
        for _ in range(joints + 1):
            free = np.flatnonzero(~locked)
            index = np.concatenate([[0, 1, 2], 3 + free])
            sub = matrix[np.ix_(index, index)]
            if constraints is None or constraints.shape[0] == 0:
                solved = np.linalg.solve(sub, forces[index])
            else:
                rows = constraints[:, index]
                width = len(index)
                kkt = np.zeros((width + rows.shape[0], width + rows.shape[0]))
                kkt[:width, :width] = sub
                kkt[:width, width:] = rows.T
                kkt[width:, :width] = rows
                rhs = np.zeros(width + rows.shape[0])
                rhs[:width] = forces[index]
                solved = np.linalg.lstsq(kkt, rhs, rcond=None)[0][:width]
            rates = np.zeros(3 + joints, dtype=np.float64)
            rates[:3] = solved[:3]
            rates[3 + free] = solved[3:]
            pushing_out = (at_high & (rates[3:] > 0.0)) | (at_low & (rates[3:] < 0.0))
            newly = pushing_out & ~locked
            if not newly.any():
                break
            locked |= newly
        return rates

    def step(self, joint_torques: np.ndarray, dt_s: float) -> None:
        """Advance the configuration by ``dt_s``.

        No inertial stability limit: the only thing the step has to resolve is
        the change in shape, so 60 Hz against a two-second gait is ample.
        """
        if self.groove:
            # Once per step: the groove depends on positions, which do not change
            # within a step, so recomputing it inside the fixed point would cost
            # 26 distance matrices for one answer.
            self._remember_track()
            if self.groove_mode == "penalty":
                self._groove_tensors = self._groove_drag()
            elif self.groove_mode == "restoring":
                self._groove_tensors = self._restoring_drag()
            else:
                self._groove_tensors = None

        rates = self._solve_rates(joint_torques)
        if not np.all(np.isfinite(rates)):
            raise FloatingPointError("force balance has no finite solution")

        self.origin = self.origin + dt_s * rates[:2]
        self.heading = float(self.heading + dt_s * rates[2])
        # The solve already holds saturated joints still, so this only trims the
        # fraction of a step by which a joint first reaches its limit. It must
        # stay: without it a joint lands just past the limit on the step that
        # reaches it and never registers as being *at* it.
        self.joint_angles = np.clip(
            self.joint_angles + dt_s * rates[3:], -JOINT_LIMIT_RAD, JOINT_LIMIT_RAD
        )

    # -- diagnostics -------------------------------------------------------

    def min_self_distance_m(self) -> float:
        """Closest approach between two non-neighbouring segments, minus their radii.

        Positive is clear, zero is touching, negative is interpenetrating. This
        solver has no self-collision, so nothing prevents a negative value; the
        point of measuring is to find out whether that permission is ever used.

        Geometry says it should rarely bind. The body is at most 3.25 mm in
        radius and a closed loop of it has a radius of 15.9 mm, so opposite sides
        of a coil clear each other by a factor of five, and even the tightest
        hairpin the 60 degree joint limit allows comes to about 4.2 mm. Only a
        shape that doubles back more than once -- a spiral rather than a coil --
        could intersect. See model_assumptions 5AF.7.

        Immediate neighbours are skipped: adjacent segments share a node and so
        always "overlap" by construction, which says nothing about the body
        folding onto itself.
        """
        centres = self.segment_centres()
        radii = self._radii
        gap = np.linalg.norm(centres[:, None, :] - centres[None, :, :], axis=-1) - (
            radii[:, None] + radii[None, :]
        )
        # Mask the diagonal and the two off-diagonals: self, and the neighbours
        # that share a node.
        index = np.arange(self.plan.n_segments)
        neighbours = np.abs(index[:, None] - index[None, :]) <= 1
        return float(gap[~neighbours].min())

    def residual_fraction(self, joint_torques: np.ndarray) -> float:
        """Net unbalanced force as a fraction of the total drag force.

        The absolute residual is the right test for the linear law, where a
        single direct solve balances to machine precision. It is the wrong test
        for a nonlinear one: the coefficients depend on the velocities they
        produce, so the balance is closed by a fixed point and carries that
        iteration's convergence error. Asserting an absolute bound against it
        just encodes the drag magnitude.

        This normalises by the summed magnitude of the per-segment forces, which
        is scale-free and says what it means -- how much of the force in the
        system fails to cancel. See model_assumptions 5AH.4.
        """
        rates = self._solve_rates(joint_torques)
        jac = self._jacobian()
        velocities = np.einsum("naj,j->na", jac, rates)
        forces = -np.einsum("nab,nb->na", self._drag_tensors(velocities), velocities)
        total = float(np.abs(forces).sum())
        if total == 0.0:
            return 0.0
        return float(np.linalg.norm(forces.sum(axis=0))) / total

    def residual_force_and_torque(self, joint_torques: np.ndarray) -> tuple[float, float]:
        """How well the solved motion balances, for tests to assert on.

        Zero to solver precision. A non-zero residual means the body is being
        pushed by something that is not the ground -- the failure that made a
        frozen body appear to rotate for an hour (model_assumptions 5AE).
        """
        rates = self._solve_rates(joint_torques)
        jac = self._jacobian()
        velocities = np.einsum("naj,j->na", jac, rates)
        drag = self._drag_tensors(velocities)
        forces = -np.einsum("nab,nb->na", drag, velocities)
        arms = self.segment_centres() - self.origin
        net_force = float(np.linalg.norm(forces.sum(axis=0)))
        net_torque = float(np.sum(arms[:, 0] * forces[:, 1] - arms[:, 1] * forces[:, 0]))
        return net_force, net_torque
