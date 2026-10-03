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


def _perp(v: np.ndarray) -> np.ndarray:
    """Rotate planar vectors 90 degrees: ``(x, y) -> (-y, x)``."""
    return np.stack([-v[..., 1], v[..., 0]], axis=-1)


@dataclass
class QuasiStaticBody:
    """A planar chain whose velocity solves the force balance each step."""

    plan: BodyPlan
    drag: DragParameters = field(default_factory=DragParameters)
    muscle: MuscleParameters = field(default_factory=MuscleParameters)

    #: First node position, body heading, and joint angles.
    origin: np.ndarray = field(init=False)
    heading: float = field(init=False, default=0.0)
    joint_angles: np.ndarray = field(init=False)

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
            return c_par[:, None, None] * along + c_perp[:, None, None] * across
        return c_par * along + c_perp * across

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

        if self.drag.is_linear:
            return self._solve_with_limits(jac, self._drag_tensors(), forces)

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
        rates = self._solve_with_limits(jac, self._drag_tensors(self._last_velocities), forces)
        velocities = np.einsum("naj,j->na", jac, rates)
        for _ in range(_NONLINEAR_PASSES):
            previous = rates
            rates = self._solve_with_limits(jac, self._drag_tensors(velocities), forces)
            velocities = np.einsum("naj,j->na", jac, rates)
            if np.max(np.abs(rates - previous)) <= _NONLINEAR_TOL * max(
                float(np.max(np.abs(rates))), 1e-30
            ):
                break
        self._last_velocities = velocities
        return rates

    def _solve_with_limits(
        self, jac: np.ndarray, drag: np.ndarray, forces: np.ndarray
    ) -> np.ndarray:
        """One force balance, with joints at their limits locked out of it."""
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
            solved = np.linalg.solve(matrix[np.ix_(index, index)], forces[index])
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
        radii = np.array(
            [self.plan.radius_at(i) for i in range(self.plan.n_segments)],
            dtype=np.float64,
        )
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
