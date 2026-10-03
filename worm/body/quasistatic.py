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

**The joints are closed separately.** Their rates come from their own
overdamped dynamics, ``qdot = (torque - k q) / b``, rather than from the full
generalised force balance including the medium's resistance to bending. That is
defensible here because the joints are heavily overdamped on their own -- damping
ratio 7 for the armature and 49 for a segment, measured -- so their own damping
dominates. It is still a simplification: the medium resists shape change as well
as translation, and that resistance is not in the joint closure. The consequence
is that the body changes shape slightly too freely, which will overestimate
amplitude and therefore thrust.

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

#: Joint limit, radians. The same ceiling the Isaac articulation enforces, kept
#: here so the two solvers constrain the body identically -- a quasi-static body
#: free to fold past 60 degrees would not be comparable with the committed runs.
from worm.isaac_limits import JOINT_LIMIT_RAD


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

    def __post_init__(self) -> None:
        self.origin = np.zeros(2, dtype=np.float64)
        self.joint_angles = np.zeros(self.plan.n_joints, dtype=np.float64)

    # -- kinematics --------------------------------------------------------

    @property
    def segment_length_m(self) -> float:
        return self.plan.total_length_m / self.plan.n_segments

    def _directions(self) -> np.ndarray:
        """Unit heading of each segment, from the chain of joint angles."""
        # Segment 0 carries the body heading; each joint adds its angle.
        angles = self.heading + np.concatenate([[0.0], np.cumsum(self.joint_angles)])
        return np.stack([np.cos(angles), np.sin(angles)], axis=1)

    def nodes(self) -> np.ndarray:
        """The ``n_segments + 1`` joint positions along the body."""
        step = self.segment_length_m * self._directions()
        return np.vstack([self.origin, self.origin + np.cumsum(step, axis=0)])

    def segment_centres(self) -> np.ndarray:
        node = self.nodes()
        return 0.5 * (node[:-1] + node[1:])

    # -- dynamics ----------------------------------------------------------

    def _joint_rates(self, torques: np.ndarray) -> np.ndarray:
        """Overdamped joint dynamics: ``b qdot = torque - k q``.

        Closed separately from the body's rigid motion; see the module
        docstring for why that is defensible here and what it costs.
        """
        restoring = self.muscle.joint_stiffness * self.joint_angles
        return (np.asarray(torques, dtype=np.float64) - restoring) / (self.muscle.joint_damping)

    def _solve_rigid_motion(self, rates: np.ndarray) -> np.ndarray:
        """``(vx, vy, omega)`` such that no net force or torque acts.

        Each segment's velocity is ``A_i u + b_i``, linear in the rigid motion
        ``u`` with ``b_i`` from the shape change. Substituting into
        ``sum -D_i v_i = 0`` and ``sum r_i x (-D_i v_i) = 0`` gives a 3x3 system.
        """
        length = self.segment_length_m
        directions = self._directions()
        normals = _perp(directions)
        centres = self.segment_centres()
        # Positions relative to the first node, which is the reference point the
        # rotation column and the torque row are both written about.
        arms = centres - self.origin

        # Shape contribution. Segment i's angular rate from bending alone is the
        # sum of the joint rates before it; the node it starts from has already
        # accumulated the swing of every segment upstream.
        segment_rates = np.concatenate([[0.0], np.cumsum(rates)])
        swing = length * segment_rates[:, None] * normals
        node_shape = np.vstack([np.zeros((1, 2)), np.cumsum(swing, axis=0)[:-1]])
        b = node_shape + 0.5 * swing

        # Rotation about the reference point contributes perp(arm) per unit omega.
        rotation = _perp(arms)

        # Anisotropic drag as a tensor per segment.
        c_par, c_perp = self.drag.tangential, self.drag.perpendicular
        d_tensor = c_par * np.einsum("ij,ik->ijk", directions, directions) + (
            c_perp * np.einsum("ij,ik->ijk", normals, normals)
        )

        # A_i columns: translation in x, translation in y, rotation.
        a_cols = np.zeros((self.plan.n_segments, 2, 3), dtype=np.float64)
        a_cols[:, 0, 0] = 1.0
        a_cols[:, 1, 1] = 1.0
        a_cols[:, :, 2] = rotation

        # Force rows and the torque row, which is perp(arm) dotted into them.
        d_a = np.einsum("nij,njk->nik", d_tensor, a_cols)
        d_b = np.einsum("nij,nj->ni", d_tensor, b)
        matrix = np.zeros((3, 3), dtype=np.float64)
        matrix[:2] = d_a.sum(axis=0)
        matrix[2] = np.einsum("ni,nik->k", rotation, d_a)
        rhs = np.zeros(3, dtype=np.float64)
        rhs[:2] = -d_b.sum(axis=0)
        rhs[2] = -np.einsum("ni,ni->", rotation, d_b)

        return np.linalg.solve(matrix, rhs)

    def step(self, joint_torques: np.ndarray, dt_s: float) -> None:
        """Advance the configuration by ``dt_s``.

        No inertial stability limit: the only thing the step has to resolve is
        the change in shape, so 60 Hz against a two-second gait is ample.
        """
        rates = self._joint_rates(joint_torques)
        rigid = self._solve_rigid_motion(rates)
        if not np.all(np.isfinite(rigid)):
            raise FloatingPointError("force balance has no finite solution")

        self.origin = self.origin + dt_s * rigid[:2]
        self.heading = float(self.heading + dt_s * rigid[2])
        limit = JOINT_LIMIT_RAD
        self.joint_angles = np.clip(self.joint_angles + dt_s * rates, -limit, limit)

    # -- diagnostics -------------------------------------------------------

    def residual_force_and_torque(self, joint_torques: np.ndarray) -> tuple[float, float]:
        """How well the solved motion balances, for tests to assert on.

        Should be zero to solver precision. A non-zero residual means the body is
        being pushed by something that is not the ground.
        """
        rates = self._joint_rates(joint_torques)
        rigid = self._solve_rigid_motion(rates)
        directions = self._directions()
        normals = _perp(directions)
        arms = self.segment_centres() - self.origin
        c_par, c_perp = self.drag.tangential, self.drag.perpendicular
        d_tensor = c_par * np.einsum("ij,ik->ijk", directions, directions) + (
            c_perp * np.einsum("ij,ik->ijk", normals, normals)
        )
        length = self.segment_length_m
        segment_rates = np.concatenate([[0.0], np.cumsum(rates)])
        swing = length * segment_rates[:, None] * normals
        node_shape = np.vstack([np.zeros((1, 2)), np.cumsum(swing, axis=0)[:-1]])
        velocities = rigid[:2] + rigid[2] * _perp(arms) + node_shape + 0.5 * swing
        forces = -np.einsum("nij,nj->ni", d_tensor, velocities)
        net_force = float(np.linalg.norm(forces.sum(axis=0)))
        net_torque = float(np.sum(arms[:, 0] * forces[:, 1] - arms[:, 1] * forces[:, 0]))
        return net_force, net_torque
