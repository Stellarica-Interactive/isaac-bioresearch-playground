"""Anisotropic ground drag — the thing that turns wriggling into crawling.

Easy to leave out, and fatal if you do. An undulating chain on an isotropic
surface goes nowhere: every push forward is matched by an equal slip sideways
and the body oscillates in place. Forward motion requires the surface to resist
sideways motion **more** than it resists motion along the body.

Real worms get this from their cuticle, which is covered in circumferential
ridges (annuli) that grip laterally and slide longitudinally. On agar the effect
is large.

The model is resistive force theory: drag on each segment is linear in its
velocity, with different coefficients along and across the body axis.

    F = -( c_par * v_par * t_hat  +  c_perp * v_perp * n_hat )

This is applied as an explicit force per segment rather than left to contact
friction, because PhysX material friction is isotropic and cannot express the
distinction. It is also the approach the standard neuromechanical models take.

The ratio is genuinely uncertain
--------------------------------

Published estimates for crawling on agar disagree by roughly a factor of four,
which is worth knowing before treating any locomotion result as quantitative:

* Rabets et al. 2014 measure the normal and tangential coefficients directly and
  report about 220 and 22, a ratio near **10**.
* The values commonly used in neuromechanical models, following Niebur and Erdös
  1991, give a ratio nearer **40**.
* For swimming in water the ratio is only about 1.4 to 2, which is why worms swim
  with a different gait than they crawl.

So :data:`DEFAULT_DRAG_RATIO` is an ASSUMED parameter sitting inside a published
range, and it belongs in any sweep of a locomotion result.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from worm.body.geometry import BodyPlan

#: Ratio of perpendicular to parallel drag. ASSUMED; published agar estimates
#: span roughly 10 to 40. See the module docstring.
DEFAULT_DRAG_RATIO = 20.0

#: Tangential drag coefficient, N s / m per unit segment length, before scaling.
#: ASSUMED: chosen with the torque scale so the body moves at a plausible speed.
DEFAULT_TANGENTIAL_DRAG = 2.0e-3


@dataclass(frozen=True, slots=True)
class DragParameters:
    tangential: float = DEFAULT_TANGENTIAL_DRAG
    ratio: float = DEFAULT_DRAG_RATIO

    @property
    def perpendicular(self) -> float:
        return self.tangential * self.ratio


class GroundDrag:
    """Per-segment anisotropic drag, computed from segment poses and velocities.

    Pure numpy and simulator-independent: it takes positions and velocities and
    returns forces, so it can be tested against analytic cases without launching
    Isaac Sim.
    """

    def __init__(self, plan: BodyPlan, params: DragParameters | None = None) -> None:
        self.plan = plan
        self.params = params or DragParameters()
        # Longer segments meet more ground, so drag scales with segment length.
        lengths = np.array([s.length_m for s in plan.segments()])
        self._weight = lengths / lengths.mean()

    def tangents(self, positions: np.ndarray) -> np.ndarray:
        """Unit vector along the body at each segment, from segment centres.

        Uses the neighbouring centres so the tangent follows the body's actual
        curve rather than each segment's own orientation. At the ends it falls
        back to the one available neighbour.
        """
        positions = np.asarray(positions, dtype=np.float64)
        if positions.ndim != 2 or positions.shape[1] != 2:
            raise ValueError(f"positions must be (n_segments, 2), got {positions.shape}")

        diffs = np.empty_like(positions)
        diffs[1:-1] = positions[2:] - positions[:-2]
        diffs[0] = positions[1] - positions[0]
        diffs[-1] = positions[-1] - positions[-2]

        norms = np.linalg.norm(diffs, axis=1, keepdims=True)
        # A degenerate segment (coincident neighbours) gets an arbitrary but
        # finite tangent rather than a division by zero.
        safe = np.where(norms > 1e-12, norms, 1.0)
        tangent = diffs / safe
        tangent[norms[:, 0] <= 1e-12] = np.array([1.0, 0.0])
        return tangent

    def forces(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
        *,
        dt_s: float | None = None,
        masses_kg: np.ndarray | None = None,
    ) -> np.ndarray:
        """Drag force on each segment, shape ``(n_segments, 2)``, in newtons.

        Always opposes motion, so it can only remove kinetic energy -- asserted by
        the test suite, because a drag model that can add energy will quietly make
        a locomotion result meaningless.

        **Pass ``dt_s`` and ``masses_kg`` when driving a simulator.** Without them
        this returns the plain linear force ``-c v``, which a physics engine
        integrates explicitly -- and explicit damping is only stable while
        ``dt < 2m/c``. With segment masses around 1e-4 kg and a perpendicular
        coefficient of 0.04 that limit is about 5 ms, uncomfortably close to a
        240 Hz timestep; cross it and the velocity reverses and grows each step
        until the position is NaN.

        Given the step and the masses, the force instead reproduces the *exact*
        solution of ``m dv/dt = -c v`` over one step::

            v(t+dt) = v exp(-c dt / m)
            F_equiv = m (v(t+dt) - v) / dt

        which is unconditionally stable and can never reverse a velocity. Same
        reasoning as the neural runtime's exponential integrator, for the same
        reason: solve the linear part exactly rather than stepping toward it.
        """
        velocities = np.asarray(velocities, dtype=np.float64)
        if velocities.shape != positions.shape:
            raise ValueError(
                f"velocities {velocities.shape} must match positions {positions.shape}"
            )

        t_hat = self.tangents(positions)
        # Rotate the tangent 90 degrees in the plane to get the normal.
        n_hat = np.stack([-t_hat[:, 1], t_hat[:, 0]], axis=1)

        v_par = np.einsum("ij,ij->i", velocities, t_hat)
        v_perp = np.einsum("ij,ij->i", velocities, n_hat)

        c_par = self.params.tangential * self._weight
        c_perp = self.params.perpendicular * self._weight

        if dt_s is None or masses_kg is None:
            return -(c_par * v_par)[:, None] * t_hat - (c_perp * v_perp)[:, None] * n_hat

        m = np.asarray(masses_kg, dtype=np.float64)
        if m.shape != (positions.shape[0],):
            raise ValueError(f"masses must be ({positions.shape[0]},), got {m.shape}")

        def equivalent_force(c: np.ndarray, v: np.ndarray) -> np.ndarray:
            # m (v e^{-c dt/m} - v) / dt, written to stay finite as c -> 0.
            return m * v * np.expm1(-c * dt_s / m) / dt_s

        return (
            equivalent_force(c_par, v_par)[:, None] * t_hat
            + equivalent_force(c_perp, v_perp)[:, None] * n_hat
        )

    def describe(self) -> str:
        return (
            f"anisotropic ground drag: tangential {self.params.tangential:.3g}, "
            f"perpendicular {self.params.perpendicular:.3g} "
            f"(ratio {self.params.ratio:g})"
        )
