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

The ratio is genuinely uncertain, and ours is above measurement
--------------------------------------------------------------

* **Rabets et al. 2014** dragged real worms across agar with a microcantilever and
  measured the anisotropy at **3 to 10**.
* The values commonly used in neuromechanical models, following Niebur and Erdös
  1991, give a ratio nearer **40**.
* For swimming in water the ratio is only about 1.4 to 2, which is why worms swim
  with a different gait than they crawl.

:data:`DEFAULT_DRAG_RATIO` is **20**, which is twice the top of the directly
measured range. It sits inside the span of values *used* by models, not inside the
span of values *measured*. Swept in §5U: at a ratio of 10 -- inside the measured
range -- the model covers 0.702 BL against 0.994 BL at 20, so the measured value
makes this model worse, which is a fact about the model rather than about agar.

An earlier version of this docstring attributed absolute coefficients of "about
220 and 22" to Rabets et al. Those numbers could not be re-verified against any
accessible source, and are not repeated here. The anisotropy range above is
confirmed; the absolute coefficients need checking against the paper before being
quoted again.

What this model leaves out, which is the reason the body slides
--------------------------------------------------------------

Linear resistive force theory is the standard choice and what Boyle et al. and
c302 also use, but Rabets et al. found it is not what agar does. Substrate
viscoelasticity introduces **nonlinearities in the force-velocity relationship**,
giving **nonconstant** drag coefficients, and the major contributing factor is the
formation of a **shallow groove** in the surface.

A real worm crawls in a groove it has made. That groove is geometric confinement:
it blocks lateral slip in a way no velocity-proportional coefficient reproduces.
This model has a flat plane, a force linear in velocity, and **no static friction
at all**, so nothing holds the body in place when it is not actively pushing.
Sliding is the visible consequence, and it is a property of the model rather than
a bug in it.

:data:`DEFAULT_TANGENTIAL_DRAG` compounds this: only the *ratio* is constrained by
anything published, while the magnitude was chosen alongside the torque scale to
give a plausible speed. Two free parameters, one constraint between them.
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
#: Only the ratio to :data:`DEFAULT_DRAG_RATIO` is constrained by published work;
#: this magnitude and the torque scale are a pair of free parameters with one
#: constraint between them. See the module docstring.
DEFAULT_TANGENTIAL_DRAG = 2.0e-3

#: Exponent of the force-velocity law: ``F ~ |v|**exponent``. 1.0 is the linear
#: resistive force theory everything in this project was measured with, and is
#: the default so that nothing changes silently.
#:
#: Values below 1 make the medium **stiffer at low speed**, which is the
#: qualitative behaviour Rabets et al. 2014 measured on agar and linear theory
#: cannot express: they found the force-velocity relation nonlinear with
#: *nonconstant* coefficients, the major contributing factor being the shallow
#: groove the animal forms in the surface. A groove is lateral confinement, so
#: its signature is resistance that rises steeply as lateral speed approaches
#: zero -- exactly what a sublinear exponent gives.
#:
#: ASSUMED, and this is the weakest number in the model: Rabets established that
#: the relation is nonlinear and that the groove dominates, but **not** a power
#: law, and **not** an exponent. A power law is the simplest function with the
#: right limiting behaviour, not a fit to their data. See model_assumptions 5AH.
DEFAULT_DRAG_EXPONENT = 1.0

#: Speed at which the nonlinear law matches the linear one, m/s. Chosen as the
#: body's own characteristic lateral speed during crawling -- a 100 mm animal
#: undulating at 0.5 Hz with 16 degrees of joint amplitude moves its midbody
#: sideways at a few mm/s -- so that the exponent rescales the law around the
#: regime the animal actually occupies rather than around 1 m/s, which it never
#: reaches. ASSUMED.
DEFAULT_REFERENCE_SPEED = 3.0e-3

#: Floor on speed inside the nonlinear law, m/s. Without it the effective
#: coefficient diverges as a segment comes to rest and the force balance becomes
#: singular. Set well below the reference speed so it does not shape the regime
#: the animal occupies, and it also bounds the stiffening at
#: ``(eps / reference)**(exponent - 1)``.
NONLINEAR_SPEED_FLOOR = 1.0e-5

#: Yield force per segment, N. 0.0 is the frictionless default everything in
#: this project was measured with.
#:
#: This is the threshold the drag model has never had. Viscous drag inverts to
#: ``v = F / c``, so *any* residual force moves the body, and a nervous system
#: whose output never settles to exactly zero slides forever -- observed in the
#: viewport as a body in a stable shape creeping steadily across the ground.
#: Agar is a viscoelastic solid with a yield stress and a real worm below that
#: stress does not translate at all.
#:
#: ASSUMED, and weakly: Rabets et al. 2014 established that agar's response is
#: nonlinear and that groove formation dominates it, which is qualitative support
#: for a threshold existing. They did not publish a yield force for a worm on
#: agar, and this project has no measurement to set the magnitude from. It is
#: therefore swept rather than fitted, and the sweep is the result -- see
#: model_assumptions 5AH.
DEFAULT_YIELD_FORCE = 0.0

#: Speed scale over which the yield force turns on, m/s. Below it the yield term
#: acts like a viscous coefficient of ``yield_force / yield_speed``, which is
#: what makes a body below the threshold effectively stuck rather than exactly
#: stuck. Smaller is a harder threshold and a worse-conditioned solve.
DEFAULT_YIELD_SPEED = 1.0e-4


@dataclass(frozen=True, slots=True)
class DragParameters:
    tangential: float = DEFAULT_TANGENTIAL_DRAG
    ratio: float = DEFAULT_DRAG_RATIO
    exponent: float = DEFAULT_DRAG_EXPONENT
    reference_speed: float = DEFAULT_REFERENCE_SPEED
    yield_force: float = DEFAULT_YIELD_FORCE
    yield_speed: float = DEFAULT_YIELD_SPEED

    @property
    def perpendicular(self) -> float:
        return self.tangential * self.ratio

    @property
    def is_linear(self) -> bool:
        """Whether the force-velocity law is the linear one.

        Checked exactly rather than with a tolerance: the linear path is the one
        every committed measurement used, and it should be taken when and only
        when both departures from it are off.
        """
        return self.exponent == 1.0 and self.yield_force == 0.0

    def yield_coefficient(self, speed: np.ndarray) -> np.ndarray:
        """Extra drag coefficient from the yield force at a given speed.

        ``F_yield / (|v| + v_reg)``, so multiplying it by ``v`` gives a force
        that saturates at ``F_yield``. Expressed as a coefficient rather than a
        force so it can be added to the drag tensor and stay inside one linear
        solve per pass.
        """
        if self.yield_force == 0.0:
            return np.zeros_like(speed)
        return self.yield_force / (np.asarray(speed, dtype=np.float64) + self.yield_speed)

    def speed_factor(self, speed: np.ndarray) -> np.ndarray:
        """Multiplier on the drag coefficients at a given speed.

        ``(speed / reference) ** (exponent - 1)``, so it is 1 at the reference
        speed and, for an exponent below 1, larger below it. Returns ones for
        the linear law, which keeps the two paths textually identical rather
        than merely numerically close.
        """
        if self.is_linear:
            return np.ones_like(speed)
        floored = np.maximum(np.asarray(speed, dtype=np.float64), NONLINEAR_SPEED_FLOOR)
        return (floored / self.reference_speed) ** (self.exponent - 1.0)


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
