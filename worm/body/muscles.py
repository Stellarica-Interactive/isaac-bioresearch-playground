"""From motor neuron activity to joint torque.

This is the transformation the brief insists on documenting: a biological motor
signal does not directly specify a torque, and every step between the two is a
modelling decision made here rather than buried in a simulator callback.

    motor neuron activation  s in [0, 1]     from the neural runtime
              |
              |  first-order activation kinetics (calcium is slow)
              v
    muscle activation        a in [0, 1]
              |
              |  antagonist difference: a muscle pulls, it cannot push
              v
    joint torque             tau (N m)

Three assumptions, all recorded in ``docs/model_assumptions.md``:

1. **Muscle activation lags neural activity.** Real muscle has calcium dynamics
   with a time constant of tens to a hundred milliseconds -- far slower than the
   neurons driving it. A first-order low-pass is the minimal model, and it
   matters: without it a jittery neural signal produces a jittery body rather
   than a smooth bend.

2. **Torque is the difference between antagonists.** A muscle generates tension
   and cannot push. Dorsal and ventral muscles oppose each other across a joint,
   so the net torque is proportional to their difference. Note what this means
   biologically: co-contraction of both sides produces no net bend but does
   stiffen the joint, which is real behaviour and falls out of the model rather
   than being added to it.

3. **Peak torque is a free parameter.** Nobody has measured the force a single
   *C. elegans* body wall muscle produces. :data:`DEFAULT_PEAK_TORQUE_SCALE` is
   set to make the model move at a plausible speed, which is fitting to a desired
   outcome and is labelled as such.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from worm.body.geometry import DORSAL, QUADRANTS, VENTRAL, BodyPlan

#: Time constant of muscle activation, in milliseconds. ASSUMED, order of
#: magnitude: C. elegans muscle calcium transients are on a scale of tens of ms.
DEFAULT_ACTIVATION_TAU_MS = 60.0

#: Torque per unit net activation, N m. ASSUMED -- tuned so the body moves at a
#: plausible speed, not measured. The most openly fitted number in the project.
#:
#: This is the value the documented W2 result was measured at: 0.109-0.119 BL/s
#: with 30.5 deg bends and body extent 0.60. It was briefly four times larger,
#: which pins every joint against its 60 deg limit and makes the body coil and
#: spin -- and, because the wave still looks superficially like a wave, that is
#: not obvious from watching it. The tell is `extent`, which falls to 0.21, and
#: `bend max`, which sits exactly on the limit.
DEFAULT_PEAK_TORQUE_SCALE = 5.0e-5

#: Passive stiffness and damping of the cuticle, resisting bending. ASSUMED.
#: Without these the chain is floppy and the wave does not hold its shape.
DEFAULT_JOINT_STIFFNESS = 1.0e-4
DEFAULT_JOINT_DAMPING = 2.0e-5


@dataclass(frozen=True, slots=True)
class MuscleParameters:
    activation_tau_ms: float = DEFAULT_ACTIVATION_TAU_MS
    peak_torque_scale: float = DEFAULT_PEAK_TORQUE_SCALE
    joint_stiffness: float = DEFAULT_JOINT_STIFFNESS
    joint_damping: float = DEFAULT_JOINT_DAMPING


class MuscleModel:
    """Muscle activation state, and the torque it produces.

    Stateful because activation lags its drive. One activation value per
    quadrant per segment, so the array is ``(4, n_segments)`` ordered by
    :data:`~worm.body.geometry.QUADRANTS`.
    """

    def __init__(self, plan: BodyPlan, params: MuscleParameters | None = None) -> None:
        self.plan = plan
        self.params = params or MuscleParameters()
        self.activation = np.zeros((len(QUADRANTS), plan.n_segments), dtype=np.float64)

    # -- state -------------------------------------------------------------

    def reset(self) -> None:
        self.activation[:] = 0.0

    def quadrant_index(self, quadrant: str) -> int:
        return QUADRANTS.index(quadrant)

    def dorsal_activation(self) -> np.ndarray:
        """Mean activation of the two dorsal quadrants, per segment."""
        return self.activation[[QUADRANTS.index(q) for q in DORSAL]].mean(axis=0)

    def ventral_activation(self) -> np.ndarray:
        return self.activation[[QUADRANTS.index(q) for q in VENTRAL]].mean(axis=0)

    # -- dynamics ----------------------------------------------------------

    def step(self, neural_drive: np.ndarray, dt_ms: float) -> None:
        """Advance muscle activation toward its neural drive.

        ``da/dt = (drive - a) / tau``, integrated with an exact exponential step
        so the result is stable at any ``dt`` -- the same reasoning as the neural
        runtime's integrator, and for the same reason.
        """
        if neural_drive.shape != self.activation.shape:
            raise ValueError(
                f"neural drive has shape {neural_drive.shape}, expected "
                f"{self.activation.shape} (quadrants x segments)"
            )
        alpha = 1.0 - np.exp(-dt_ms / self.params.activation_tau_ms)
        self.activation += alpha * (np.clip(neural_drive, 0.0, 1.0) - self.activation)

    # -- output ------------------------------------------------------------

    def joint_torques(
        self, joint_angles: np.ndarray | None = None, joint_velocities: np.ndarray | None = None
    ) -> np.ndarray:
        """Torque at each joint, in N m. Length ``plan.n_joints``.

        Positive torque bends the joint toward the dorsal side.

        A joint sits between segments ``i`` and ``i+1``, so it is driven by the
        muscles spanning it: the mean of the two segments' activations. Passive
        cuticle stiffness and damping are subtracted when joint state is supplied,
        which is what stops the chain behaving like a wet noodle.
        """
        dorsal = self.dorsal_activation()
        ventral = self.ventral_activation()
        net = dorsal - ventral

        # A joint is spanned by the muscles of the segments on either side.
        spanning = 0.5 * (net[:-1] + net[1:])
        torque = self.params.peak_torque_scale * spanning

        if joint_angles is not None:
            torque -= self.params.joint_stiffness * np.asarray(joint_angles, dtype=np.float64)
        if joint_velocities is not None:
            torque -= self.params.joint_damping * np.asarray(joint_velocities, dtype=np.float64)
        return torque

    def stiffness(self) -> np.ndarray:
        """Co-contraction per joint: how hard both antagonists pull at once.

        Not used to produce torque -- by construction it cancels -- but it is real
        and worth exposing. A worm that co-contracts is stiffening rather than
        bending, and the neural model can produce that without being asked to.
        """
        both = self.dorsal_activation() + self.ventral_activation()
        return 0.5 * (both[:-1] + both[1:])


def sine_wave_drive(
    plan: BodyPlan,
    t_ms: float,
    *,
    frequency_hz: float = 0.5,
    wavelength_fraction: float = 0.65,
    amplitude: float = 1.0,
) -> np.ndarray:
    """A hand-written travelling wave, for testing the body without neurons.

    **This is deliberately not biology.** It is the control: if the body crawls
    under a scripted wave and then fails to crawl under the connectome, the fault
    is in the neural model rather than the mechanics. It gets deleted from any
    experiment that claims a biological result.

    Defaults are drawn from observed forward crawling on agar -- roughly 0.5 Hz
    undulation with a little over one wavelength along the body -- so that the
    control is at least in the right regime.
    """
    phase_per_segment = 2.0 * np.pi / (wavelength_fraction * plan.n_segments)
    segments = np.arange(plan.n_segments)
    # Negative segment phase makes the wave travel head to tail, which is the
    # direction that drives forward locomotion.
    phase = 2.0 * np.pi * frequency_hz * (t_ms / 1000.0) - phase_per_segment * segments
    signal = amplitude * np.sin(phase)

    drive = np.zeros((len(QUADRANTS), plan.n_segments))
    # A muscle can only pull, so each side gets the half-wave rectified signal.
    for q in DORSAL:
        drive[QUADRANTS.index(q)] = np.clip(signal, 0.0, None)
    for q in VENTRAL:
        drive[QUADRANTS.index(q)] = np.clip(-signal, 0.0, None)
    return drive
