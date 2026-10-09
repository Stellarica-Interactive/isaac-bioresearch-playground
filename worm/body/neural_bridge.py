"""Between the connectome and the body: muscle drive out, proprioception back in.

This is the boundary the whole project is about, so every transformation across
it is named and justified here rather than buried in a simulation loop.

    connectome muscle cells --> quadrant + segment --> MuscleModel drive
    joint angles --> body curvature --> injected current into B-type neurons

Two directions, two different kinds of claim.

**Outward** is nearly free of assumption. Body wall muscles are named cells in
the connectome, their names encode quadrant and row, and the neuromuscular
junctions that drive them are measured synapses carrying a sign from measured
physiology. The only step we add is reading a muscle cell's synaptic activation
as its contraction.

**Inward** is an assumption, and a load-bearing one. *C. elegans* forward
locomotion is neuromechanical: there is no central pattern generator producing
the travelling wave on its own. B-type motor neurons are themselves
proprioceptive, sensing the curvature of the body immediately anterior to them,
and that feedback is what recruits each next segment and propagates the wave
backwards (Wen et al. 2012, Neuron 76:750-761). Without it the wave has nothing
to propagate it. With it, the body becomes part of the computation rather than a
display hung off the end of the nervous system.

What is genuinely derived rather than assumed is *where each motor neuron sits*.
Rather than assigning neurons to segments by hand, each one's body position is
computed from the muscles it actually innervates in the connectome, weighted by
synapse count. VB7's position comes from VB7's own measured synapses.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from common.data.schemas import CellCategory, Connectome, SynapseType
from worm.body.geometry import QUADRANTS, BodyPlan, parse_muscle, segment_of_muscle

#: Motor neuron classes carrying proprioceptive feedback. B-type neurons drive
#: forward locomotion and are the ones shown to be stretch sensitive; A-type drive
#: backward locomotion and are left out of a forward-crawling model.
PROPRIOCEPTIVE_CLASSES = ("VB", "DB")

#: How far anterior of its own position a neuron starts sensing, in segments.
#: ASSUMED. Wen et al. describe B-type neurons responding to the curvature of the
#: region anterior to them; the size of that region in segment units is our choice.
DEFAULT_SENSING_OFFSET = 2.0

#: The share of a cell's neuromuscular synapses that must go to one side for it
#: to count as driving that side. Every B-type, A-type and SMD cell is above 0.8;
#: the lateral RMDL and RMDR, at 0.5 and 0.67, are not.
SIDED_FRACTION = 0.75

#: How far the receptive field extends, as a fraction of the body. One joint.
#:
#: Boyle, Berri & Cohen 2012 use ``N_SR = M/2`` -- a stretch receptor field
#: spanning half the animal -- and their model undulates where this one does not,
#: so it was worth trying. **Measured here, it is worse:** half the body latches
#: this model, amp 0.00 against 0.70 degrees, and so does their dorsal asymmetry
#: independently. Either alone is enough to stop it.
#:
#: That is not evidence against their model. Their receptive field and asymmetry
#: are fitted together with binary hysteretic B-class neurons, their own muscle
#: model, a spring-rod body and their drag; lifting two components out of a tuned
#: system has no reason to work. The alternative reading -- that 0.70 degrees of
#: amplitude is a marginal instability rather than a mechanism -- is at least as
#: likely. Both settings stay reachable because the comparison is the useful part.
#: See docs/model_assumptions.md 5M.
DEFAULT_RECEPTIVE_FRACTION = 1 / 24

#: Directional asymmetry of the stretch response, from the same source, their
#: equation 13: ventral receptors respond linearly, dorsal ones are suppressed
#: when stretched and amplified when compressed. **Their fitted values, not
#: measured ones.** Included because a symmetric feedback law cannot distinguish
#: bending one way from the other, and the asymmetry is part of what breaks the
#: symmetry their model needs.
DORSAL_STRETCH_GAIN = 0.8
DORSAL_COMPRESS_GAIN = 1.2
VENTRAL_GAIN = 1.0

#: Current injected per radian of sensed curvature, pA. ASSUMED, and the single
#: parameter that decides whether the loop oscillates at all.
DEFAULT_PROPRIOCEPTIVE_GAIN = 400.0

#: Fraction of the proprioceptive signal that responds to the *rate of change* of
#: curvature rather than to curvature itself. **0.7 since model_assumptions 5AO.5**:
#: 70% rate, 30% curvature. 0.0, the pure curvature law, is what every result
#: before 5AO.5 was measured with (``--legacy-defaults`` in the runner).
#:
#: Why 0.7. The rate term is what lets the loop crawl forward at all (5AN.14-
#: 5AN.24): it makes the loop pass bends tailward at the gait frequency and
#: removes the static latch. The curvature term is *measured*: Wen et al. 2012
#: show B-type neurons hold their response while a bend is held. 30% is the
#: largest curvature share that still crawls robustly to perturbed starts (5AN.24).
#:
#: Why it exists: the curvature law is monotone positive feedback, so a static
#: bend is a stable fixed point of the loop -- which is what §5C.4's latch and
#: §5P.3's hand-over both are. A rate term has no fixed point at a static bend,
#: because a body that has stopped moving produces no drive. §5C.5 lists this as
#: option 3 and it had not been tried.
#:
#: ASSUMED, and a modelling choice rather than a parameter. Real mechanoreceptors
#: usually have both a phasic and a tonic component, so a blend is more
#: defensible than either pure law -- but **the fraction is ours and nothing
#: measures it.** See model_assumptions 5AK.
DEFAULT_RATE_FRACTION = 0.7

#: Lag over which curvature rate is measured, ms. A difference between
#: consecutive physics steps is mostly solver jitter; comparing against the body
#: a while ago is the same thing smoothed, without needing any state.
#:
#: ARBITRARY ENGINEERING: short against a 2 s gait so the rate is still local in
#: time, long against a 4.17 ms step so it is not reading noise.
RATE_LAG_MS = 50.0

#: Frequency at which the rate term is scaled to match the curvature term, Hz.
#:
#: Curvature is in radians and its rate in radians per second, so blending them
#: raw would change the signal's magnitude as well as its character: measured, a
#: pure rate law delivers 1117 pA where the tonic law delivers 195 for the same
#: motion. A run that then oscillated could not tell the law from a 5.7-fold gain
#: increase, which is the mistake §5Q records making with a binarised output.
#:
#: For a bend oscillating at ``f``, ``|dk/dt| = 2 pi f |k|``, so dividing the rate
#: by ``2 pi f`` puts it in the same units and at the same magnitude. 0.5 Hz is
#: the gait frequency used throughout (`sine_wave_drive`), drawn from observed
#: crawling on agar.
RATE_REFERENCE_HZ = 0.5


def muscle_slots(plan: BodyPlan, cell_ids: tuple[str, ...]) -> dict[str, tuple[int, int]]:
    """Map each body wall muscle cell to ``(quadrant index, segment index)``.

    Cells that are not body wall muscles are simply absent from the result, which
    is how pharyngeal, vulval and uterine muscles stay out of the body model
    without a special case anywhere.
    """
    out: dict[str, tuple[int, int]] = {}
    for cell_id in cell_ids:
        parsed = parse_muscle(cell_id)
        if parsed is None:
            continue
        quadrant, _ = parsed
        segment = segment_of_muscle(cell_id, plan)
        if segment is not None:
            out[cell_id] = (QUADRANTS.index(quadrant), segment)
    return out


def neuron_body_positions(connectome: Connectome, plan: BodyPlan) -> dict[str, float]:
    """Where along the body each motor neuron acts, derived from its own synapses.

    A neuron's position is the synapse-count-weighted mean segment of the body
    wall muscles it innervates. That is measured anatomy rather than an assignment
    we invent.

    Neurons innervating no body wall muscle are absent from the result.
    """
    slots = muscle_slots(plan, tuple(c.id for c in connectome.cells))
    weighted: dict[str, list[tuple[float, int]]] = {}
    for e in connectome.connections:
        if e.synapse_type is not SynapseType.CHEMICAL:
            continue
        slot = slots.get(e.post)
        if slot is None:
            continue
        weighted.setdefault(e.pre, []).append((float(e.weight), slot[1]))

    return {
        cell: float(np.average([s for _, s in rows], weights=[w for w, _ in rows]))
        for cell, rows in weighted.items()
        if rows
    }


@dataclass(frozen=True)
class MuscleDrive:
    """Reads muscle contraction out of the neural state.

    A muscle cell's synaptic activation is taken as its contraction. That is the
    one modelling step here: in the neuron model that variable means transmitter
    release, and for a muscle we reinterpret it as contraction. Both are a
    saturating function of membrane potential, so the shape is right; the
    identification is still ours.
    """

    plan: BodyPlan
    indices: np.ndarray
    """Index into the runtime's cell array, shape ``(4, n_segments)``; -1 where no
    muscle occupies that slot."""

    @classmethod
    def build(cls, plan: BodyPlan, cell_ids: tuple[str, ...]) -> MuscleDrive:
        indices = np.full((len(QUADRANTS), plan.n_segments), -1, dtype=np.int64)
        position = {cell: i for i, cell in enumerate(cell_ids)}
        for cell, (quadrant, segment) in muscle_slots(plan, cell_ids).items():
            indices[quadrant, segment] = position[cell]
        return cls(plan=plan, indices=indices)

    @property
    def covered(self) -> int:
        return int((self.indices >= 0).sum())

    def drive(self, activations: np.ndarray) -> np.ndarray:
        """Neural state to a ``(4, n_segments)`` drive for :class:`MuscleModel`.

        Slots with no muscle -- the animal genuinely lacks ``MVL24`` -- take their
        nearest neighbour's value along the body rather than zero, so a missing
        cell does not read as a commanded relaxation.
        """
        activations = np.asarray(activations, dtype=np.float64)
        out = np.zeros(self.indices.shape, dtype=np.float64)
        for q in range(self.indices.shape[0]):
            row = self.indices[q]
            present = np.flatnonzero(row >= 0)
            if present.size == 0:
                continue
            nearest = present[np.abs(np.arange(row.size)[:, None] - present).argmin(axis=1)]
            out[q] = activations[row[nearest]]
        return out


@dataclass(frozen=True)
class Proprioception:
    """Body curvature back into the nervous system.

    Without this the model has no way to propagate a wave: the animal has no
    central pattern generator for forward crawling, and the bend of one region is
    what recruits the next. See the module docstring.
    """

    targets: tuple[str, ...]
    """B-type motor neurons, ordered to match :attr:`sensed_segment`."""

    sensed_segment: np.ndarray
    """First joint of each target's receptive field, anterior-most."""

    gain_pa_per_rad: float = DEFAULT_PROPRIOCEPTIVE_GAIN

    receptive_joints: int = 1
    """How many joints each neuron integrates over, starting at
    :attr:`sensed_segment` and running posteriorly.

    One reproduces the original single-joint reading. Boyle, Berri & Cohen use
    half the body, and their model undulates where this one does not."""

    asymmetric: bool = False
    """Whether dorsal receptors use the directional gains of
    :data:`DORSAL_STRETCH_GAIN` and :data:`DORSAL_COMPRESS_GAIN`.

    Off by default: measured here, it latches the model. See
    :data:`DEFAULT_RECEPTIVE_FRACTION`."""

    rate_fraction: float = DEFAULT_RATE_FRACTION
    """How much of the signal responds to curvature *rate* instead of curvature.

    0 is the pure curvature law, which cannot help latching; 1 is purely phasic
    and cannot report a held posture at all. See :data:`DEFAULT_RATE_FRACTION`."""

    dorsal_targets: frozenset[str] | None = None
    """Targets that drive dorsal muscle, so are excited by dorsal bending.

    :meth:`build` fills this from each cell's own neuromuscular synapses in the
    connectome -- measured anatomy -- rather than from its name. ``None`` falls
    back to the name (``DB``, ``DA``, ``SMDD``), which is what a hand-built
    instance in a test uses."""

    @classmethod
    def build(
        cls,
        connectome: Connectome,
        plan: BodyPlan,
        *,
        offset: float = DEFAULT_SENSING_OFFSET,
        gain_pa_per_rad: float = DEFAULT_PROPRIOCEPTIVE_GAIN,
        classes: tuple[str, ...] = PROPRIOCEPTIVE_CLASSES,
        receptive_fraction: float = DEFAULT_RECEPTIVE_FRACTION,
        asymmetric: bool = False,
        rate_fraction: float = DEFAULT_RATE_FRACTION,
    ) -> Proprioception:
        positions = neuron_body_positions(connectome, plan)
        by_class = {
            c.id: c.class_name
            for c in connectome.cells
            if c.category is CellCategory.NEURON and c.class_name
        }
        chosen = sorted(
            cell
            for cell, position in positions.items()
            if by_class.get(cell) in classes and np.isfinite(position)
        )
        # Which side each cell drives, from its own synapses onto body-wall
        # muscle. A cell without a clear side -- RMDL drives 13 dorsal and 13
        # ventral -- cannot be excited by "bending toward its own side", so it is
        # left out rather than assigned one. Every B-type, A-type and SMD cell has
        # a side, and it agrees with its name.
        dorsal_syn: dict[str, int] = {c: 0 for c in chosen}
        ventral_syn: dict[str, int] = {c: 0 for c in chosen}
        for e in connectome.connections:
            if e.pre in dorsal_syn and e.synapse_type is SynapseType.CHEMICAL:
                if e.post.startswith("MD"):
                    dorsal_syn[e.pre] += e.weight
                elif e.post.startswith("MV"):
                    ventral_syn[e.pre] += e.weight
        chosen = [
            c
            for c in chosen
            if max(dorsal_syn[c], ventral_syn[c])
            >= SIDED_FRACTION * (dorsal_syn[c] + ventral_syn[c])
        ]
        dorsal = frozenset(c for c in chosen if dorsal_syn[c] > ventral_syn[c])
        # Sense anterior to self; clamp at the head, where there is nothing ahead.
        sensed = np.clip(
            np.array([positions[c] - offset for c in chosen]).round().astype(int),
            0,
            plan.n_joints - 1,
        )
        receptive = max(1, int(round(receptive_fraction * plan.n_joints)))
        return cls(
            tuple(chosen),
            sensed,
            gain_pa_per_rad,
            receptive_joints=receptive,
            asymmetric=asymmetric,
            rate_fraction=rate_fraction,
            dorsal_targets=dorsal,
        )

    def currents(
        self,
        joint_angles_rad: np.ndarray,
        *,
        earlier_angles_rad: np.ndarray | None = None,
        elapsed_ms: float | None = None,
    ) -> dict[str, float]:
        """Curvature to injected current, one entry per target neuron.

        With :attr:`rate_fraction` above zero the sensed quantity is a blend of
        curvature and its rate of change, and the rate needs a past posture to
        compare against: pass ``earlier_angles_rad`` from ``elapsed_ms`` ago, by
        preference :data:`RATE_LAG_MS`.

        Omitting them leaves the *rate* at zero, so the blend returns only its
        tonic share -- nothing at all for a purely phasic law. That is deliberate
        rather than a degenerate case: a phasic receptor with no motion history
        has genuinely nothing to report, and silently substituting the curvature
        law for the first :data:`RATE_LAG_MS` of a run would mean the loop briefly
        obeys a law the run was not configured with.

        Each neuron integrates curvature over a stretch of body starting anterior
        to itself, rather than reading a single joint. DB neurons drive dorsal
        muscle and VB ventral, so the two read opposite signs of the same bend:
        each is excited by the body bending toward its own side, which is what
        makes the feedback recruit the next segment rather than fight it.

        Dorsal receptors respond asymmetrically to stretch and compression
        (Boyle, Berri & Cohen 2012, eq. 13). A symmetric law cannot tell bending
        one way from the other.
        """
        angles = np.asarray(joint_angles_rad, dtype=np.float64)
        rate = np.zeros_like(angles)
        if (
            self.rate_fraction > 0.0
            and earlier_angles_rad is not None
            and elapsed_ms is not None
            and elapsed_ms > 0.0
        ):
            earlier = np.asarray(earlier_angles_rad, dtype=np.float64)
            if earlier.shape == angles.shape:
                # Scaled by 1 / (2 pi f) so a bend oscillating at the gait
                # frequency gives the same magnitude as the curvature it
                # replaces. Without this the blend would be a gain change as
                # well as a change of law, and the two could not be separated.
                rate = (angles - earlier) / (elapsed_ms / 1000.0)
                rate = rate / (2.0 * np.pi * RATE_REFERENCE_HZ)

        out: dict[str, float] = {}
        for cell, start in zip(self.targets, self.sensed_segment, strict=True):
            stop = int(start) + self.receptive_joints
            window = angles[int(start) : stop]
            rate_window = rate[int(start) : stop]
            if window.size == 0:
                window = angles[int(start) : int(start) + 1]
                rate_window = rate[int(start) : int(start) + 1]
            if self.rate_fraction > 0.0:
                # Blended in the sensed quantity rather than in the current, so
                # the asymmetric gains below apply to whatever is being sensed.
                window = (1.0 - self.rate_fraction) * window + self.rate_fraction * rate_window
            # DB, DA and SMDD drive dorsal muscle; VB, VA and SMDV ventral. Only
            # "DB" was tested while B-type cells were the only targets; an A-type
            # DA cell would otherwise have been read as ventral and given the
            # wrong sign. SMDD's side is measured twice over: its synapses in Cook
            # 2019 are 80 and 76 dorsal against 14 and 11 ventral, and its calcium
            # rises with dorsal head bending (Yeon et al. 2018).
            dorsal = (
                cell in self.dorsal_targets
                if self.dorsal_targets is not None
                else cell.startswith(("DB", "DA", "SMDD"))
            )
            sign = 1.0 if dorsal else -1.0
            if self.asymmetric and dorsal:
                gain = np.where(window > 0.0, DORSAL_STRETCH_GAIN, DORSAL_COMPRESS_GAIN)
            else:
                gain = np.full_like(window, VENTRAL_GAIN)
            out[cell] = float(self.gain_pa_per_rad * sign * np.mean(gain * window))
        return out
