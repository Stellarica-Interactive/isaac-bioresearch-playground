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

#: How far anterior of its own position a neuron senses curvature, in segments.
#: ASSUMED. Wen et al. describe B-type neurons responding to the curvature of the
#: region anterior to them; the size of that region in segment units is our
#: choice, and it largely sets the wavelength the loop settles into.
DEFAULT_SENSING_OFFSET = 2.0

#: Current injected per radian of sensed curvature, pA. ASSUMED, and the single
#: parameter that decides whether the loop oscillates at all.
DEFAULT_PROPRIOCEPTIVE_GAIN = 400.0


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
    """Which joint's curvature each target neuron reads."""

    gain_pa_per_rad: float = DEFAULT_PROPRIOCEPTIVE_GAIN

    @classmethod
    def build(
        cls,
        connectome: Connectome,
        plan: BodyPlan,
        *,
        offset: float = DEFAULT_SENSING_OFFSET,
        gain_pa_per_rad: float = DEFAULT_PROPRIOCEPTIVE_GAIN,
        classes: tuple[str, ...] = PROPRIOCEPTIVE_CLASSES,
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
        # Sense anterior to self; clamp at the head, where there is nothing ahead.
        sensed = np.clip(
            np.array([positions[c] - offset for c in chosen]).round().astype(int),
            0,
            plan.n_joints - 1,
        )
        return cls(tuple(chosen), sensed, gain_pa_per_rad)

    def currents(self, joint_angles_rad: np.ndarray) -> dict[str, float]:
        """Curvature to injected current, one entry per target neuron.

        DB neurons drive dorsal muscle and VB ventral, so the two read opposite
        signs of the same bend: each is excited by the body bending toward its own
        side, which is what makes the feedback recruit the next segment rather
        than fight it.
        """
        angles = np.asarray(joint_angles_rad, dtype=np.float64)
        out: dict[str, float] = {}
        for cell, segment in zip(self.targets, self.sensed_segment, strict=True):
            sign = 1.0 if cell.startswith("DB") else -1.0
            out[cell] = self.gain_pa_per_rad * sign * float(angles[segment])
        return out
