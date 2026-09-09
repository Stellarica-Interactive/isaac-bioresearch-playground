"""The connectome-to-body boundary.

These tests care about one thing above all: that a motor neuron's position along
the body is *derived from its own measured synapses* and not assigned by us. The
strongest evidence for that is ordering -- DB1 anterior to DB7, VB1 anterior to
VB11 -- because nothing in the code knows that the number in a cell's name means
anything. If the ordering comes out right, it came out of the anatomy.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from common.data.schemas import (
    Cell,
    CellCategory,
    Connection,
    Connectome,
    Provenance,
    Scope,
    SynapseType,
)
from worm.body.geometry import QUADRANTS, BodyPlan
from worm.body.neural_bridge import (
    MuscleDrive,
    Proprioception,
    muscle_slots,
    neuron_body_positions,
)
from worm.loader import load

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")


# -- synthetic fixtures, so the unit behaviour is checkable without a dataset ---

PROV = Provenance(
    source_id="toy",
    kind="connectome",
    citation="synthetic fixture, not a real animal",
    url="",
    license="n/a",
)


def _toy() -> Connectome:
    """Two motor neurons, one clearly anterior, one clearly posterior."""
    muscles = ["MDL01", "MDL02", "MDL20", "MDL21"]
    cells = [
        Cell(id="NA", category=CellCategory.NEURON, class_name="NA"),
        Cell(id="NP", category=CellCategory.NEURON, class_name="NP"),
        *[Cell(id=m, category=CellCategory.MUSCLE) for m in muscles],
    ]
    connections = [
        Connection(pre="NA", post="MDL01", synapse_type=SynapseType.CHEMICAL, weight=3),
        Connection(pre="NA", post="MDL02", synapse_type=SynapseType.CHEMICAL, weight=1),
        Connection(pre="NP", post="MDL20", synapse_type=SynapseType.CHEMICAL, weight=1),
        Connection(pre="NP", post="MDL21", synapse_type=SynapseType.CHEMICAL, weight=1),
    ]
    return Connectome(
        id="toy",
        organism="C. elegans",
        sex="hermaphrodite",
        stage="adult",
        scope=Scope.WHOLE_ANIMAL,
        cells=tuple(cells),
        connections=tuple(connections),
        provenance=PROV,
    )


def test_muscle_slots_ignores_non_body_wall_cells() -> None:
    plan = BodyPlan()
    slots = muscle_slots(plan, ("MDL01", "AVAL", "pm1", "um1AL", "MVR12"))
    assert set(slots) == {"MDL01", "MVR12"}
    quadrant, segment = slots["MDL01"]
    assert QUADRANTS[quadrant] == "DL"
    assert segment == 0


def test_position_is_the_synapse_weighted_mean() -> None:
    """Not the mean of the muscles -- the mean weighted by how many synapses."""
    positions = neuron_body_positions(_toy(), BodyPlan())
    # MDL01 in segment 0 with weight 3, MDL02 in segment 1 with weight 1.
    assert positions["NA"] == pytest.approx(0.25)
    assert positions["NP"] > positions["NA"]


def test_electrical_connections_do_not_place_a_neuron() -> None:
    """A gap junction to a muscle is not a command, so it must not count."""
    base = _toy()
    with_gap = replace(
        base,
        connections=(
            *base.connections,
            Connection(pre="NA", post="MDL21", synapse_type=SynapseType.ELECTRICAL, weight=50),
        ),
    )
    assert neuron_body_positions(with_gap, BodyPlan())["NA"] == pytest.approx(0.25)


# -- MuscleDrive ---------------------------------------------------------------


def test_drive_routes_each_muscle_to_its_own_slot() -> None:
    plan = BodyPlan()
    cell_ids = ("AVAL", "MDL01", "MVR01")
    drive = MuscleDrive.build(plan, cell_ids)
    out = drive.drive(np.array([0.9, 0.4, 0.7]))

    assert out.shape == (len(QUADRANTS), plan.n_segments)
    assert out[QUADRANTS.index("DL"), 0] == pytest.approx(0.4)
    assert out[QUADRANTS.index("VR"), 0] == pytest.approx(0.7)
    # AVAL is not a muscle and must not appear anywhere in the body drive.
    assert not np.any(np.isclose(out, 0.9))


def test_empty_quadrant_stays_zero_rather_than_borrowing() -> None:
    plan = BodyPlan()
    drive = MuscleDrive.build(plan, ("MDL01",))
    out = drive.drive(np.array([0.8]))
    assert np.all(out[QUADRANTS.index("DL")] == pytest.approx(0.8))
    assert np.all(out[QUADRANTS.index("VR")] == 0.0)


def test_drive_output_is_in_range_for_activations_in_range() -> None:
    plan = BodyPlan()
    connectome, _ = load("cook_2019_herm")
    cell_ids = tuple(c.id for c in connectome.cells)
    drive = MuscleDrive.build(plan, cell_ids)
    rng = np.random.default_rng(0)
    out = drive.drive(rng.uniform(0.0, 1.0, size=len(cell_ids)))
    assert out.min() >= 0.0 and out.max() <= 1.0


# -- Proprioception ------------------------------------------------------------


def test_feedback_sign_is_opposite_for_dorsal_and_ventral() -> None:
    """DB and VB read opposite signs of the same bend.

    They drive antagonist muscle, so identical signs would mean the feedback
    fights itself rather than recruiting the next segment.
    """
    proprio = Proprioception(
        targets=("DB3", "VB4"), sensed_segment=np.array([5, 5]), gain_pa_per_rad=100.0
    )
    angles = np.zeros(BodyPlan().n_joints)
    angles[5] = 0.2
    currents = proprio.currents(angles)
    assert currents["DB3"] == pytest.approx(20.0)
    assert currents["VB4"] == pytest.approx(-20.0)


def test_a_straight_body_injects_nothing() -> None:
    """The loop must be quiescent at rest, or the wave is coming from an offset
    we added rather than from the body."""
    proprio = Proprioception(("DB1",), np.array([0]), 400.0)
    assert proprio.currents(np.zeros(BodyPlan().n_joints))["DB1"] == 0.0


def test_sensed_segment_is_always_a_valid_joint() -> None:
    plan = BodyPlan()
    connectome, _ = load("cook_2019_herm")
    proprio = Proprioception.build(connectome, plan)
    assert proprio.sensed_segment.min() >= 0
    assert proprio.sensed_segment.max() < plan.n_joints


# -- the real dataset ----------------------------------------------------------


def test_derived_positions_are_ordered_head_to_tail() -> None:
    """The load-bearing result.

    Nothing in the code parses the number out of a motor neuron's name. If DB1
    lands anterior to DB7 it is because DB1's measured synapses go to anterior
    muscles, which is the claim the whole bridge rests on.
    """
    connectome, _ = load("cook_2019_herm")
    positions = neuron_body_positions(connectome, BodyPlan())

    for prefix, count in (("DB", 7), ("VB", 11), ("DA", 9), ("VA", 12)):
        members = [f"{prefix}{i}" for i in range(1, count + 1)]
        present = [(c, positions[c]) for c in members if c in positions]
        assert len(present) >= count - 1, f"{prefix} mostly missing from the connectome"
        ordered = [p for _, p in present]
        assert np.corrcoef(np.arange(len(ordered)), ordered)[0, 1] > 0.9, (
            f"{prefix} positions do not run head to tail: {present}"
        )


def test_motor_neurons_span_the_body() -> None:
    connectome, _ = load("cook_2019_herm")
    positions = neuron_body_positions(connectome, BodyPlan())
    values = np.array(list(positions.values()))
    assert values.min() < 4.0, "no motor neuron near the head"
    assert values.max() > 19.0, "no motor neuron near the tail"


def test_body_wall_muscle_coverage_is_essentially_complete() -> None:
    """95 of 96 slots. The animal genuinely lacks a 24th ventral-left muscle:
    MVL has 23 cells where the other three quadrants have 24."""
    plan = BodyPlan()
    connectome, _ = load("cook_2019_herm")
    drive = MuscleDrive.build(plan, tuple(c.id for c in connectome.cells))
    assert drive.covered == 95

    missing = np.argwhere(drive.indices < 0)
    assert len(missing) == 1
    quadrant, segment = missing[0]
    assert QUADRANTS[quadrant] == "VL"
    assert segment == plan.n_segments - 1


def test_proprioceptive_targets_are_b_type_only() -> None:
    """A-type motor neurons drive backward locomotion and must not be in a
    forward-crawling feedback loop."""
    connectome, _ = load("cook_2019_herm")
    proprio = Proprioception.build(connectome, BodyPlan())
    assert proprio.targets
    assert all(c.startswith(("DB", "VB")) for c in proprio.targets)
    assert len(proprio.targets) == len(proprio.sensed_segment)


def test_sensing_offset_moves_the_sensed_region_anteriorly() -> None:
    connectome, plan = load("cook_2019_herm")[0], BodyPlan()
    near = Proprioception.build(connectome, plan, offset=0.0)
    far = Proprioception.build(connectome, plan, offset=4.0)
    assert near.targets == far.targets
    # Clamping at the head means some entries tie; none may move posteriorly.
    assert np.all(far.sensed_segment <= near.sensed_segment)
    assert np.any(far.sensed_segment < near.sensed_segment)
