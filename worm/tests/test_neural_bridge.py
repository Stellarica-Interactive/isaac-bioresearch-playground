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
    DEFAULT_PROPRIOCEPTIVE_GAIN,
    RATE_LAG_MS,
    RATE_REFERENCE_HZ,
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
        targets=("DB3", "VB4"),
        sensed_segment=np.array([5, 5]),
        gain_pa_per_rad=100.0,
        asymmetric=False,
        rate_fraction=0.0,  # the curvature law's sign, not the default blend's
    )
    angles = np.zeros(BodyPlan().n_joints)
    angles[5] = 0.2
    currents = proprio.currents(angles)
    assert currents["DB3"] == pytest.approx(20.0)
    assert currents["VB4"] == pytest.approx(-20.0)


def test_dorsal_receptors_respond_asymmetrically() -> None:
    """Suppressed when stretched, amplified when compressed.

    Boyle, Berri & Cohen 2012 eq. 13, and their fitted values rather than measured
    ones. A symmetric law cannot tell bending one way from the other, and the
    asymmetry is part of what lets their model break symmetry and undulate.
    """
    proprio = Proprioception(
        targets=("DB3", "VB4"),
        sensed_segment=np.array([5, 5]),
        gain_pa_per_rad=100.0,
        asymmetric=True,
        rate_fraction=0.0,  # the curvature law's asymmetry, not the default blend's
    )
    plan = BodyPlan()
    bent, folded = np.zeros(plan.n_joints), np.zeros(plan.n_joints)
    bent[5], folded[5] = 0.2, -0.2

    assert proprio.currents(bent)["DB3"] == pytest.approx(100.0 * 0.8 * 0.2)
    assert proprio.currents(folded)["DB3"] == pytest.approx(-100.0 * 1.2 * 0.2)
    # Ventral receptors stay linear, so the two sides are not mirror images.
    assert proprio.currents(bent)["VB4"] == pytest.approx(-20.0)
    assert abs(proprio.currents(bent)["DB3"]) != pytest.approx(abs(proprio.currents(bent)["VB4"]))


def test_receptive_field_spans_a_stretch_of_body() -> None:
    """Half the body, following Boyle et al.'s N_SR = M/2.

    This model previously read one joint two segments ahead. In a feedback loop
    that is not a detail: one segment reports local curvature, half a body reports
    the shape of the wave.
    """
    connectome, _ = load("cook_2019_herm")
    plan = BodyPlan()
    wide = Proprioception.build(connectome, plan, receptive_fraction=0.5)
    assert wide.receptive_joints == pytest.approx(plan.n_joints // 2, abs=1)

    narrow = Proprioception.build(connectome, plan)
    assert narrow.receptive_joints == 1, "the default is one joint; see 5M"

    # A bend outside the window must not be felt; one inside it must.
    angles = np.zeros(plan.n_joints)
    angles[wide.sensed_segment[0] + 1] = 0.3
    assert wide.currents(angles)[wide.targets[0]] != 0.0
    assert narrow.currents(angles)[narrow.targets[0]] == 0.0


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


# -- the phasic proprioceptive term -----------------------------------------
#
# The curvature law is monotone positive feedback: a dorsal bend excites DB, DB
# deepens the bend, and a static bend is a stable fixed point of the loop. That
# fixed point is the latch of model_assumptions 5C.4 and the settled body of
# 5P.3. A law responding to curvature *rate* has no fixed point there, because a
# body that has stopped moving produces no drive. 5C.5 lists this as option 3.


def _probe(rate_fraction: float) -> Proprioception:
    """Two targets reading the same joints, so only the law is under test."""
    return Proprioception(
        ("DB1", "VB1"),
        np.array([5, 5]),
        DEFAULT_PROPRIOCEPTIVE_GAIN,
        receptive_joints=3,
        rate_fraction=rate_fraction,
    )


def test_the_default_law_is_mostly_rate_with_a_measured_tonic_part() -> None:
    """The committed law since model_assumptions 5AO.5, and deliberately pinned.

    This test used to require 0.0, the pure curvature law, because every result
    up to then assumed it. Changing it was a decision, recorded in 5AO.5: 70% rate
    so the loop can crawl, 30% curvature because Wen et al. 2012 measured a tonic
    response. A future change should be as deliberate, so it fails here first.
    """
    assert Proprioception(("DB1",), np.array([0]), 1.0).rate_fraction == 0.7


def test_a_held_bend_drives_the_tonic_law_and_not_the_phasic_one() -> None:
    """The fixed point being removed, which is the whole point.

    A body holding a bend reports it forever under the curvature law -- that is
    what sustains the latch -- and reports nothing at all under the rate law.
    """
    held = np.full(23, np.radians(20.0))
    tonic = _probe(0.0).currents(held, earlier_angles_rad=held, elapsed_ms=RATE_LAG_MS)
    phasic = _probe(1.0).currents(held, earlier_angles_rad=held, elapsed_ms=RATE_LAG_MS)
    assert abs(tonic["DB1"]) > 100.0
    assert phasic["DB1"] == pytest.approx(0.0)


def test_the_rate_term_is_gain_matched_at_the_gait_frequency() -> None:
    """The blend must change the law without changing the gain.

    Raw, the rate term is 5.7x the tonic one for the same motion, because
    radians per second are numerically larger than radians. A phasic run that
    then oscillated could not be told from a 5.7-fold gain increase -- which is
    the mistake 5Q records making with a binarised output. Scaling by
    1/(2 pi f) puts them in the same units at the gait's own frequency.
    """
    bend = np.full(23, np.radians(20.0))
    step = bend * np.sin(2.0 * np.pi * RATE_REFERENCE_HZ * RATE_LAG_MS / 1000.0)
    currents = [
        _probe(f).currents(bend + step, earlier_angles_rad=bend, elapsed_ms=RATE_LAG_MS)["DB1"]
        for f in (0.0, 0.5, 1.0)
    ]
    # Within 20%: exact equality is not expected, since a finite lag on a sine
    # is not the derivative, but a factor of 5.7 would be a different experiment.
    assert max(currents) < 1.2 * min(currents), f"not gain-matched: {currents}"


def test_without_history_only_the_tonic_share_is_reported() -> None:
    """A phasic receptor with no motion history has nothing to report.

    Written first as "degrades to the tonic law", which is what it felt like it
    should do and is not what it does: with no past posture the rate is zero, so
    the blend returns only its tonic share -- nothing at all when that share is
    zero. The behaviour is the right one. Substituting the curvature law for the
    first RATE_LAG_MS would mean the loop briefly obeys a law the run was not
    configured with, which is worse than a short transient.
    """
    bend = np.full(23, np.radians(15.0))
    tonic = _probe(0.0).currents(bend)["DB1"]

    assert _probe(1.0).currents(bend)["DB1"] == pytest.approx(0.0)
    assert _probe(1.0).currents(bend, earlier_angles_rad=bend)["DB1"] == pytest.approx(0.0)
    # A half-and-half law keeps half of it, which is the same rule.
    assert _probe(0.5).currents(bend)["DB1"] == pytest.approx(0.5 * tonic)
    # And nothing raises for the missing arguments.
    assert _probe(0.5).currents(bend, elapsed_ms=RATE_LAG_MS)["DB1"] == pytest.approx(0.5 * tonic)


def test_a_type_cells_take_the_side_they_drive() -> None:
    """DA drives dorsal muscle as DB does, VA ventral as VB does.

    The sign used to be chosen by the prefix "DB" alone, which was right while
    B-type cells were the only targets and would have read every DA cell as
    ventral -- excited by the wrong bend -- the moment A-type cells were added
    (model_assumptions 5AO.2).
    """
    proprio = Proprioception(
        targets=("DA1", "VA1", "DB1", "VB1"),
        sensed_segment=np.array([5, 5, 5, 5]),
        gain_pa_per_rad=1.0,
    )
    angles = np.zeros(23)
    angles[5] = 0.1  # dorsal
    out = proprio.currents(angles)
    assert out["DA1"] > 0.0 and out["DB1"] > 0.0
    assert out["VA1"] < 0.0 and out["VB1"] < 0.0
    assert out["DA1"] == out["DB1"] and out["VA1"] == out["VB1"]


def test_build_takes_each_cells_side_from_its_own_synapses() -> None:
    """Measured anatomy, not names, decides which bend excites a cell.

    The name rule was right for B-type cells and would have been silently wrong
    for any head neuron whose name does not start with DB; model_assumptions
    5AO.7. RMDL and RMDR drive both sides about equally, so they have no "own
    side" and are left out rather than assigned one.
    """
    connectome, _ = load("cook_2019_herm")
    plan = BodyPlan()
    b_type = Proprioception.build(connectome, plan)
    assert b_type.dorsal_targets == frozenset(c for c in b_type.targets if c.startswith("DB"))

    head = Proprioception.build(connectome, plan, classes=("SMD", "RMD"))
    assert {"SMDDL", "SMDDR", "RMDDL", "RMDDR"} <= head.dorsal_targets
    assert {"SMDVL", "SMDVR", "RMDVL", "RMDVR"} <= set(head.targets) - head.dorsal_targets
    assert "RMDL" not in head.targets and "RMDR" not in head.targets


def test_ventral_scale_weights_only_the_ventral_cells() -> None:
    """The probe of the head's dorsal-ventral balance is a gain on the ventral side.

    At 1.0 (the default) both sides respond equally; any other value scales the
    ventral cells' current and leaves the dorsal cells' untouched, so the
    parameter cannot change which bend excites a cell (model_assumptions 5AO.9).
    """
    angles = np.zeros(23)
    angles[5] = 0.1
    out = {}
    for scale in (1.0, 0.5):
        proprio = Proprioception(
            targets=("SMDDL", "SMDVL"),
            sensed_segment=np.array([5, 5]),
            gain_pa_per_rad=1.0,
            dorsal_targets=frozenset({"SMDDL"}),
            ventral_scale=scale,
        )
        out[scale] = proprio.currents(angles)
    assert out[1.0]["SMDDL"] == pytest.approx(-out[1.0]["SMDVL"])
    assert out[0.5]["SMDDL"] == out[1.0]["SMDDL"]
    assert out[0.5]["SMDVL"] == pytest.approx(0.5 * out[1.0]["SMDVL"])
