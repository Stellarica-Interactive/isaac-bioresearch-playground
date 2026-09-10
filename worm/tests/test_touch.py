"""Mechanosensory receptive fields.

The claims worth testing here are the ones taken from Chalfie et al. 1985: that a
touch near the head reaches ALM/AVM and a touch near the tail reaches PLM/PVM, and
that the ventral receptors answer only to ventral contact. The exact fractional
extents are our assumption and are deliberately *not* asserted to specific values
-- pinning those down would be testing our own guess.
"""

from __future__ import annotations

import numpy as np
import pytest

from worm.body.geometry import BodyPlan
from worm.body.touch import (
    HARSH_TOUCH_CELLS,
    NOSE_TOUCH_CELLS,
    TOUCH_RECEPTORS,
    TouchField,
)
from worm.loader import load


@pytest.fixture(scope="module")
def field() -> TouchField:
    return TouchField.build(BodyPlan())


def test_every_touch_receptor_exists_in_the_connectome() -> None:
    """A receptive field naming a cell the network does not contain would silently
    drop that whole branch of the circuit."""
    connectome, _ = load("cook_2019_herm")
    present = {c.id for c in connectome.cells}
    for group in (
        tuple(r.cell for r in TOUCH_RECEPTORS),
        HARSH_TOUCH_CELLS,
        NOSE_TOUCH_CELLS,
    ):
        assert set(group) <= present, f"missing from the connectome: {set(group) - present}"


def test_head_touch_reaches_the_anterior_receptors(field: TouchField) -> None:
    currents = field.currents(contact_fraction=0.15)
    assert set(currents) >= {"ALML", "ALMR"}
    assert "PLML" not in currents and "PLMR" not in currents


def test_tail_touch_reaches_the_posterior_receptors(field: TouchField) -> None:
    currents = field.currents(contact_fraction=0.9)
    assert set(currents) >= {"PLML", "PLMR"}
    assert "ALML" not in currents and "ALMR" not in currents


def test_midbody_recruits_both_and_weakly(field: TouchField) -> None:
    """The animal is least touch-sensitive around mid-body, because a contact there
    reaches both fields and the two responses oppose."""
    mid = field.currents(contact_fraction=0.5)
    head = field.currents(contact_fraction=0.15)
    assert {"ALML", "PLML"} <= set(mid)
    assert mid["ALML"] < head["ALML"]


def test_ventral_receptors_need_ventral_contact(field: TouchField) -> None:
    assert "AVM" not in field.currents(contact_fraction=0.2, ventral=False)
    assert "AVM" in field.currents(contact_fraction=0.2, ventral=True)
    assert "PVM" in field.currents(contact_fraction=0.7, ventral=True)


def test_response_is_graded_across_a_field(field: TouchField) -> None:
    """Peaks in the middle of the field and falls to zero at the edge, rather than
    switching on at an invented boundary."""
    values = [field.currents(contact_fraction=f).get("ALML", 0.0) for f in np.linspace(0, 0.6, 13)]
    assert values[0] < max(values)
    assert values[-1] == 0.0
    assert max(values) == pytest.approx(field.current_pa, rel=0.2)


def test_force_scales_the_response(field: TouchField) -> None:
    light = field.currents(contact_fraction=0.2, force=0.25)["ALML"]
    hard = field.currents(contact_fraction=0.2, force=1.0)["ALML"]
    assert hard == pytest.approx(4.0 * light)


def test_no_contact_injects_nothing(field: TouchField) -> None:
    """The probe sitting off to one side must be indistinguishable from no probe."""
    plan = BodyPlan()
    assert field.currents_from_segments(np.zeros(plan.n_segments)) == {}


def test_segment_contacts_take_the_strongest_not_the_sum(field: TouchField) -> None:
    """A blunt object covering six segments must not read as a bigger stimulus than
    a sharp one pressing just as hard on one."""
    plan = BodyPlan()
    sharp = np.zeros(plan.n_segments)
    sharp[3] = 1.0
    blunt = np.zeros(plan.n_segments)
    blunt[1:7] = 1.0

    one = field.currents_from_segments(sharp)["ALML"]
    many = field.currents_from_segments(blunt)["ALML"]

    # The blunt contact may read *higher* than this particular sharp one, because
    # spanning six segments includes one nearer the centre of ALML's field. What it
    # must never do is accumulate: six contacts are not six times the stimulus.
    individual = [
        field.currents_from_segments(np.eye(plan.n_segments)[i]).get("ALML", 0.0)
        for i in range(1, 7)
    ]
    assert many == pytest.approx(max(individual))
    assert many < 0.5 * sum(individual), "contacts are being summed, not maxed"
    assert many <= field.current_pa
    assert one > 0.0


def test_segment_contacts_reject_the_wrong_shape(field: TouchField) -> None:
    with pytest.raises(ValueError, match="per segment"):
        field.currents_from_segments(np.zeros(5))


def test_build_drops_receptors_absent_from_the_network() -> None:
    field = TouchField.build(BodyPlan(), cells=("ALML", "ALMR", "AVAL"))
    assert set(field.cells) == {"ALML", "ALMR"}
    assert "PLML" not in field.currents(contact_fraction=0.9)


def test_segment_fraction_spans_the_body() -> None:
    plan = BodyPlan()
    field = TouchField.build(plan)
    assert field.segment_fraction(0) < 0.1
    assert field.segment_fraction(plan.n_segments - 1) > 0.9
