"""Chemosensation.

The claim these tests defend is that the model responds to *change* in
concentration rather than to concentration, because that is the mechanism the
animal actually uses and the easiest thing to get wrong.
"""

from __future__ import annotations

import numpy as np
import pytest

from worm.body.chemotaxis import (
    CHEMOSENSORS,
    MODALITIES,
    Chemosensation,
    FoodSource,
)
from worm.loader import load


@pytest.fixture
def sense() -> Chemosensation:
    """All six sensors. The default is one modality (see the module); these
    tests are checking the ON/OFF machinery, not which pathway is selected."""
    return Chemosensation.build(modality=None)


# -- the gradient --------------------------------------------------------------


def test_concentration_peaks_on_the_source() -> None:
    food = FoodSource(x_m=0.2, y_m=0.0, strength=1.0)
    assert food.concentration(np.array([0.2, 0.0])) == pytest.approx(1.0)
    assert food.concentration(np.array([0.2, 0.05])) < 1.0


def test_concentration_falls_with_distance() -> None:
    food = FoodSource(x_m=0.0, y_m=0.0)
    values = [food.concentration(np.array([d, 0.0])) for d in (0.0, 0.1, 0.2, 0.4)]
    assert values == sorted(values, reverse=True)


def test_decay_length_means_what_it_says() -> None:
    """One decay length out, concentration is 1/e. Worth pinning: the default is
    three body lengths, so the field is still a quarter strength four body lengths
    away -- which is a deliberately broad gradient, not a bug."""
    food = FoodSource(x_m=0.0, y_m=0.0, decay_m=0.3)
    assert food.concentration(np.array([0.3, 0.0])) == pytest.approx(1 / np.e, rel=1e-6)
    assert food.concentration(np.array([0.4, 0.0])) == pytest.approx(0.264, abs=0.01)


# -- the mechanism -------------------------------------------------------------


def test_a_steady_field_produces_no_response(sense: Chemosensation) -> None:
    """The point of the whole module.

    A worm sitting in a uniform concentration is not being told anything, however
    strong that concentration is. A model keyed to absolute concentration would
    inject current forever and describe a different animal.
    """
    for _ in range(2000):
        sense.step(0.8, dt_ms=5.0)
    assert sense.currents() == {} or max(abs(v) for v in sense.currents().values()) < 0.01


def test_arriving_in_a_field_is_not_itself_a_stimulus(sense: Chemosensation) -> None:
    """Starting adapted, so the first sample does not read as a huge step."""
    sense.step(0.9, dt_ms=5.0)
    assert not sense.currents()


def test_rising_concentration_excites_the_on_cells(sense: Chemosensation) -> None:
    sense.step(0.1, dt_ms=5.0)
    for _ in range(40):
        sense.step(0.9, dt_ms=5.0)
    currents = sense.currents()
    assert currents["AWAL"] > 0.0, "AWA is an ON cell"
    assert currents["ASEL"] > 0.0, "ASEL responds to increases"
    assert currents["AWCL"] < 0.0, "AWC is an OFF cell and must be inhibited here"
    assert currents["ASER"] < 0.0, "ASER responds to decreases"


def test_falling_concentration_reverses_every_sign(sense: Chemosensation) -> None:
    """AWC firing on removal is the measured behaviour, and the asymmetry between
    ASEL and ASER is the reason both are in the model."""
    sense.step(0.9, dt_ms=5.0)
    for _ in range(40):
        sense.step(0.1, dt_ms=5.0)
    currents = sense.currents()
    assert currents["AWCL"] > 0.0, "AWC fires on odour removal"
    assert currents["ASER"] > 0.0
    assert currents["AWAL"] < 0.0
    assert currents["ASEL"] < 0.0


def test_response_decays_as_the_cell_adapts(sense: Chemosensation) -> None:
    sense.step(0.1, dt_ms=5.0)
    sense.step(0.9, dt_ms=5.0)
    immediate = abs(sense.currents()["AWAL"])
    for _ in range(3000):
        sense.step(0.9, dt_ms=5.0)
    later = abs(sense.currents().get("AWAL", 0.0))
    assert later < 0.1 * immediate


def test_bigger_change_gives_bigger_response(sense: Chemosensation) -> None:
    small = Chemosensation.build(modality=None)
    large = Chemosensation.build(modality=None)
    for s, target in ((small, 0.3), (large, 0.9)):
        s.step(0.1, dt_ms=5.0)
        s.step(target, dt_ms=5.0)
    assert abs(large.currents()["AWAL"]) > abs(small.currents()["AWAL"])


# -- wiring --------------------------------------------------------------------


def test_every_chemosensor_exists_in_the_connectome() -> None:
    connectome, _ = load("cook_2019_herm")
    present = {c.id for c in connectome.cells}
    missing = {s.cell for s in CHEMOSENSORS} - present
    assert not missing, f"not in the connectome: {missing}"


def test_on_and_off_cells_are_both_represented() -> None:
    """A model with only ON cells could not report a decrease, which is half the
    information a biased random walk runs on."""
    signs = {s.sign for s in CHEMOSENSORS}
    assert signs == {+1.0, -1.0}


def test_build_drops_sensors_absent_from_the_network() -> None:
    sense = Chemosensation.build(cells=("AWAL", "AWCL", "AVAL"), modality=None)
    assert set(sense.cells) == {"AWAL", "AWCL"}


def test_per_cell_currents_override_the_flat_scale() -> None:
    sense = Chemosensation.build(per_cell_pa={"AWAL": 100.0}, modality=None)
    sense.step(0.0, dt_ms=5.0)
    sense.step(1.0, dt_ms=5.0)
    currents = sense.currents()
    assert currents["AWAL"] == pytest.approx(100.0 * sense.deviation)
    assert currents["AWAR"] != currents["AWAL"]


def test_reset_forgets_the_history(sense: Chemosensation) -> None:
    sense.step(0.1, dt_ms=5.0)
    sense.step(0.9, dt_ms=5.0)
    assert sense.currents()
    sense.reset()
    sense.step(0.9, dt_ms=5.0)
    assert not sense.currents()


# -- modality ------------------------------------------------------------------


def test_default_modality_is_one_pathway_not_all_of_them() -> None:
    """One concentration scalar is one molecule.

    AWC answers volatile odorants and ASE water-soluble ones. Driving both from a
    single "food" value makes the ON and OFF cells fight, and measurably inverts
    the downstream answer: AWC alone gives AIY-AIB of +0.148 mV on a rising
    gradient, which is a run, while all six together give -0.057 mV, a turn
    (docs/model_assumptions.md 5J). The default therefore selects a pathway.
    """
    assert set(Chemosensation.build().cells) == {"AWCL", "AWCR"}
    assert len(Chemosensation.build(modality=None).cells) == 6
    assert set(Chemosensation.build(modality="soluble").cells) == {"ASEL", "ASER"}
    assert set(Chemosensation.build(modality="volatile").cells) == {
        "AWAL",
        "AWAR",
        "AWCL",
        "AWCR",
    }


def test_an_unknown_modality_is_refused() -> None:
    with pytest.raises(KeyError, match="modality must be one of"):
        Chemosensation.build(modality="taste")


def test_every_modality_names_cells_that_exist() -> None:
    connectome, _ = load("cook_2019_herm")
    present = {c.id for c in connectome.cells}
    for name, cells in MODALITIES.items():
        assert set(cells) <= present, f"{name} names cells not in the connectome"
    known = {s.cell for s in CHEMOSENSORS}
    for name, cells in MODALITIES.items():
        assert set(cells) <= known, f"{name} names a cell with no response direction"
