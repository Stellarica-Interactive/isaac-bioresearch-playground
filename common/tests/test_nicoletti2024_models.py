"""The seven Nicoletti 2024 cells, as runnable conductance models.

These are the first models in this project assembled by a tool rather than
written by hand, and the assembly crosses two sources -- channel kinetics from
the ``.mod`` files, conductances from the per-cell drivers. Most of what went
wrong during the import produced a model that loaded and was wrong, so these
check the properties that a wrong one would fail: settling to a real
equilibrium, carrying the channels its own source lists, and leaving the 2019
models untouched.

They do **not** check agreement with the published resting potentials. That
comparison needs the authors' own figures and has not been done -- see
docs/model_assumptions.md 5X.
"""

from __future__ import annotations

import numpy as np
import pytest

from common.neural.conductance import (
    MODEL_CHANNELS,
    ConductanceModel,
    channels_for,
    state_names_for,
)

CELLS_2024 = ("aval", "avar", "rim", "va5", "vd5", "aiy", "vb6")

#: What each cell's own driver assigns, recovered from the clamp modules rather
#: than from the comments above the vectors -- three of which are wrong (5W).
EXPECTED_CHANNELS: dict[str, int] = {
    "aval": 4,
    "avar": 5,
    "rim": 7,
    "va5": 7,
    "vd5": 8,
    "aiy": 7,
    "vb6": 13,
}


def _model(cell: str) -> ConductanceModel:
    return ConductanceModel.load(f"{cell}_nicoletti2024")


@pytest.mark.parametrize("cell", CELLS_2024)
def test_each_cell_loads_with_the_channels_its_source_lists(cell: str) -> None:
    assert len(channels_for(f"{cell}_nicoletti2024")) == EXPECTED_CHANNELS[cell]


@pytest.mark.parametrize("cell", CELLS_2024)
def test_each_cell_settles_to_a_real_equilibrium(cell: str) -> None:
    """Currents summing to zero is not enforced anywhere.

    Nothing in the code makes this true -- it is what a model that genuinely
    relaxed to a fixed point looks like, and the cheapest evidence that the
    kinetics and the conductances were assembled consistently. A model with a
    channel at the wrong density still loads, still runs, and does not balance.
    """
    model = _model(cell)
    rest = model.resting_state()
    assert np.all(np.isfinite(rest))
    total = float(sum(np.atleast_1d(v)[0] for v in model.currents(rest).values()))
    assert abs(total) < 0.05, f"{cell} rests with {total:+.4f} pA unbalanced"


@pytest.mark.parametrize("cell", CELLS_2024)
def test_resting_potential_is_in_a_plausible_range(cell: str) -> None:
    """Deliberately wide. This catches a model that has diverged or inverted a
    sign, not one that is merely a few millivolts off -- which is exactly the
    error this cannot detect, and why the published-value comparison still
    matters."""
    v = float(_model(cell).resting_state()[0][0])
    assert -100.0 < v < 0.0, f"{cell} rests at {v:+.2f} mV"


@pytest.mark.parametrize("cell", CELLS_2024)
def test_every_state_has_kinetics_driving_it(cell: str) -> None:
    """A state nothing writes to would sit at its initial value forever, and the
    model would still run. ``slo1iso`` was imported with the phantom states
    ``FROM`` and ``TO`` at one point, which is exactly that failure."""
    mid = f"{cell}_nicoletti2024"
    model = _model(cell)
    names = state_names_for(mid)
    inf, tau = model._gates(model.resting_state())
    for name in names:
        if name in ("v", "ca_intra1"):
            continue  # integrated directly, not through a gate
        assert name in inf, f"{cell}: no steady state for {name}"
        assert name in tau, f"{cell}: no time constant for {name}"


@pytest.mark.parametrize("cell", ("va5", "vd5", "aiy", "vb6"))
def test_the_2024_calcium_pool_respects_its_floor(cell: str) -> None:
    """cadiff.mod clamps calcium at 1e-4 mM -- 100 nM, as its comment says -- and
    starts there. A pool that drifted below it would mean the clamp was lost in
    translation, and the SLO iso variants read this value."""
    model = _model(cell)
    rest = model.resting_state()
    ca = float(rest[model.index("ca_intra1")][0])
    assert ca >= 1e-4 - 1e-12, f"{cell} calcium fell to {ca:.3e} mM"


@pytest.mark.parametrize("cell", CELLS_2024)
def test_each_cell_responds_monotonically_to_injected_current(cell: str) -> None:
    """A depolarising current must depolarise. Cheap, and it fails for a model
    with an inverted reversal potential or a sign error in a current term."""
    model = _model(cell)
    rest = model.resting_state()
    voltages = [
        float(model.run(rest.copy(), duration_ms=2000.0, i_ext_pa=i)[0][0])
        for i in (-5.0, 0.0, 5.0)
    ]
    assert voltages[0] < voltages[1] < voltages[2], f"{cell}: {voltages}"


def test_vb6_is_the_cell_that_was_called_unrecoverable() -> None:
    """5G.4 named VB6 the most wanted datum in the project, 5I concluded its
    channel mapping could not be recovered, and an email was sent to the authors
    asking for it. The mapping was in the clamp driver the whole time (5W)."""
    channels = channels_for("vb6_nicoletti2024")
    assert len(channels) == 13
    # All six SLO variants, which is what made it the largest of the seven.
    assert sum(c.startswith("slo") for c in channels) == 6
    assert float(_model("vb6").resting_state()[0][0]) < 0.0


@pytest.mark.parametrize("cell", CELLS_2024)
def test_generated_models_are_not_in_the_hand_written_table(cell: str) -> None:
    """They declare their own channels, so MODEL_CHANNELS needs no entry -- a
    step that could be forgotten, and whose absence would read as an unknown
    model rather than a missing line."""
    assert f"{cell}_nicoletti2024" not in MODEL_CHANNELS


def test_the_2019_models_are_unaffected() -> None:
    """Three rounds of refactoring went into making a subset of channels
    loadable. These two values are the regression test for all of it."""
    assert float(ConductanceModel.load("rmd_nicoletti2019").resting_state()[0][0]) == (
        pytest.approx(-69.49, abs=0.01)
    )
    assert float(ConductanceModel.load("awc_nicoletti2019").resting_state()[0][0]) == (
        pytest.approx(-69.18, abs=0.01)
    )


#: Resting potentials given as model outputs in Nicoletti et al. 2024
#: (PLOS ONE 19(3):e0298105, CC-BY). The paper states no value for AVAR, AIY or
#: RIM, so three of the seven cannot be checked this way -- and for AVAR it notes
#: a discrepancy between current-clamp and voltage-clamp in the experimental data
#: itself.
PUBLISHED_RESTING_MV: dict[str, float] = {
    "va5": -75.20,
    "vb6": -53.19,
    "vd5": -44.61,
    "aval": -25.40,
}


@pytest.mark.parametrize("cell", ("va5", "vb6", "vd5"))
def test_resting_potential_matches_the_paper(cell: str) -> None:
    """The check that makes this an import rather than a plausible imitation.

    Any one of the name normalisations going wrong -- SLO-2's vestigial ``1``
    suffixes, ``bkg`` against ``fondo``, ``shift`` standing in for four of our
    names, ``pi = 3.14`` rather than the real constant, the scale factors written
    as bare literals inside formulas -- would move these by more than the
    tolerance. VB6 in particular carries all six SLO variants and 323 parameters.

    2 mV, because the paper quotes two decimals and the drivers read their value
    at t = 50-60 ms while ours is the asymptote. AVAL is excluded; see
    :func:`test_aval_does_not_match_the_paper`.
    """
    observed = float(_model(cell).resting_state()[0][0])
    assert observed == pytest.approx(PUBLISHED_RESTING_MV[cell], abs=2.0)


@pytest.mark.xfail(
    reason="AVAL rests at -39.4 mV against the paper's -25.4, a 14 mV gap with no "
    "diagnosis. It is the simplest of the seven -- four channels, NCA at zero, a "
    "comment that agrees with its driver -- and shares every constant with cells "
    "that do reproduce. Reaching -25.4 needs about +1.7 pA of standing inward "
    "current and its IClamp fires at 1023 ms, long after the measurement window. "
    "Recorded as a known failure rather than fitted away: see 5AD.3.",
    strict=True,
)
def test_aval_does_not_match_the_paper() -> None:
    observed = float(_model("aval").resting_state()[0][0])
    assert observed == pytest.approx(PUBLISHED_RESTING_MV["aval"], abs=2.0)


@pytest.mark.parametrize(
    "model_id",
    [f"{c}_nicoletti2024" for c in CELLS_2024] + ["rmd_nicoletti2019", "awc_nicoletti2019"],
)
def test_derivatives_is_the_field_step_integrates(model_id: str) -> None:
    """``derivatives`` exists to linearise a cell, so it must be the same model.

    Off equilibrium -- 8 mV above rest, under 3 pA and 0.5 nS of external
    input -- so every term is non-zero. As the step shrinks, ``step`` becomes
    one Euler step of the true field (the gates' exact exponential update
    differs from Euler by ``h / 2 tau``), so the two must agree.
    """
    model = ConductanceModel.load(model_id)
    state = model.resting_state(1).copy()
    state[0] += 8.0
    h = 1.0e-6
    stepped = (model.step(state, h, 3.0, 0.5) - state) / h
    field = model.derivatives(state, 3.0, 0.5)
    assert np.allclose(field, stepped, rtol=1e-3, atol=1e-7)
