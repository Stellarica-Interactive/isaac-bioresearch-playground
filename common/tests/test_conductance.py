"""The conductance-based RMD model.

The test that matters is :func:`test_plateau_is_bistable`. Everything else here
guards the transcription; that one asserts the property the whole model was
imported for, and the property the graded leaky integrator provably lacks.
"""

from __future__ import annotations

import numpy as np
import pytest

from common.neural.conductance import (
    DEFAULT_DT_MS,
    RMD_STATE_NAMES,
    ConductanceModel,
    load_parameters,
    model_provenance,
)

MODEL = "rmd_nicoletti2019"


@pytest.fixture(scope="module")
def model() -> ConductanceModel:
    return ConductanceModel.load(MODEL)


@pytest.fixture(scope="module")
def rest(model: ConductanceModel) -> np.ndarray:
    return model.settle(model.initial_state(1), duration_ms=3000.0)


# -- provenance ----------------------------------------------------------------


def test_parameters_come_from_the_published_source() -> None:
    """Spot-check against RMD.ode. If the importer ever silently changes what it
    extracts, these are the values that would move."""
    p = load_parameters(MODEL)
    assert p["c"] == 1.2  # membrane capacitance, pF
    assert p["gleak"] == 0.4  # nS
    assert p["ek"] == -80.0
    assert p["eca"] == 60.0
    assert p["ena"] == 30.0
    assert p["gshal"] == 2.48
    assert p["vol"] == 5.65  # um^3, from NeuroMorpho
    # Declared as a bare assignment rather than `par` in the source; the importer
    # promotes numeric literals so it is reachable here.
    assert p["t_sk"] == 6.3


def test_provenance_records_the_correction() -> None:
    """A reader must be able to discover that the paper's printed equations are
    wrong without already knowing it."""
    prov = model_provenance(MODEL)
    assert prov["doi"] == "10.1371/journal.pone.0218738"
    assert prov["correction_doi"] == "10.1371/journal.pone.0256930"
    assert prov["sha256"]
    assert prov["file"] == "RMD.ode"


def test_every_state_variable_in_the_source_is_modelled() -> None:
    assert len(RMD_STATE_NAMES) == 22
    assert RMD_STATE_NAMES[0] == "v", "voltage must be first, to match the runtime"
    assert len(set(RMD_STATE_NAMES)) == len(RMD_STATE_NAMES)


# -- resting behaviour ---------------------------------------------------------


def test_rests_near_minus_seventy(rest: np.ndarray) -> None:
    """RMD rests around -70 mV in the animal (Lockery & Goodman 2009), and the
    source's own initial condition is -70. Reaching it from those initials without
    being told to is evidence the channel transcription is right."""
    assert rest[0, 0] == pytest.approx(-69.5, abs=1.5)


def test_resting_calcium_is_near_background(model: ConductanceModel, rest: np.ndarray) -> None:
    ca = rest[model.index("ca_intra1"), 0]
    assert 0.0 < ca < 1e-3


def test_gating_variables_stay_in_range(model: ConductanceModel, rest: np.ndarray) -> None:
    """A gate outside [0, 1] means an equation is wrong, even if the voltage
    trace happens to look plausible."""
    state = model.run(rest.copy(), duration_ms=200.0, i_ext_pa=20.0)
    gates = [i for i, n in enumerate(RMD_STATE_NAMES) if n not in ("v", "ca_intra1")]
    assert state[gates].min() >= -1e-9
    assert state[gates].max() <= 1.0 + 1e-9


# -- the point of the exercise -------------------------------------------------


def test_plateau_is_bistable(model: ConductanceModel, rest: np.ndarray) -> None:
    """Two stable states at zero input, reachable in both directions.

    This is what a plateau potential *is*, and it is exactly what the graded
    leaky integrator cannot produce -- see docs/model_assumptions.md 5C.9, where
    the whole network's sustained variation measured 4e-5 mV.

    A brief depolarising pulse leaves the cell elevated after the current is
    removed, and it stays there. A brief hyperpolarising pulse puts it back.
    Neither final state has any input at all.
    """
    resting_v = rest[0, 0]

    up = model.run(rest.copy(), duration_ms=30.0, i_ext_pa=10.0)
    up = model.run(up, duration_ms=1000.0, i_ext_pa=0.0)
    assert up[0, 0] > resting_v + 15.0, "no plateau: it fell straight back to rest"

    held = model.run(up.copy(), duration_ms=2000.0, i_ext_pa=0.0)
    assert held[0, 0] == pytest.approx(up[0, 0], abs=1.0), "plateau did not hold"

    down = model.run(up.copy(), duration_ms=30.0, i_ext_pa=-15.0)
    down = model.run(down, duration_ms=1000.0, i_ext_pa=0.0)
    assert down[0, 0] == pytest.approx(resting_v, abs=2.0), "could not switch back off"


def test_depolarisation_is_regenerative_not_ohmic(
    model: ConductanceModel, rest: np.ndarray
) -> None:
    """Equal and opposite currents must not give equal and opposite voltages.

    A passive cell answers +5 pA and -5 pA symmetrically. This one does not, and
    the asymmetry is the active conductance doing its job.
    """
    hyper = model.run(rest.copy(), duration_ms=50.0, i_ext_pa=-5.0)[0, 0]
    depol = model.run(rest.copy(), duration_ms=50.0, i_ext_pa=+5.0)[0, 0]
    rise = depol - rest[0, 0]
    fall = rest[0, 0] - hyper
    assert rise > 3.0 * fall, f"response looks passive: +{rise:.1f} mV vs -{fall:.1f} mV"


def test_hyperpolarisation_is_monotonic(model: ConductanceModel, rest: np.ndarray) -> None:
    voltages = [
        model.run(rest.copy(), duration_ms=50.0, i_ext_pa=amp)[0, 0]
        for amp in (-30.0, -15.0, -5.0, 0.0)
    ]
    assert voltages == sorted(voltages)


# -- mechanics -----------------------------------------------------------------


def test_cells_are_independent_when_driven_differently(model: ConductanceModel) -> None:
    """One instance serves every RMD cell, so the vectorisation must not leak
    state between columns."""
    state = model.settle(model.initial_state(3), duration_ms=2000.0)
    driven = model.run(state, duration_ms=100.0, i_ext_pa=np.array([0.0, 0.0, 20.0]))
    assert driven[0, 0] == pytest.approx(driven[0, 1])
    assert driven[0, 2] > driven[0, 0] + 10.0


def test_vectorised_matches_single_cell(model: ConductanceModel) -> None:
    one = model.run(
        model.settle(model.initial_state(1), duration_ms=500.0),
        duration_ms=100.0,
        i_ext_pa=8.0,
    )
    many = model.run(
        model.settle(model.initial_state(4), duration_ms=500.0),
        duration_ms=100.0,
        i_ext_pa=8.0,
    )
    assert np.allclose(many[0], one[0, 0])


def test_is_deterministic(model: ConductanceModel, rest: np.ndarray) -> None:
    a = model.run(rest.copy(), duration_ms=200.0, i_ext_pa=7.0)
    b = model.run(rest.copy(), duration_ms=200.0, i_ext_pa=7.0)
    assert np.array_equal(a, b)


def test_currents_sum_is_finite_across_the_voltage_range(model: ConductanceModel) -> None:
    state = model.initial_state(1)
    for v in (-100.0, -70.0, -40.0, 0.0, 30.0):
        state[0] = v
        total = sum(model.currents(state).values())
        assert np.all(np.isfinite(total))


def test_endpoints_agree_at_every_timestep_and_prove_nothing(
    model: ConductanceModel, rest: np.ndarray
) -> None:
    """The trap at :data:`DEFAULT_DT_MS`, asserted so nobody re-derives it.

    Rest and plateau are both equilibria, so their voltages are essentially
    independent of the timestep. A test that compared only final voltages would
    pass at any dt and would be measuring the fixed point, not the integration.
    """
    ends = [
        model.run(rest.copy(), duration_ms=120.0, i_ext_pa=10.0, dt_ms=dt)[0, 0]
        for dt in (0.01, 0.1, 1.0)
    ]
    assert max(ends) - min(ends) < 0.1, "endpoints should agree -- that is the point"


def test_trajectory_error_is_what_actually_constrains_the_timestep(
    model: ConductanceModel, rest: np.ndarray
) -> None:
    """Measured on the path, the timestep does matter, and monotonically so."""

    def path(dt: float) -> np.ndarray:
        _, trace = model.run_recorded(rest.copy(), duration_ms=60.0, i_ext_pa=10.0, dt_ms=dt)
        trace = np.asarray(trace).ravel()
        grid = np.arange(0.5, 59.0, 0.5)
        return np.interp(grid, np.arange(1, len(trace) + 1) * dt, trace)

    reference = path(0.01)
    errors = [np.abs(path(dt) - reference).max() for dt in (DEFAULT_DT_MS, 1.0, 2.0)]
    assert errors == sorted(errors), "coarser steps must not be more accurate"
    assert errors[0] < 0.5, "the default timestep should track the reference closely"
    assert errors[-1] > 2.0, "a very coarse step should be visibly worse"


def test_recorded_trace_matches_the_plain_run(model: ConductanceModel, rest: np.ndarray) -> None:
    plain = model.run(rest.copy(), duration_ms=20.0, i_ext_pa=10.0)
    final, trace = model.run_recorded(rest.copy(), duration_ms=20.0, i_ext_pa=10.0)
    assert np.array_equal(plain, final)
    assert trace[-1] == pytest.approx(final[0, 0])
