"""The A-type motor neurons as slow oscillators (model_assumptions 5AO.18).

Phenomenological, so what is tested is that it reproduces the two numbers it was
chosen to: an isolated cell bursts every 50 to 90 s for 2 to 3 s (Gao et al.
2018), and that it runs inside the network the way any intrinsic model does.
"""

from __future__ import annotations

import argparse

import numpy as np
import pytest

from worm.neural.a_oscillator import A_TYPE_CELLS, AOscillator, burst_statistics


def test_an_isolated_cell_bursts_as_gao_et_al_measured() -> None:
    model = AOscillator()
    trace = model.run(model.resting_state(1), 300_000.0, dt_ms=1.0)[:, 0]
    interval, length = burst_statistics(trace, 1.0)
    assert 50.0 <= interval <= 90.0
    assert 2.0 <= length <= 3.0
    assert trace.max() > -31.0  # a burst reaches past the release threshold


def test_it_is_quiet_most_of_the_time() -> None:
    """Held down between bursts: what lets a command decide when it runs."""
    model = AOscillator()
    trace = model.run(model.resting_state(1), 300_000.0, dt_ms=1.0)[:, 0]
    assert np.mean(trace < -60.0) > 0.9


def test_the_cells_are_the_twenty_one_a_type_motor_neurons() -> None:
    assert len(A_TYPE_CELLS) == 21
    assert all(c[:2] in ("DA", "VA") for c in A_TYPE_CELLS)


def test_the_loop_runs_the_a_cells_on_it() -> None:
    from tools.analyse_loop import Loop

    loop = Loop(
        argparse.Namespace(
            unknown_sign="exclude",
            physics_hz=240.0,
            neural_dt_ms=1.0,
            torque_scale=1e-3,
            proprioceptive_mv=400.0,
            proprio_rate=0.7,
            command_mv=20.0,
            settle_s=0.0,
            open_loop=False,
            param=[],
            intrinsic=AOscillator.model_id,
        )
    )
    assert loop.hybrid is not None
    assert set(loop.hybrid.cells) == set(A_TYPE_CELLS)
    for _ in range(24):
        loop.step()
    v = loop.runtime.state[0, loop.hybrid.indices]
    np.testing.assert_array_equal(v, loop.hybrid.state[0])


def test_a_model_object_needs_its_cells_named() -> None:
    from common.neural.hybrid import HybridRuntime
    from worm.neural.config import build_runtime

    runtime, _ = build_runtime("cook_2019_herm", unknown_sign="exclude")
    with pytest.raises(ValueError):
        HybridRuntime.build(runtime, model=AOscillator())
