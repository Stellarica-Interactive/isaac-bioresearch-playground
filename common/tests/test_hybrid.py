"""Running a few named cells on their own measured biophysics.

The load-bearing test is :func:`test_other_cells_are_untouched`. Promoting six
cells out of 397 must not perturb the other 391, or every earlier result becomes
incomparable with every later one.
"""

from __future__ import annotations

import numpy as np
import pytest

from common.data.schemas import CellCategory
from common.neural.hybrid import INTRINSIC_MODELS, HybridRuntime
from common.neural.runtime import NeuralRuntime
from common.neural.synapses import UnknownSignPolicy
from worm.importers.naming import body_wall_muscle_ids
from worm.loader import load
from worm.neural.config import RUNTIME_OVERLAYS, build_runtime

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def _runtime() -> NeuralRuntime:
    connectome, _ = load("cook_2019_herm", annotations=RUNTIME_OVERLAYS)
    muscles = body_wall_muscle_ids()
    cells = tuple(
        c.id for c in connectome.cells if c.category is CellCategory.NEURON or c.id in muscles
    )
    runtime, _ = build_runtime(
        "cook_2019_herm",
        unknown_sign=UnknownSignPolicy.EXCLUDE,
        cells=cells,
        connectome=connectome,
        dt_ms=1.0,
    )
    return runtime


@pytest.fixture(scope="module")
def hybrid() -> HybridRuntime:
    return HybridRuntime.build(_runtime())


# -- what it is --------------------------------------------------------------


def test_promotes_only_the_cells_with_a_published_model(hybrid: HybridRuntime) -> None:
    """Nicoletti characterised RMD. Nobody characterised the other 391 cells, and
    extending a channel complement to them would be the largest silent assumption
    in the project."""
    assert set(hybrid.cells) == set(INTRINSIC_MODELS["rmd_nicoletti2019"])
    assert len(hybrid.cells) == 6


def test_intrinsic_cells_rest_where_the_measurement_says(hybrid: HybridRuntime) -> None:
    """Before coupling, they must sit at the conductance model's own resting
    potential rather than the graded model's."""
    fresh = HybridRuntime.build(_runtime())
    assert all(v == pytest.approx(-69.5, abs=1.5) for v in fresh.voltages().values())


def test_the_network_sees_the_measured_resting_potential(hybrid: HybridRuntime) -> None:
    fresh = HybridRuntime.build(_runtime())
    for cell in fresh.cells:
        assert fresh.runtime.voltage(cell) == pytest.approx(-69.5, abs=1.5)


# -- the test that matters ---------------------------------------------------


def test_other_cells_are_untouched() -> None:
    """Promoting six cells must not disturb the other 391.

    Otherwise every result measured before this change becomes incomparable with
    every result measured after it, and there would be no way to attribute a
    difference to the intrinsic dynamics rather than to a side effect.

    RMD's own partners are excluded from the comparison: they are *supposed* to
    change, because RMD now rests 48 mV away from where the graded model put it.
    """
    plain = _runtime()
    hybrid = HybridRuntime.build(_runtime())

    plain.run(200.0)
    hybrid.run(200.0)

    network = hybrid.runtime.network
    promoted = set(hybrid.cells)
    # Anything wired directly to an RMD cell legitimately differs.
    coupled = set(promoted)
    for cell in promoted:
        i = network.index(cell)
        touching = np.flatnonzero(
            (network.g_gap[i] > 0) | (network.g_syn[i] > 0) | (network.g_syn[:, i] > 0)
        )
        coupled |= {network.cell_ids[int(j)] for j in touching}

    far = [c for c in network.cell_ids if c not in coupled]
    assert len(far) > 200, "too few uncoupled cells left to make this test meaningful"

    idx = np.array([network.index(c) for c in far])
    difference = np.abs(plain.state[0][idx] - hybrid.runtime.state[0][idx]).max()
    assert difference < 5.0, (
        f"cells not connected to RMD moved by {difference:.2f} mV; promoting six "
        "cells should not reach them this strongly in 200 ms"
    )


# -- coupling ----------------------------------------------------------------


def test_network_drive_returns_a_conductance_not_just_a_current(
    hybrid: HybridRuntime,
) -> None:
    """The restoring term. Handing the conductance model only a frozen current
    removes the part of gap coupling that opposes the cell's own movement, and RMD
    carries up to 3.8 nS against 1.2 pF -- it diverged within 9 ms."""
    current, conductance = hybrid.network_drive()
    assert current.shape == (len(hybrid.cells),)
    assert conductance.shape == (len(hybrid.cells),)
    assert np.all(conductance > 0.0), "RMD is gap-coupled; its conductance cannot be 0"


def test_stays_finite_embedded_in_the_network() -> None:
    """The regression that motivated the conductance term."""
    hybrid = HybridRuntime.build(_runtime())
    hybrid.run(500.0)
    assert np.all(np.isfinite(hybrid.state))
    assert np.all(np.isfinite(hybrid.runtime.state))
    assert np.abs(hybrid.state[0]).max() < 200.0


def test_a_missing_cell_set_is_refused() -> None:
    with pytest.raises(ValueError, match="none of"):
        HybridRuntime.build(_runtime(), cells=("NOT_A_CELL",))


def test_substepping_is_finer_than_the_network_step(hybrid: HybridRuntime) -> None:
    """The conductance model needs dt <= 0.05 ms; the network runs at 1 ms."""
    assert hybrid.substep_ms <= 0.1
    assert hybrid.runtime.dt_ms > hybrid.substep_ms
