"""Bistable output for named cells.

These tests guard the mechanism, not a result. The module exists to test a
hypothesis (docs/model_assumptions.md 5C.9, 5M, 5P) and must never become a
default, so the test that matters most is
:func:`test_a_cell_holds_its_state_inside_the_dead_band` -- holding state is the
entire point, and a Schmitt trigger with no memory is just a threshold.
"""

from __future__ import annotations

import numpy as np
import pytest

from common.data.schemas import (
    Cell,
    CellCategory,
    Connection,
    Connectome,
    Provenance,
    Scope,
    Sign,
    SynapseType,
)
from common.neural.hysteresis import DEFAULT_WIDTH_MV, Hysteresis
from common.neural.runtime import NeuralRuntime
from common.neural.synapses import UnknownSignPolicy

PROV = Provenance(
    source_id="toy",
    kind="connectome",
    citation="synthetic fixture, not a real animal",
    url="",
    license="n/a",
)


def _runtime() -> NeuralRuntime:
    cells = tuple(Cell(id=i, category=CellCategory.NEURON) for i in ("A", "B", "C"))
    connectome = Connectome(
        id="toy",
        organism="C. elegans",
        sex="hermaphrodite",
        stage="adult",
        scope=Scope.WHOLE_ANIMAL,
        cells=cells,
        connections=(
            Connection(
                pre="A",
                post="B",
                synapse_type=SynapseType.CHEMICAL,
                weight=3,
                sign=Sign.EXCITATORY,
            ),
        ),
        provenance=PROV,
    )
    return NeuralRuntime.build(connectome, unknown_sign=UnknownSignPolicy.EXCLUDE)


def test_a_cell_holds_its_state_inside_the_dead_band() -> None:
    """The whole point. Without memory this is a threshold, not a latch, and it
    could not hold a phase difference between neighbouring segments -- which is
    the property 5P says the model is missing.
    """
    runtime = _runtime()
    switch = Hysteresis.build(runtime, ("A",), width_mv=10.0, binary_output=True)
    centre = float(switch.on_mv[0] - 5.0)

    runtime.state[0, switch.indices] = centre + 20.0  # well above the on threshold
    switch.apply(runtime)
    assert runtime.state[1, switch.indices][0] == 1.0

    runtime.state[0, switch.indices] = centre  # back into the dead band
    switch.apply(runtime)
    assert runtime.state[1, switch.indices][0] == 1.0, "latched high must stay high"

    runtime.state[0, switch.indices] = centre - 20.0  # below the off threshold
    switch.apply(runtime)
    assert runtime.state[1, switch.indices][0] == 0.0

    runtime.state[0, switch.indices] = centre  # dead band again
    switch.apply(runtime)
    assert runtime.state[1, switch.indices][0] == 0.0, "latched low must stay low"


def test_output_takes_exactly_two_values() -> None:
    """Two-valued, but not 0 and 1 by default -- see
    :func:`test_the_two_output_levels_match_the_graded_model_they_replace`."""
    runtime = _runtime()
    switch = Hysteresis.build(runtime, ("A", "B"), width_mv=4.0)
    for offset in (-30.0, -1.0, 0.0, 1.0, 30.0):
        runtime.state[0, switch.indices] = switch.on_mv + offset
        switch.apply(runtime)
        for i, j in enumerate(switch.indices):
            assert runtime.state[1, j] in (switch.off_output[i], switch.on_output[i])


def test_the_two_output_levels_match_the_graded_model_they_replace() -> None:
    """Otherwise this is a gain experiment wearing a bistability label.

    The graded activation of these cells never leaves the bottom tail of its
    sigmoid, so emitting Boyle's literal 0 and 1 would multiply the muscle drive
    several-fold at the same moment as it introduced memory, and neither effect
    could be separated from the other afterwards.
    """
    runtime = _runtime()
    runtime.run(200.0)
    switch = Hysteresis.build(runtime, ("A", "B"), width_mv=6.0)

    assert np.all(switch.on_output > switch.off_output)
    assert np.all(switch.on_output < 1.0), "must not reach the binary scale"

    # The level emitted on latching high is what the graded model would settle at
    # if the cell were held at the upper switch point. Voltage-clamped, because an
    # unclamped cell relaxes off the switch point and settles somewhere else --
    # which is the whole reason the latch has to state its output explicitly.
    # 20_000 steps at dt 0.1 ms is 2 s, ten activation time constants (a_d = 5e-3
    # per ms). The error plateaus at 3.4e-6 there and does not improve with ten
    # times more, against a tolerance of 1e-3 -- the original 200_000 cost 50 s of
    # suite time to buy nothing.
    for _ in range(20_000):
        runtime.state[0, switch.indices] = switch.on_mv
        runtime.step()
    assert np.allclose(runtime.state[1, switch.indices], switch.on_output, atol=1e-3)


def test_binary_output_is_available_for_comparison() -> None:
    runtime = _runtime()
    switch = Hysteresis.build(runtime, ("A", "B"), width_mv=4.0, binary_output=True)
    assert np.array_equal(switch.on_output, np.ones(2))
    assert np.array_equal(switch.off_output, np.zeros(2))


def test_the_band_sits_where_the_cells_currently_are() -> None:
    """Not on their activation threshold, which is where it was first put.

    That threshold is v_threshold_offset_mv above rest -- 30 mV since 5N -- and
    the B-type neurons operate far below it, so nothing ever crossed and the
    population was silenced rather than made bistable. Two dead-band widths then
    produced byte-identical runs, which is what exposed it.
    """
    runtime = _runtime()
    runtime.run(200.0)
    switch = Hysteresis.build(runtime, ("A", "B"), width_mv=6.0)
    assert np.allclose(switch.on_mv - switch.off_mv, 6.0)
    assert np.allclose(0.5 * (switch.on_mv + switch.off_mv), runtime.state[0, switch.indices])

    threshold = np.asarray(runtime.model.v_threshold_mv)[switch.indices]  # type: ignore[attr-defined]
    legacy = Hysteresis.build(runtime, ("A", "B"), width_mv=6.0, centre_on_threshold=True)
    assert np.allclose(0.5 * (legacy.on_mv + legacy.off_mv), threshold)


def test_cells_latch_independently() -> None:
    """Two cells in different states at the same moment is exactly what a
    travelling wave requires and what this model could not previously produce."""
    runtime = _runtime()
    switch = Hysteresis.build(runtime, ("A", "B"), width_mv=4.0)
    runtime.state[0, switch.indices] = switch.on_mv + np.array([10.0, -10.0])
    switch.apply(runtime)
    assert switch.latched[0] and not switch.latched[1]
    assert runtime.state[1, switch.indices[0]] == switch.on_output[0]
    assert runtime.state[1, switch.indices[1]] == switch.off_output[1]
    assert switch.fraction_on == pytest.approx(0.5)


def test_starts_latched_low() -> None:
    runtime = _runtime()
    switch = Hysteresis.build(runtime, ("A", "B"))
    assert switch.fraction_on == 0.0


def test_only_the_named_cells_are_touched() -> None:
    runtime = _runtime()
    switch = Hysteresis.build(runtime, ("A",), width_mv=4.0)
    others = [i for i in range(runtime.network.n) if i not in switch.indices]
    before = runtime.state[1, others].copy()
    runtime.state[0, switch.indices] = switch.on_mv + 50.0
    switch.apply(runtime)
    assert np.array_equal(runtime.state[1, others], before)


def test_an_absent_cell_set_is_refused() -> None:
    with pytest.raises(ValueError, match="none of"):
        Hysteresis.build(_runtime(), ("NOT_A_CELL",))


def test_default_width_is_documented_as_assumed() -> None:
    """Boyle et al.'s epsilon is on a unit input scale with no principled
    translation into millivolts, so a null result has to be swept before it means
    anything."""
    assert DEFAULT_WIDTH_MV > 0.0
