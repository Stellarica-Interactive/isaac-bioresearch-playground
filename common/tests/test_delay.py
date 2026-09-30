"""Conduction delay on chemical transmission.

The delay experiment exists because c302's travelling wave is produced by a
per-connection delay ladder (docs/model_assumptions.md §5Z), and the honest
version of that is one uniform number. These tests guard the two things that
would quietly make it meaningless: delaying the wrong term, and a zero delay
that is not actually a no-op.
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
from common.neural.delay import DEFAULT_DELAY_MS, ConductionDelay
from common.neural.runtime import Integrator, NeuralRuntime
from common.neural.synapses import UnknownSignPolicy

PROV = Provenance(
    source_id="toy",
    kind="connectome",
    citation="synthetic fixture, not a real animal",
    url="",
    license="n/a",
)


def _runtime(
    *, electrical: bool = False, integrator: Integrator = Integrator.EXPONENTIAL
) -> NeuralRuntime:
    cells = tuple(Cell(id=i, category=CellCategory.NEURON) for i in ("A", "B"))
    kind = SynapseType.ELECTRICAL if electrical else SynapseType.CHEMICAL
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
                synapse_type=kind,
                weight=50,
                sign=None if electrical else Sign.EXCITATORY,
            ),
        ),
        provenance=PROV,
    )
    return NeuralRuntime.build(
        connectome, unknown_sign=UnknownSignPolicy.EXCLUDE, integrator=integrator
    )


def test_zero_delay_is_a_no_op() -> None:
    """Otherwise the control condition is not a control.

    Every comparison in the delay experiment is against no delay, so a zero
    delay that perturbed anything would make the whole sweep unreadable.
    """
    plain = _runtime()
    plain.run(500.0)

    delayed = _runtime()
    delayed.s_pre = ConductionDelay.build(delayed, delay_ms=0.0)
    delayed.run(500.0)

    assert delayed.s_pre.steps == 0
    assert np.allclose(plain.state, delayed.state)


def test_the_buffer_starts_primed_not_empty() -> None:
    """A buffer of zeros silences every synapse for the first delay_ms.

    That is a transient nobody asked for, and in a model whose first twenty
    seconds are already transient (§5Q.2) it would be indistinguishable from a
    result.
    """
    runtime = _runtime()
    runtime.run(200.0)
    resting = runtime.state[1].copy()
    delay = ConductionDelay.build(runtime, delay_ms=20.0)
    assert np.allclose(delay.delayed, resting)
    assert not np.allclose(delay.delayed, 0.0)


def test_the_delay_is_quantised_to_the_timestep_and_says_so() -> None:
    runtime = _runtime()
    delay = ConductionDelay.build(runtime, delay_ms=7.3)
    assert delay.quantised_delay_ms == pytest.approx(delay.steps * runtime.dt_ms, abs=1e-12)
    assert abs(delay.quantised_delay_ms - 7.3) <= runtime.dt_ms


def test_activation_arrives_late_by_the_requested_amount() -> None:
    runtime = _runtime()
    delay = ConductionDelay.build(runtime, delay_ms=10.0)
    runtime.s_pre = delay

    marks = []
    for step in range(delay.steps * 3):
        runtime.state[1, :] = float(step)
        delay.advance(runtime.state[1])
        marks.append(float(delay.delayed[0]))

    # What arrives now is what was recorded `steps` advances ago.
    assert marks[delay.steps * 2] == pytest.approx(delay.steps * 2 - delay.steps)


@pytest.mark.parametrize("integrator", [Integrator.EXPONENTIAL, Integrator.RK4, Integrator.EULER])
def test_a_delay_changes_a_chemical_network(integrator: Integrator) -> None:
    """The mechanism has to do something, or the null result is about nothing.

    Parametrised over every integrator because the first version of this test
    used only the default RK4, while the project runs on EXPONENTIAL -- which
    bypasses ``_derivs`` and computes its own currents. The delay reached the
    path the test exercised and not the path the simulator uses, so four delays
    from 0 to 50 ms gave byte-identical trajectories and the test still passed.
    """
    fresh = _runtime(integrator=integrator)
    fresh.i_ext_pa[fresh.network.index("A")] = 15.0
    fresh.run(30.0)

    delayed = _runtime(integrator=integrator)
    delayed.s_pre = ConductionDelay.build(delayed, delay_ms=50.0)
    delayed.i_ext_pa[delayed.network.index("A")] = 15.0
    delayed.run(30.0)

    post = delayed.network.index("B")
    assert not np.isclose(fresh.state[0, post], delayed.state[0, post])


def test_gap_junctions_are_not_delayed() -> None:
    """An electrical synapse is a resistive pore: it has no transmission delay.

    c302's ``DelayedGapJunction`` is not modelling conduction, it is injecting a
    phase gradient (§5Z). Delaying gap junctions here would make the experiment
    easier to fit and impossible to defend, so a purely electrical network must
    be completely unaffected.
    """
    plain = _runtime(electrical=True)
    plain.i_ext_pa[plain.network.index("A")] = 15.0
    plain.run(100.0)

    delayed = _runtime(electrical=True)
    delayed.s_pre = ConductionDelay.build(delayed, delay_ms=50.0)
    delayed.i_ext_pa[delayed.network.index("A")] = 15.0
    delayed.run(100.0)

    assert np.allclose(plain.state, delayed.state)


def test_a_negative_delay_is_refused() -> None:
    with pytest.raises(ValueError, match="must not be negative"):
        ConductionDelay(n_cells=2, dt_ms=1.0, delay_ms=-1.0)


def test_the_default_is_documented_as_assumed() -> None:
    """Nothing in this project measures a synaptic delay for C. elegans, and
    c302's per-connection values run to 1500 ms -- three orders larger than any
    synaptic delay, which is what makes them a phase gradient rather than one."""
    assert 0.0 < DEFAULT_DELAY_MS < 50.0
