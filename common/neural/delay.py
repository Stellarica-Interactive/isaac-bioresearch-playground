"""Conduction delay on synaptic transmission.

    delay = ConductionDelay.build(runtime, delay_ms=50.0)
    runtime.s_pre = delay           # the runtime reads presynaptic activation here

Where it came from
------------------

Every hypothesis tested in ``docs/negative_result.md`` is a mechanism that would
let a population *hold* or *express* a phase difference, and all seven were
rejected. c302, OpenWorm's neural model, names a mechanism none of them covers:
its coupling carries **delays**. From ``c302/parameters_C2.py``, AVB to the
B-type motor neurons:

    AVBR_to_DB1_elec_syn_delay      0 ms
    AVBL_to_DB2_elec_syn_delay    250 ms
    ...
    AVBL_to_DB7_elec_syn_delay   1500 ms

A 250 ms per segment ladder from the command interneuron **is** a travelling
wave: the phase gradient is in the parameter file, not computed by the circuit.
Importing it would produce a gait belonging to those eighteen numbers, and
reporting it as the connectome's would be dishonest -- see §5Z.

So this is deliberately **one** number, uniform across every connection, rather
than a per-connection gradient. That makes it a statement about the mechanism:
if a uniform delay produces propagation, delay is sufficient; if it does not,
then c302's wave depends on its gradient being graded, which is a sharper result
than anything in §5T.

Is it biology?
--------------

Partly, and the parts should not be confused.

**Chemical transmission delay is real.** Vesicle fusion, diffusion across the
cleft and receptor binding take time -- of order a millisecond at a typical
synapse, and this project has no measurement of it for *C. elegans*.

**Gap-junction delay is not.** An electrical synapse is a resistive pore; it has
no transmission delay, only the RC filtering the membrane already provides.
c302's ``DelayedGapJunction`` is not modelling conduction, it is injecting a
phase gradient. This module therefore delays **chemical transmission only**, and
leaves the gap-junction term instantaneous. A flag to delay gap junctions too
would make the experiment easier to fit and impossible to defend.

:data:`DEFAULT_DELAY_MS` is ASSUMED. Nothing measures it here, which is the whole
reason it is swept rather than set.

How it works
------------

A ring buffer of past activation vectors. The postsynaptic current uses
activation from ``delay_ms`` ago; each cell's own ``ds/dt`` uses its present
value, because the delay is in transmission and not in the cell's own kinetics.

The buffer advances once per integration step rather than once per Runge-Kutta
stage, so the delayed value is held constant across the stages of a step. That
is the standard treatment of a delay under a fixed-step integrator and it means
the effective delay is quantised to the timestep -- reported by
:attr:`quantised_delay_ms`, which is what to quote rather than the requested
value.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

#: Delay applied to chemical transmission, ms. ASSUMED -- there is no measurement
#: of it for *C. elegans* in this project. Chemical synapses do have a real delay
#: of roughly this order (vesicle fusion, cleft diffusion, receptor binding), but
#: the value here is chosen to be swept, not believed. c302 uses per-connection
#: delays up to 1500 ms, which is three orders larger than any synaptic delay and
#: is a phase gradient wearing a delay's name.
DEFAULT_DELAY_MS = 5.0


@dataclass
class ConductionDelay:
    """Presynaptic activation, as it was ``delay_ms`` ago."""

    n_cells: int
    dt_ms: float
    delay_ms: float = DEFAULT_DELAY_MS
    _buffer: np.ndarray = field(init=False, repr=False)
    _cursor: int = field(init=False, default=0, repr=False)

    def __post_init__(self) -> None:
        if self.delay_ms < 0.0:
            raise ValueError("delay_ms must not be negative")
        if self.dt_ms <= 0.0:
            raise ValueError("dt_ms must be positive")
        # One slot per step of delay, plus the slot being written. A zero delay
        # therefore still allocates one slot and returns the present value, so
        # the zero case is a genuine no-op rather than a special path.
        self._steps = int(round(self.delay_ms / self.dt_ms))
        self._buffer = np.zeros((self._steps + 1, self.n_cells), dtype=np.float64)

    @classmethod
    def build(cls, runtime: object, *, delay_ms: float = DEFAULT_DELAY_MS) -> ConductionDelay:
        """Sized from the runtime, and primed with its current activation.

        Primed rather than zeroed: a buffer of zeros would mean every synapse in
        the network falling silent for the first ``delay_ms``, which is a
        transient nobody asked for and which looks like a result.
        """
        state = np.asarray(runtime.state)  # type: ignore[attr-defined]
        delay = cls(
            n_cells=state.shape[1],
            dt_ms=float(runtime.dt_ms),  # type: ignore[attr-defined]
            delay_ms=delay_ms,
        )
        delay._buffer[:] = state[1]
        return delay

    @property
    def steps(self) -> int:
        """Buffer depth. Zero means no delay."""
        return self._steps

    @property
    def quantised_delay_ms(self) -> float:
        """The delay actually applied, which is the requested one rounded to a
        whole timestep. Quote this, not what was asked for."""
        return self._steps * self.dt_ms

    def advance(self, s: np.ndarray) -> None:
        """Record the present activation. Call once per integration step."""
        self._cursor = (self._cursor + 1) % self._buffer.shape[0]
        self._buffer[self._cursor] = s

    @property
    def delayed(self) -> np.ndarray:
        """Activation as of ``quantised_delay_ms`` ago."""
        oldest = (self._cursor + 1) % self._buffer.shape[0]
        return self._buffer[oldest]
