"""A-type motor neurons as slow intrinsic oscillators, after Gao et al. 2018.

    from common.neural.hybrid import HybridRuntime
    from worm.neural.a_oscillator import AOscillator, A_TYPE_CELLS

    hybrid = HybridRuntime.build(runtime, model=AOscillator(), cells=A_TYPE_CELLS)

What is measured, and what is not
---------------------------------

Gao et al. 2018 showed that the A-type motor neurons, which drive backward
locomotion, are rhythm generators in their own right:

* Isolated from the premotor interneurons, they produce bursts on their own: in
  dissected preparations about every 90 s, each burst lasting 2 to 3 s; in glued
  animals DA9's calcium cycles about every 50 s.
* In crawling animals the same cells cycle about five times faster, which the
  authors attribute to proprioceptive feedback.
* The P/Q/N-type calcium channel UNC-2 underlies it: loss of function slows and
  weakens the oscillation, gain of function speeds and strengthens it.
* AVA potentiates them, mainly through its chemical synapses, and the premotor
  interneurons otherwise attenuate their rhythmic activity.

None of their recordings gives the conductances, kinetics or voltages a model
needs, and the one measured A-type cell the project holds (Nicoletti et al.
2024's VA5) has no UNC-2 and no oscillation at all (model_assumptions 5AN.10,
5W). So this is **phenomenological**: the smallest relaxation oscillator with the
measured shape --

* a fast regenerative inward current standing in for UNC-2, ``g_ca * m_inf(V)``,
  which starts each burst;
* a slow recovery current, ``g_k * w``, that builds during a burst and ends it, and
  decays slowly between bursts, so the cell is held down most of the time;
* the graded model's own membrane: 1 pF, 10 pS of leak at -70 mV.

Every one of its parameters is ASSUMED, chosen so that an isolated cell bursts
every 50 to 90 s for 2 to 3 s, Gao et al.'s two measured numbers. What it is for
is model_assumptions 5AO.18: whether A cells that are held down between bursts let
the command interneurons choose the direction of crawling.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

#: The A-type motor neurons of the ventral cord.
A_TYPE_CELLS: tuple[str, ...] = tuple(
    [f"DA{i}" for i in range(1, 10)] + [f"VA{i}" for i in range(1, 13)]
)


@dataclass(frozen=True)
class AOscillatorParameters:
    """All ASSUMED; see the module docstring for the two measured targets."""

    c_pf: float = 1.0
    g_leak_ns: float = 0.01
    e_leak_mv: float = -70.0
    g_ca_ns: float = 3.0
    """Regenerative inward conductance, the UNC-2 stand-in."""
    v_ca_mv: float = -55.0
    k_ca_mv: float = 4.0
    e_ca_mv: float = 50.0
    g_k_ns: float = 10.0
    """Slow recovery conductance: builds during a burst and ends it."""
    e_k_mv: float = -80.0
    v_w_mv: float = -60.0
    """Recovery is fully on across the whole burst, so it always ends one: with its
    midpoint at the burst voltage instead, a cell settles into a stable middle
    state and never bursts again."""
    k_w_mv: float = 3.0
    tau_up_ms: float = 1400.0
    """How fast recovery builds while the cell is up: sets the burst length."""
    tau_down_ms: float = 23000.0
    """How slowly it decays while the cell is down: sets the interval. With these
    values an isolated cell bursts every 66 s for 2.2 s."""


@dataclass
class AOscillator:
    """A model with the interface :class:`~common.neural.hybrid.HybridRuntime` drives.

    State per cell: ``[V (mV), w]``. ``step`` takes the network's drive in the same
    split form the hybrid hands every intrinsic model, ``i_ext - g_ext * V``, and
    integrates the membrane exponentially with the conductances held across the
    step, as the graded runtime does, so it is stable at the hybrid's substep.
    """

    params: AOscillatorParameters = field(default_factory=AOscillatorParameters)
    model_id: str = "a_oscillator_gao2018"

    def m_inf(self, v: np.ndarray) -> np.ndarray:
        p = self.params
        return 1.0 / (1.0 + np.exp(-(v - p.v_ca_mv) / p.k_ca_mv))

    def w_inf(self, v: np.ndarray) -> np.ndarray:
        p = self.params
        return 1.0 / (1.0 + np.exp(-(v - p.v_w_mv) / p.k_w_mv))

    def tau_w(self, v: np.ndarray) -> np.ndarray:
        p = self.params
        return p.tau_down_ms + (p.tau_up_ms - p.tau_down_ms) * self.w_inf(v)

    def resting_state(self, n_cells: int = 1, **_: float) -> np.ndarray:
        """The quiet phase: at the leak reversal, recovery partly built.

        An oscillator has no rest; this is where a cell sits between bursts, so a
        run starts in the long quiet phase rather than in a burst.
        """
        return np.vstack([np.full(n_cells, self.params.e_leak_mv), np.full(n_cells, 0.5)])

    def step(
        self,
        state: np.ndarray,
        dt_ms: float,
        i_ext_pa: np.ndarray | float = 0.0,
        g_ext_ns: np.ndarray | float = 0.0,
    ) -> np.ndarray:
        p = self.params
        v, w = state[0], state[1]
        m = self.m_inf(v)
        g_tot = p.g_leak_ns + p.g_ca_ns * m + p.g_k_ns * w + np.asarray(g_ext_ns)
        drive = (
            p.g_leak_ns * p.e_leak_mv
            + p.g_ca_ns * m * p.e_ca_mv
            + p.g_k_ns * w * p.e_k_mv
            + np.asarray(i_ext_pa)
        )
        v_inf = drive / g_tot
        v_next = v_inf + (v - v_inf) * np.exp(-dt_ms * g_tot / p.c_pf)
        tau = self.tau_w(v)
        w_next = self.w_inf(v) + (w - self.w_inf(v)) * np.exp(-dt_ms / tau)
        return np.vstack([v_next, w_next])

    def run(self, state: np.ndarray, duration_ms: float, dt_ms: float = 0.5) -> np.ndarray:
        """Voltage trace of isolated cells, sampled every step: ``(steps, cells)``."""
        out = []
        for _ in range(int(round(duration_ms / dt_ms))):
            state = self.step(state, dt_ms)
            out.append(state[0].copy())
        return np.array(out)


def burst_statistics(
    trace: np.ndarray, dt_ms: float, threshold_mv: float = -50.0
) -> tuple[float, float]:
    """Mean interval between burst onsets and mean burst length, in seconds."""
    up = trace > threshold_mv
    onsets = np.flatnonzero(np.diff(up.astype(int)) == 1)
    offsets = np.flatnonzero(np.diff(up.astype(int)) == -1)
    if len(onsets) < 2:
        return float("nan"), float("nan")
    interval = float(np.mean(np.diff(onsets))) * dt_ms / 1000.0
    lengths = [
        (offsets[offsets > a][0] - a) * dt_ms / 1000.0 for a in onsets if (offsets > a).any()
    ]
    return interval, float(np.mean(lengths)) if lengths else float("nan")
