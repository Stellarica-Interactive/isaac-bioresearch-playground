"""Single-neuron dynamics.

One model is implemented: :class:`GradedLeakyIntegrator`. It sits behind the
:class:`NeuronModel` protocol so a conductance-based model can be added later for
the handful of cells whose ion channels have actually been characterised, without
inventing channel kinetics for the other three hundred.

Why graded rather than spiking
------------------------------

Most *C. elegans* neurons do not fire action potentials. Membrane voltage varies
continuously and transmitter release varies continuously with it, because the
animal is small enough that a neuron is nearly isopotential end to end and needs
no regenerating spike to carry a signal (Goodman et al. 1998, who found no
classical action potentials in ASER or in 42 further neurons surveyed).

This is a simplification, not a fact about every neuron. AWA fires calcium-based
all-or-none action potentials (Liu et al. 2018), RMD has plateau potentials, and
AVL fires compound action potentials. Those are recorded as known deviations in
``docs/model_assumptions.md`` rather than quietly ignored.

Honesty about the parameters
----------------------------

In this lineage of models -- Wicks et al. 1996, then Kunert et al. 2014 -- the
**connectivity is the only measured quantity**. Every biophysical constant is an
order-of-magnitude assumption applied identically to all neurons. Only whole-cell
capacitance (0.5-3 pF) and input resistance (gigaohm range) have measured
support, and even those are single figures for the whole nervous system rather
than per-cell values. See ``worm/neural/parameters.toml``, where each value
carries its own provenance, and ``docs/neural_runtime.md`` for the derivation.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Protocol, runtime_checkable

import numpy as np

from common.neural.synapses import NetworkMatrices


@runtime_checkable
class NeuronModel(Protocol):
    """Dynamics of a population of cells sharing one set of equations."""

    @property
    def state_names(self) -> tuple[str, ...]:
        """Names of the per-cell state variables, in order."""
        ...

    def initial_state(self, network: NetworkMatrices) -> np.ndarray:
        """Starting state, shape ``(len(state_names), n_cells)``."""
        ...

    def derivatives(
        self, state: np.ndarray, i_ext: np.ndarray, network: NetworkMatrices
    ) -> np.ndarray:
        """Time derivative of ``state``, same shape. Units per millisecond."""
        ...


@dataclass(frozen=True, slots=True)
class GradedLeakyIntegratorParameters:
    """Biophysical constants.

    Defaults follow the Kunert et al. 2014 parameterisation of the Wicks et al.
    1996 model. **Almost none of them are measurements.** Provenance for each is in
    ``worm/neural/parameters.toml``; that file is the authority and this class is
    just its in-memory form.
    """

    c_m_pf: float = 1.0
    """Whole-cell capacitance. The one value with real measured support: Goodman
    et al. 1998 report 0.5-3 pF, so 1 pF is inside the measured range -- though
    applying a single figure to every neuron is still an assumption."""

    g_leak_ps: float = 10.0
    """Leak conductance. ASSUMED. Gives a membrane time constant of 100 ms in
    isolation, which is consistent with the gigaohm input resistances reported by
    Goodman et al. 1998 but is not a per-neuron measurement."""

    e_leak_mv: float = -35.0
    """Leak reversal potential. ASSUMED, order-of-magnitude."""

    g_syn_ps: float = 100.0
    """Conductance of a single chemical synapse. ASSUMED."""

    g_gap_ps: float = 100.0
    """Conductance of a single gap junction. ASSUMED. Rectification -- some worm
    gap junctions pass current asymmetrically -- is not modelled."""

    e_exc_mv: float = 0.0
    """Excitatory reversal potential. ASSUMED."""

    e_inh_mv: float = -45.0
    """Inhibitory reversal potential. ASSUMED."""

    beta_per_mv: float = 0.125
    """Steepness of the transmitter-release sigmoid. ASSUMED."""

    a_r_per_ms: float = 1.0e-3
    """Synaptic activation rise rate. ASSUMED. Published as 1 s^-1."""

    a_d_per_ms: float = 5.0e-3
    """Synaptic activation decay rate. ASSUMED. Published as 5 s^-1."""

    @property
    def g_leak_ns(self) -> float:
        return self.g_leak_ps / 1000.0

    @property
    def s_at_half_activation(self) -> float:
        """Steady-state synaptic activation when the release sigmoid sits at 0.5.

        Solving ``a_r phi (1 - s) = a_d s`` with ``phi = 0.5``.
        """
        half_rise = 0.5 * self.a_r_per_ms
        return half_rise / (half_rise + self.a_d_per_ms)

    def with_overrides(self, **kwargs: float) -> GradedLeakyIntegratorParameters:
        return replace(self, **kwargs)


def sigmoid(v: np.ndarray, v_half: np.ndarray, beta: float) -> np.ndarray:
    """Numerically safe logistic ``1 / (1 + exp(-beta (v - v_half)))``."""
    x = beta * (v - v_half)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    ex = np.exp(x[~pos])
    out[~pos] = ex / (1.0 + ex)
    return out


@dataclass(frozen=True)
class GradedLeakyIntegrator:
    """Graded, non-spiking single-compartment neurons.

    Two state variables per cell: membrane potential ``V`` (mV) and synaptic
    activation ``s`` (dimensionless, 0 to 1). The equations, in the units described
    in :mod:`common.neural.synapses`::

        C dV_i/dt = -G_leak (V_i - E_leak)                        leak
                    + sum_j Ggap[i,j] (V_j - V_i)                 gap junctions
                    + sum_j Gsyn[i,j] s_j (E[i,j] - V_i)          chemical synapses
                    + I_ext_i                                     injected current

        ds_i/dt   = a_r phi(V_i) (1 - s_i) - a_d s_i
        phi(V)    = 1 / (1 + exp(-beta (V - Vth_i)))

    Note ``E[i,j]``: the reversal potential is per *connection*, not per
    presynaptic neuron. The receptor that sets a synapse's sign sits on the
    postsynaptic cell, so one neuron can excite one target and inhibit another.
    Most published implementations collapse this to one value per presynaptic
    neuron; we do not, because the polarity data is per connection.

    ``Vth_i`` is not a free parameter. It is solved for, per neuron, so that the
    network sits at equilibrium with every sigmoid at its midpoint -- see
    :meth:`threshold_potentials`. That is the Kunert et al. 2014 procedure, and it
    is itself a modelling assumption: it asserts that a resting worm's neurons sit
    in the responsive middle of their range.
    """

    params: GradedLeakyIntegratorParameters = GradedLeakyIntegratorParameters()
    v_threshold_mv: np.ndarray | None = None
    """Per-neuron sigmoid midpoint. Filled in by :func:`prepare`."""

    @property
    def state_names(self) -> tuple[str, ...]:
        return ("v_mv", "s")

    # -- equilibrium -------------------------------------------------------

    def threshold_potentials(
        self, network: NetworkMatrices, i_ext: np.ndarray | None = None
    ) -> np.ndarray:
        """Solve for the resting potential of every cell, used as the sigmoid midpoint.

        With all sigmoids pinned at 0.5, synaptic activation settles at a known
        constant and the voltage equation becomes linear, so the resting state is a
        single linear solve rather than something to integrate towards::

            A V = b
            A = diag(G_leak + rowsum(Ggap) + s* rowsum(Gsyn)) - Ggap
            b = G_leak E_leak + s* rowsum(Gsyn o E) + I_ext

        Setting each neuron's threshold to its own resting potential centres every
        neuron in the responsive part of its sigmoid. **This is an assumption**, not
        a measurement -- nobody has measured resting potentials across the worm's
        nervous system.
        """
        p = self.params
        n = network.n
        s_star = p.s_at_half_activation
        i_ext = np.zeros(n) if i_ext is None else np.asarray(i_ext, dtype=np.float64)

        diag = p.g_leak_ns + network.gap_row_sum + s_star * network.g_syn.sum(axis=1)
        a = np.diag(diag) - network.g_gap
        b = p.g_leak_ns * p.e_leak_mv + s_star * network.g_syn_e_rev.sum(axis=1) + i_ext
        return np.linalg.solve(a, b)

    def initial_state(self, network: NetworkMatrices) -> np.ndarray:
        """Rest: every neuron at its solved resting potential, sigmoids at midpoint.

        There is no measured resting state of a *C. elegans* nervous system to
        initialise from. This one is self-consistent and reproducible, which is the
        most that can be claimed for it.
        """
        v_th = self._thresholds(network)
        s = np.full(network.n, self.params.s_at_half_activation)
        return np.vstack([v_th.copy(), s])

    def _thresholds(self, network: NetworkMatrices) -> np.ndarray:
        if self.v_threshold_mv is not None:
            return self.v_threshold_mv
        return self.threshold_potentials(network)

    # -- dynamics ----------------------------------------------------------

    def currents(
        self, v: np.ndarray, s: np.ndarray, i_ext: np.ndarray, network: NetworkMatrices
    ) -> dict[str, np.ndarray]:
        """The four current terms separately, in pA. Useful for inspection and teaching.

        Written as matrix-vector products against precomputed constants rather than
        the more obvious ``(g_syn * s * e_rev).sum(axis=1)``: the latter allocates an
        ``n x n`` temporary on every one of the millions of evaluations a run needs.
        The algebra is identical --

            sum_j Gsyn[i,j] s_j E[i,j]  ==  (Gsyn o E) @ s
            sum_j Gsyn[i,j] s_j         ==  Gsyn @ s
        """
        p = self.params
        i_leak = -p.g_leak_ns * (v - p.e_leak_mv)
        i_gap = network.g_gap @ v - v * network.gap_row_sum
        i_syn = network.g_syn_e_rev @ s - v * (network.g_syn @ s)
        return {"leak": i_leak, "gap": i_gap, "chemical": i_syn, "external": i_ext}

    def derivatives(
        self, state: np.ndarray, i_ext: np.ndarray, network: NetworkMatrices
    ) -> np.ndarray:
        p = self.params
        v, s = state[0], state[1]
        terms = self.currents(v, s, i_ext, network)
        dv = (terms["leak"] + terms["gap"] + terms["chemical"] + terms["external"]) / p.c_m_pf
        phi = sigmoid(v, self._thresholds(network), p.beta_per_mv)
        ds = p.a_r_per_ms * phi * (1.0 - s) - p.a_d_per_ms * s
        return np.vstack([dv, ds])

    # -- stability ---------------------------------------------------------

    def total_conductance(self, s: np.ndarray, network: NetworkMatrices) -> np.ndarray:
        """Instantaneous total membrane conductance per cell, nS."""
        return self.params.g_leak_ns + network.gap_row_sum + network.g_syn @ s

    def fastest_time_constant_ms(self, network: NetworkMatrices) -> float:
        """Shortest membrane time constant in the network, at rest.

        This sets the timestep. Explicit integration needs ``dt`` well below it; a
        hub neuron with a hundred partners can be more than a thousand times faster
        than an isolated one, so the busiest cell decides.
        """
        s = np.full(network.n, self.params.s_at_half_activation)
        g_total = self.total_conductance(s, network)
        return float(self.params.c_m_pf / np.max(g_total))


def prepare(model: GradedLeakyIntegrator, network: NetworkMatrices) -> GradedLeakyIntegrator:
    """Solve and cache the threshold potentials so they are computed once."""
    return replace(model, v_threshold_mv=model.threshold_potentials(network))
