"""Turning a connectome into the matrices a dynamical model needs.

A :class:`~common.data.schemas.Connectome` is a list of measured contacts. A
simulation needs conductance matrices. Every step of that conversion is a
modelling assumption, and this module is where they are made explicit.

Units
-----

The whole runtime uses one consistent set, chosen so the numbers stay readable:

===========  ======  ==================================================
Quantity     Unit    Note
===========  ======  ==================================================
voltage      mV
time         ms
capacitance  pF      whole-cell, 0.5-3 pF measured (Goodman et al. 1998)
conductance  nS      published values are in pS; converted on load
current      pA
===========  ======  ==================================================

These are mutually consistent: ``C dV/dt`` in ``pF·mV/ms`` is ``pA``, and ``G·V``
in ``nS·mV`` is also ``pA``.

The three assumptions made here
-------------------------------

1. **Synapse count maps linearly to conductance.** ``g = g_unit x n_synapses``.
   Nobody has measured this relationship. A saturating map is equally plausible;
   :class:`WeightScaling` exists so the choice can be varied and its effect
   measured rather than assumed away.

2. **Sign lives on the connection, not on the presynaptic neuron.** Most published
   whole-network models give each *neuron* one reversal potential. That is
   backwards: the receptor that decides the sign is on the *postsynaptic* cell, so
   the same neuron can excite one target and inhibit another. We keep a full
   reversal-potential matrix, which is what the Fenyves data actually provides.

3. **Unknown signs need an explicit policy.** Roughly half of chemical connections
   have no predicted sign, and every neuromuscular junction is unsigned. There is
   no defensible default, so :class:`UnknownSignPolicy` is a required argument.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np

from common.data.schemas import Connectome, Sign, SynapseType

PS_PER_NS = 1000.0


class UnknownSignPolicy(StrEnum):
    """What to do with a chemical synapse whose sign nobody has predicted.

    There is no right answer, which is exactly why this must be chosen explicitly
    and reported alongside any result.

    ``EXCLUDE``
        Drop the connection. The network loses real anatomy, but nothing is
        invented. The most conservative option and the recommended default for a
        first experiment.

    ``EXCITATORY`` / ``INHIBITORY``
        Assume one sign for all of them. Fast, and badly wrong somewhere. Among
        connections Fenyves *can* call, excitatory outnumber inhibitory about
        3.1:1 by connection count, so ``EXCITATORY`` is the less unreasonable
        guess -- but it is still a guess about ~1900 connections.

    ``NEUTRAL``
        Keep the connection with its reversal potential set to the leak potential,
        so it shunts (increases conductance, pulls toward rest) without pushing
        either way. Arguably the least committal way to keep the anatomy.
    """

    EXCLUDE = "exclude"
    EXCITATORY = "excitatory"
    INHIBITORY = "inhibitory"
    NEUTRAL = "neutral"


class WeightScaling(StrEnum):
    """How a synapse count becomes a conductance.

    ``LINEAR``
        ``g = g_unit x n``. The usual choice, and an assumption.

    ``SQRT``
        ``g = g_unit x sqrt(n)``. A cheap stand-in for saturation; included so the
        sensitivity of a result to this choice can be measured.

    ``BINARY``
        ``g = g_unit`` for any n >= 1. Discards weight entirely, which is a useful
        control: if a behaviour survives this, it did not depend on synapse counts.
    """

    LINEAR = "linear"
    SQRT = "sqrt"
    BINARY = "binary"


def _scale(counts: np.ndarray, scaling: WeightScaling) -> np.ndarray:
    if scaling is WeightScaling.LINEAR:
        return counts.astype(np.float64)
    if scaling is WeightScaling.SQRT:
        return np.sqrt(counts.astype(np.float64))
    return (counts > 0).astype(np.float64)


@dataclass(frozen=True)
class SignCensus:
    """How many connections got their sign from data, and how many from a policy."""

    excitatory: int
    inhibitory: int
    mixed: int
    from_policy: int
    excluded: int

    @property
    def total(self) -> int:
        return self.excitatory + self.inhibitory + self.mixed + self.from_policy + self.excluded

    @property
    def measured_fraction(self) -> float:
        """Fraction whose sign came from the polarity data rather than our policy.

        Named ``measured`` only in contrast to ``from_policy``. Even these are
        *predictions* from gene expression; nothing here is a measurement of sign.
        """
        known = self.excitatory + self.inhibitory + self.mixed
        return known / self.total if self.total else 0.0

    def summary(self) -> str:
        return (
            f"{self.excitatory} excitatory, {self.inhibitory} inhibitory, "
            f"{self.mixed} mixed, {self.from_policy} assigned by policy, "
            f"{self.excluded} excluded "
            f"({self.measured_fraction:.1%} from the polarity data)"
        )


@dataclass(frozen=True)
class NetworkMatrices:
    """Everything a :class:`~common.neural.runtime.NeuralRuntime` needs about wiring.

    All matrices are indexed ``[post, pre]`` -- row ``i`` collects the inputs *to*
    cell ``i``. That is the direction current flows, and it makes every update a
    plain matrix-vector product.
    """

    cell_ids: tuple[str, ...]
    g_gap: np.ndarray
    """Symmetric gap junction conductance, nS. ``[i, j]`` couples i and j."""

    g_syn: np.ndarray
    """Chemical synapse conductance, nS. ``[i, j]`` is the synapse from j onto i."""

    e_rev: np.ndarray
    """Reversal potential per connection, mV. Only meaningful where ``g_syn > 0``."""

    census: SignCensus
    unknown_sign_policy: UnknownSignPolicy
    weight_scaling: WeightScaling

    # Precomputed constants. The derivative is evaluated four times per RK4 step and
    # millions of times per experiment, so anything that does not change with state
    # is computed once here. Without this, each evaluation allocates several n x n
    # temporaries and the whole runtime is roughly thirty times slower.
    gap_row_sum: np.ndarray = field(init=False, repr=False)
    """``sum_j Ggap[i,j]``, nS. Constant."""

    g_syn_e_rev: np.ndarray = field(init=False, repr=False)
    """``Gsyn o E``, so the synaptic driving term is a plain matrix-vector product."""

    def __post_init__(self) -> None:
        object.__setattr__(self, "gap_row_sum", self.g_gap.sum(axis=1))
        object.__setattr__(self, "g_syn_e_rev", self.g_syn * self.e_rev)

    def index(self, cell_id: str) -> int:
        try:
            return self.cell_ids.index(cell_id)
        except ValueError:
            raise KeyError(f"cell not in this network: {cell_id!r}") from None

    @property
    def n(self) -> int:
        return len(self.cell_ids)

    def describe(self) -> str:
        gap_pairs = int(np.count_nonzero(np.triu(self.g_gap)))
        return (
            f"{self.n} cells, {int(np.count_nonzero(self.g_syn))} chemical synapses, "
            f"{gap_pairs} gap junctions\n"
            f"  sign: {self.census.summary()}\n"
            f"  unknown-sign policy: {self.unknown_sign_policy}, "
            f"weight scaling: {self.weight_scaling}"
        )


def build_matrices(
    connectome: Connectome,
    *,
    unknown_sign: UnknownSignPolicy,
    g_syn_ps: float,
    g_gap_ps: float,
    e_exc_mv: float,
    e_inh_mv: float,
    e_leak_mv: float,
    weight_scaling: WeightScaling = WeightScaling.LINEAR,
    cells: tuple[str, ...] | None = None,
) -> NetworkMatrices:
    """Build conductance and reversal-potential matrices from a connectome.

    ``unknown_sign`` has no default on purpose. Around half of chemical connections
    carry no predicted sign, and every neuromuscular junction is unsigned, so this
    choice materially changes the network and must be a conscious one.

    A ``MIXED`` sign -- the postsynaptic cell expresses both excitatory and
    inhibitory receptors for that transmitter -- is treated as a shunt at the leak
    potential. That is a modelling choice, not a finding: the source is telling us
    the net effect is undetermined, and a shunt is the least committal thing to do
    with it while keeping the anatomy.
    """
    ids = tuple(cells) if cells is not None else tuple(connectome.cell_ids())
    index = {c: i for i, c in enumerate(ids)}
    n = len(ids)

    gap_counts = np.zeros((n, n), dtype=np.int64)
    syn_counts = np.zeros((n, n), dtype=np.int64)
    e_rev = np.zeros((n, n), dtype=np.float64)

    excitatory = inhibitory = mixed = from_policy = excluded = 0

    for e in connectome.connections:
        i, j = index.get(e.post), index.get(e.pre)
        if i is None or j is None:
            continue

        if e.synapse_type is SynapseType.ELECTRICAL:
            if i == j:
                # A gap junction from a cell to itself is electrically inert in a
                # single-compartment model: the driving force (V_j - V_i) is zero.
                # Kept in the data, dropped here, and said so rather than silently.
                continue
            gap_counts[i, j] += e.weight
            gap_counts[j, i] += e.weight
            continue

        if e.sign is Sign.EXCITATORY:
            reversal = e_exc_mv
            excitatory += 1
        elif e.sign is Sign.INHIBITORY:
            reversal = e_inh_mv
            inhibitory += 1
        elif e.sign is Sign.MIXED:
            reversal = e_leak_mv
            mixed += 1
        else:
            if unknown_sign is UnknownSignPolicy.EXCLUDE:
                excluded += 1
                continue
            if unknown_sign is UnknownSignPolicy.EXCITATORY:
                reversal = e_exc_mv
            elif unknown_sign is UnknownSignPolicy.INHIBITORY:
                reversal = e_inh_mv
            else:
                reversal = e_leak_mv
            from_policy += 1

        syn_counts[i, j] += e.weight
        e_rev[i, j] = reversal

    g_gap = _scale(gap_counts, weight_scaling) * (g_gap_ps / PS_PER_NS)
    g_syn = _scale(syn_counts, weight_scaling) * (g_syn_ps / PS_PER_NS)

    return NetworkMatrices(
        cell_ids=ids,
        g_gap=g_gap,
        g_syn=g_syn,
        e_rev=e_rev,
        census=SignCensus(excitatory, inhibitory, mixed, from_policy, excluded),
        unknown_sign_policy=unknown_sign,
        weight_scaling=weight_scaling,
    )
