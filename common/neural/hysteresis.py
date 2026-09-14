"""Bistable output for named cells, as a hypothesis test rather than a model.

    switch = Hysteresis.build(runtime, cells=proprio.targets, width_mv=4.0)
    runtime.step()
    switch.apply(runtime)      # after every step

**This is not biology and must never become a default.** It is Boyle, Berri &
Cohen's modelling choice, imported to find out whether it is the missing
ingredient in this model -- and the answer is worth having either way.

Why it is here
--------------

Three independent lines of evidence in ``docs/model_assumptions.md`` point at the
same gap:

* §5C.9 -- the graded leaky integrator has no limit cycle anywhere, under any
  drive. It is a pure relaxation system.
* §5M -- Boyle et al. build a neuromechanical model with no central pattern
  generator that undulates, and its oscillator is **binary B-class motor neurons
  with hysteresis** (`S = 1 if I > 0.5 + ε(0.5 − S)`, switching on at 0.75 and off
  at 0.25), reset by antagonistic D-class inhibition.
* §5P -- this model does not merely fail to start a wave, it destroys a seeded
  one within five seconds. Its feedback restores toward uniformity, and nothing in
  it can hold a phase difference between neighbouring segments.

Hysteresis is exactly a mechanism for holding such a difference: a latched cell
stays latched while its neighbour sits in the other state. Whether real B-type
neurons have it is unknown to this project -- VB6 in Nicoletti et al. 2024 is the
one published conductance-based model of a B-class cell, and its channel mapping is
the ambiguity of §5I.

What the two outcomes mean
--------------------------

**If a gait appears**, the hypothesis survives and the VB6 datum becomes the thing
that decides whether the rhythm can be attributed to biology or must be attributed
to this modelling choice. The gait itself would be Boyle's, not the connectome's,
and could not honestly be reported as the latter.

**If no gait appears**, three converging lines of evidence are wrong together, and
the missing ingredient is something nobody in this document has yet identified.
That is the more valuable outcome and the reason to run it.

How it works here
-----------------

The graded model's output variable is the synaptic activation ``s``, a saturating
function of membrane potential. This replaces that function, for named cells only,
with a Schmitt trigger on voltage: rise above ``centre + width/2`` and the cell
latches on, fall below ``centre − width/2`` and it latches off, and in between it
holds whatever state it already had.

Two placement choices decide whether the test measures anything, and both were
wrong in the first version:

* **Where the band sits.** Centred on each cell's activation threshold -- the
  obvious choice -- it sat 30 mV above where the cells ever go (§5N), so nothing
  switched and every cell stayed latched low. That silences the population; it
  does not make it bistable. The band is now centred on the driven operating
  point. ``fraction_on`` pinned at 0 **or** 1 still means the band is misplaced.
* **What the two states emit.** Boyle's binary neurons output 0 and 1. These
  cells' graded activation never leaves ``[0, 0.165]``, because the sigmoid is
  centred 30 mV above rest and even their most depolarised excursion is in its
  bottom tail. Emitting 0/1 would multiply the muscle drive roughly sixfold at
  the same moment as it introduced memory, and no later measurement could
  separate the two. Each state therefore emits the graded activation the cell
  would settle at if held at the band edge it last crossed, leaving bistability
  as the only difference. ``binary_output=True`` asks the gain question instead.

Applied *after* each integration step rather than inside the derivative, so the
runtime stays exactly as it was and the override is visible at the call site. That
also means it is a discrete map composed with a continuous flow, which is not a
faithful way to model a bistable membrane -- the conductance-based route of
:mod:`common.neural.conductance` is. It is adequate for asking whether bistability
in these cells changes the outcome.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from common.neural.neuron_models import sigmoid
from common.neural.runtime import NeuralRuntime

#: Width of the hysteresis band, mV. Boyle et al. use ``ε_hys = 0.5`` on a unit
#: input scale, switching on at 0.75 and off at 0.25 -- half the full range. There
#: is no principled translation of that into millivolts here, so this is ASSUMED
#: and is the parameter to sweep before concluding anything from a null result.
DEFAULT_WIDTH_MV = 4.0


@dataclass
class Hysteresis:
    """Latches named cells' output high or low, with a dead band between."""

    indices: np.ndarray
    on_mv: np.ndarray
    off_mv: np.ndarray
    #: Activation to emit in each state. Not 0 and 1: see :meth:`build`.
    on_output: np.ndarray = field(default_factory=lambda: np.array([1.0]))
    off_output: np.ndarray = field(default_factory=lambda: np.array([0.0]))
    cells: tuple[str, ...] = ()
    latched: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.latched = np.zeros(self.indices.shape, dtype=bool)

    @classmethod
    def build(
        cls,
        runtime: NeuralRuntime,
        cells: tuple[str, ...],
        *,
        width_mv: float = DEFAULT_WIDTH_MV,
        centre_on_threshold: bool = False,
        binary_output: bool = False,
    ) -> Hysteresis:
        """Centre the switch where the cells **currently sit**, not on their threshold.

        This was originally centred on each cell's activation threshold, which is
        the obvious choice and is wrong. That threshold sits
        ``v_threshold_offset_mv`` above rest -- 30 mV, since §5N -- and the B-type
        neurons operate far below it, so nothing ever crossed and every cell stayed
        latched low. The effect was to *silence* the population rather than make it
        bistable, and two different dead-band widths produced byte-identical runs
        because in both the output was uniformly zero.

        Centring on the present membrane potential puts the band where the cells
        actually live. Build it after settling the network, so "present" means the
        operating point rather than an initial condition.

        ``centre_on_threshold`` restores the original behaviour, which is useful
        only for demonstrating the failure.
        """
        present = tuple(c for c in cells if c in runtime.network.cell_ids)
        if not present:
            raise ValueError(f"none of {cells} are in this network")
        idx = np.array([runtime.network.index(c) for c in present], dtype=np.int64)
        centre = (
            np.asarray(
                runtime.model.v_threshold_mv,  # type: ignore[attr-defined]
                dtype=np.float64,
            )[idx]
            if centre_on_threshold
            else np.asarray(runtime.state[0, idx], dtype=np.float64).copy()
        )
        on_mv = centre + 0.5 * width_mv
        off_mv = centre - 0.5 * width_mv

        # Emit the *graded* activation the cell would have had at the band edge it
        # last crossed, rather than Boyle's literal 0 and 1.
        #
        # This matters more than it looks. Measured over a baseline run, these
        # cells' graded activation never leaves [0, 0.165] -- the sigmoid is
        # centred 30 mV above rest (5N), so even the most depolarised excursion
        # sits in its bottom tail. Substituting 0/1 would therefore multiply the
        # muscle drive by about six at the same time as it introduced memory, and
        # any change in gait could be attributed to either. Matching the output to
        # the two switch points leaves bistability as the only thing that changed,
        # which is the question being asked. ``binary_output`` restores Boyle's
        # scale for comparison, and answers the gain question instead.
        if binary_output:
            on_output = np.ones_like(on_mv)
            off_output = np.zeros_like(off_mv)
        else:
            on_output = _steady_activation(runtime, idx, on_mv)
            off_output = _steady_activation(runtime, idx, off_mv)

        return cls(
            indices=idx,
            on_mv=on_mv,
            off_mv=off_mv,
            on_output=on_output,
            off_output=off_output,
            cells=present,
        )

    def apply(self, runtime: NeuralRuntime) -> None:
        """Overwrite the latched cells' activation. Call after every step."""
        v = runtime.state[0, self.indices]
        self.latched = np.where(
            v > self.on_mv, True, np.where(v < self.off_mv, False, self.latched)
        )
        runtime.state[1, self.indices] = np.where(self.latched, self.on_output, self.off_output)

    @property
    def fraction_on(self) -> float:
        """How much of the population is currently latched high.

        Worth watching: a value pinned at 0 or 1 means the dead band is in the
        wrong place and nothing is switching, which would produce a null result
        for a reason that has nothing to do with the hypothesis.
        """
        return float(self.latched.mean()) if self.latched.size else 0.0


def _steady_activation(runtime: NeuralRuntime, idx: np.ndarray, v_mv: np.ndarray) -> np.ndarray:
    """Synaptic activation a cell settles at while held at ``v_mv``.

    Solving ``a_r phi(V) (1 - s) = a_d s`` for ``s``, with the same sigmoid and
    the same rate constants the runtime integrates. Sharing the formula rather
    than restating it is deliberate: a latch whose output scale drifted away from
    the model it replaces would be a gain experiment wearing a bistability label.
    """
    params = runtime.model.params
    thresholds = np.asarray(
        runtime.model.v_threshold_mv,  # type: ignore[attr-defined]
        dtype=np.float64,
    )[idx]
    phi = sigmoid(v_mv, thresholds, params.beta_per_mv)
    a_r, a_d = params.a_r_per_ms, params.a_d_per_ms
    return np.asarray(a_r * phi / (a_d + a_r * phi), dtype=np.float64)
