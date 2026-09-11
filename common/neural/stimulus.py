"""Express an input as the depolarisation it should cause, not as a current.

    drive = depolarising_current(runtime, ("AVBL", "AVBR"), 20.0)   # 20 mV worth
    runtime.inject_many(drive)

Why not just pick picoamps
--------------------------

Because a picoamp is not interpretable and a millivolt is.

A neuron's voltage lives between roughly -80 and +30 mV. Given a target in
millivolts, anybody can say whether it is plausible. Given a target in picoamps,
nobody can, because the answer depends on the cell's total conductance -- and in
this network that varies **fourteenfold**, from 1.6 nS at ALM to 23.3 nS at AVA.
The same 20 pA is a 12 mV event at ALM and a 5 mV event at PLM.

That is not a hypothetical. Every injected current in this project was chosen as a
current, and every one of them turned out wrong in a way that took a measurement
to notice:

* The touch current started at 400 pA by analogy with the command drive, and drove
  ALM to **+443 mV** (``worm/body/touch.py``).
* The proprioceptive gain of 400 pA/rad implies **124 mV per radian** at DB1 under
  the current biophysics -- a half-radian bend would swing it 62 mV.
* Part of the head-versus-tail asymmetry reported for touch in
  ``docs/model_assumptions.md`` §5D came from ALM and PLM having different
  conductances, not from the circuit.

The deeper problem is that a current chosen against one set of biophysical
parameters silently becomes wrong when those parameters change. §5G changed two of
them, which would have invalidated all three inputs at once. Expressing the
stimulus as a depolarisation and converting through the cell's own conductance
means the biophysics can move and the inputs follow.

What is still assumed
---------------------

The *size* of each target. Nothing measures how many millivolts a command
interneuron depolarises when the animal decides to move forward. What changes is
that the assumption is now stated in a unit where it can be judged, and where
being wrong by a factor of ten is obvious rather than invisible.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from common.neural.runtime import NeuralRuntime


def total_conductance(runtime: NeuralRuntime) -> np.ndarray:
    """Each cell's instantaneous total membrane conductance, nS.

    Leak plus gap junctions plus whatever chemical synapses are currently open, so
    it reflects the network's actual state rather than a resting approximation.
    """
    return np.asarray(
        runtime.model.total_conductance(  # type: ignore[attr-defined]
            runtime.state[1], runtime.network
        )
    )


def depolarising_current(
    runtime: NeuralRuntime, cells: Sequence[str], delta_mv: float
) -> dict[str, float]:
    """Current that shifts each of ``cells`` by roughly ``delta_mv``, in pA.

    ``I = dV * G`` -- nS times mV is pA, so the units already agree.

    **This is the open-loop answer**, computed with the rest of the network held
    fixed, and the network does not hold still. Measured: the current sized for a
    20 mV shift in AVB produces **+74 mV** once everything settles, because with
    ``e_exc`` at 0 mV and the network resting near -31 mV, excitation feeds back
    positively and amplifies. The obvious intuition -- that a cell embedded in a
    network is harder to move than an isolated one -- is wrong here, and was
    written into this docstring before being checked.

    Use :func:`solve_for_depolarisation` when the achieved shift is what matters.
    """
    conductance = total_conductance(runtime)
    return {
        cell: float(delta_mv * conductance[runtime.network.index(cell)])
        for cell in cells
        if cell in runtime.network.cell_ids
    }


def scale_for_depolarisation(
    runtime: NeuralRuntime, cells: Sequence[str], delta_mv: float
) -> float:
    """One current, in pA, sized for the *median* cell in ``cells``.

    For stimuli that must be a single number rather than a per-cell mapping -- a
    proprioceptive gain, say, shared across every B-type neuron. Uses the median so
    one unusually well-connected cell does not set the scale for the rest.
    """
    conductance = total_conductance(runtime)
    indices = [runtime.network.index(c) for c in cells if c in runtime.network.cell_ids]
    if not indices:
        raise ValueError(f"none of {tuple(cells)} are in this network")
    return float(delta_mv * np.median(conductance[indices]))


def solve_for_depolarisation(
    runtime: NeuralRuntime,
    cells: Sequence[str],
    delta_mv: float,
    *,
    settle_ms: float = 2000.0,
    tolerance_mv: float = 0.5,
    max_iterations: int = 12,
) -> dict[str, float]:
    """Current that actually shifts ``cells`` by ``delta_mv`` once the network settles.

    :func:`depolarising_current` answers the open-loop question and can be out by a
    factor of three or more, in either direction: excitatory feedback amplifies,
    heavy shunting attenuates, and which dominates depends on the cell.

    So this measures instead of predicting. It injects, lets the network settle,
    compares the achieved shift against the target and rescales, which converges in
    a handful of iterations because the relationship is close to monotonic.

    Expensive -- each iteration is a full settle -- so call it once when a run is
    configured, not inside a loop.
    """
    present = [c for c in cells if c in runtime.network.cell_ids]
    if not present:
        raise ValueError(f"none of {tuple(cells)} are in this network")

    baseline_state = runtime.state.copy()
    baseline_input = runtime.i_ext_pa.copy()
    baseline = {c: runtime.voltage(c) for c in present}

    current = depolarising_current(runtime, present, delta_mv)
    try:
        for _ in range(max_iterations):
            runtime.state = baseline_state.copy()
            runtime.i_ext_pa = baseline_input.copy()
            runtime.inject_many(current)
            runtime.run(settle_ms)

            achieved = float(np.median([runtime.voltage(c) - baseline[c] for c in present]))
            if abs(achieved - delta_mv) <= tolerance_mv:
                break
            if abs(achieved) < 1e-9:
                current = {c: v * 2.0 for c, v in current.items()}
                continue
            # Damped, because the response is not exactly linear and an undamped
            # correction can oscillate around the target instead of reaching it.
            factor = 1.0 + 0.6 * (delta_mv / achieved - 1.0)
            current = {c: v * float(np.clip(factor, 0.2, 5.0)) for c, v in current.items()}
    finally:
        runtime.state = baseline_state
        runtime.i_ext_pa = baseline_input

    return current
