"""Test the model's assumed biophysics against a cell whose real values are known.

    python tools/calibrate_biophysics.py                 # where we stand
    python tools/calibrate_biophysics.py --scan          # search the parameter grid

Why this exists
---------------

Nine of the ten parameters in ``worm/neural/parameters.toml`` are ASSUMED --
order-of-magnitude values from Kunert et al. 2014, in a modelling lineage where
connectivity is the only measured quantity. ``docs/model_assumptions.md`` §5A.1
showed the answer moves fourteenfold across defensible choices of them, and until
recently there was nothing to do about that: no measured quantity in the model
could tell one choice from another.

§5F supplied one, and AWC^on a second. Their resting potentials are **−69.5 mV**
and **−69.2 mV**, from separate channel complements fitted to separate data, and a
real RMD or AWC sits there *while embedded in a real animal*. So the assumed network has a
number it must reproduce, and it does not: it drags RMD to −22 mV, delivering
hundreds of picoamps into a cell whose measured bistability survives about four.

That is a constraint, and this tool applies it.

What it measures
----------------

Hold RMD where the measurement says it sits, let the rest of the network settle
around it, and ask how much current the network pushes into it. Below roughly 4 pA
the measured cell can still do what it was measured doing; far above that, its
biophysics is overwhelmed by our assumptions about everything else.

The evaluation is deliberately cheap -- no conductance substepping -- so a wide
grid is affordable.

What this is not
----------------

**It measures a clamped state, not an equilibrium.** The number below is the
current needed to *hold* a measured cell where its measurement says it sits. It is
not evidence that the network would rest there, and it must not be read as one:
free of any clamp the network settles at -31.3 mV with the committed parameters,
some 38 mV away. See docs/model_assumptions.md 5K, which corrects an earlier
reading of exactly this output. A correct optimisation target is the *free*
resting potential.

**It is not a fit.** Two constraints cannot determine nine parameters, and this
tool does not pretend otherwise: it reports where in the grid the constraint is
satisfied, and leaves the choice, and the argument for it, to a human.

It also cannot be read in isolation from the injected currents. The command drive
(500 pA), the touch current (20 pA) and the proprioceptive gain (400 pA/rad) were
all chosen against the *present* conductances. Lower ``g_syn`` tenfold and the
same 500 pA drives AVB to +100 mV, which is not a neuron. Any change here
invalidates those numbers and they must be re-derived, not carried over.
"""

from __future__ import annotations

import argparse
import itertools
import sys

import numpy as np

from common.data.schemas import CellCategory
from common.neural.synapses import UnknownSignPolicy
from worm.importers.naming import body_wall_muscle_ids
from worm.loader import load
from worm.neural.config import RUNTIME_OVERLAYS, build_runtime

#: Cells whose resting potential we have from a conductance-based model, and the
#: value that model settles to with no input.
#:
#: Two independent measurements, which matters more than two numbers. RMD and
#: AWC^on carry different channel complements -- EGL-36 against EGL-2, KVS-1 and
#: KQT-3 -- and were fitted to separate voltage-clamp data. That they both settle
#: within a third of a millivolt of -69 mV is not something either fit was told to
#: do, and it makes "the network should rest near -69 mV" a far stronger claim than
#: one cell could support.
MEASURED_CELLS: dict[str, tuple[tuple[str, ...], float]] = {
    "rmd_nicoletti2019": (("RMDDL", "RMDDR", "RMDVL", "RMDVR", "RMDL", "RMDR"), -69.5),
    "awc_nicoletti2019": (("AWCL", "AWCR"), -69.2),
}

#: Kept for callers that predate the second cell.
RMD_CELLS = MEASURED_CELLS["rmd_nicoletti2019"][0]
MEASURED_REST_MV = MEASURED_CELLS["rmd_nicoletti2019"][1]

#: Standing bias RMD's plateau bistability survives, pA. Measured in §5F.2: two
#: stable states at 0 and 2 pA, one only from 4 pA upward.
BISTABILITY_BUDGET_PA = 4.0

#: The grid. Deliberately coarse and centred on the committed values, because the
#: point is to find out whether a consistent region exists at all, not to optimise.
GRID: dict[str, tuple[float, ...]] = {
    "g_leak_ps": (10.0, 100.0, 500.0),
    "g_syn_ps": (100.0, 30.0, 10.0),
    "g_gap_ps": (100.0, 30.0, 10.0),
    "e_leak_mv": (-35.0, -60.0, -70.0),
}


def evaluate(settle_ms: float = 3000.0, **overrides: float) -> tuple[float, float]:
    """``(network resting median mV, largest current into an RMD cell pA)``.

    RMD is clamped at its measured resting potential throughout, which is the
    situation the real animal is in: the cell holds that voltage, and everything
    else arranges itself around it.
    """
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
        parameter_overrides=overrides,
    )
    clamped: list[tuple[np.ndarray, float]] = []
    for cells_, rest in MEASURED_CELLS.values():
        present = [c for c in cells_ if c in runtime.network.cell_ids]
        if present:
            clamped.append((np.array([runtime.network.index(c) for c in present]), rest))

    for _ in range(int(settle_ms)):
        for idx, rest in clamped:
            runtime.state[0, idx] = rest
        runtime.step()
    for idx, rest in clamped:
        runtime.state[0, idx] = rest

    net = runtime.network
    v, s = runtime.state[0], runtime.state[1]
    # The network's drive decomposes as constant - conductance * V; see
    # common.neural.hybrid.HybridRuntime.network_drive.
    constant = net.g_gap @ v + net.g_syn_e_rev @ s
    conductance = net.gap_row_sum + net.g_syn @ s

    worst = 0.0
    held: set[int] = set()
    for idx, rest in clamped:
        current = constant[idx] - conductance[idx] * rest
        worst = max(worst, float(np.abs(current).max()))
        held |= set(idx.tolist())

    others = [i for i in range(net.n) if i not in held]
    return float(np.median(v[others])), worst


def evaluate_free(settle_ms: float = 4000.0, **overrides: float) -> tuple[float, float]:
    """``(resting median mV, spread mV)`` with **nothing clamped**.

    The objective :func:`evaluate` should have been. Clamping a cell and measuring
    the current needed to hold it there asks how far the network is from the
    measurement, which is a diagnostic; a parameter set can minimise it while the
    free network barely moves. This settles everything and reads off where it
    actually rests, which is the quantity that has to equal about -69 mV.

    Slower, because it cannot stop early.
    """
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
        parameter_overrides=overrides,
    )
    runtime.run(settle_ms)
    v = runtime.state[0]
    if not np.all(np.isfinite(v)):
        return float("nan"), float("nan")
    return float(np.median(v)), float(np.percentile(v, 90) - np.percentile(v, 10))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--scan", action="store_true", help="search the grid")
    parser.add_argument("--top", type=int, default=8)
    args = parser.parse_args(argv)

    median, current = evaluate()
    names = ", ".join(
        f"{'/'.join(c[:2])}{'...' if len(c) > 2 else ''} at {r:g} mV"
        for c, r in MEASURED_CELLS.values()
    )
    print(f"Constraints: {names}")
    print("With the committed parameters:")
    print(f"  network resting median       {median:8.1f} mV")
    print(f"  current into RMD held at rest {current:7.1f} pA")
    print(f"  budget before RMD loses its bistability {BISTABILITY_BUDGET_PA:.0f} pA")
    print(
        f"  -> over budget by {current / BISTABILITY_BUDGET_PA:.0f}x"
        if current > BISTABILITY_BUDGET_PA
        else "  -> within budget"
    )
    if not args.scan:
        print("\nRun with --scan to search for parameter values that satisfy it.")
        return 0

    print(f"\nScanning {np.prod([len(v) for v in GRID.values()]):.0f} combinations ...\n")
    results = []
    for values in itertools.product(*GRID.values()):
        overrides = dict(zip(GRID.keys(), values, strict=True))
        med, cur = evaluate(**overrides)
        results.append((cur, med, overrides))

    results.sort(key=lambda r: r[0])
    print(f"{'into RMD':>9s} {'net median':>11s}  parameters")
    for cur, med, overrides in results[: args.top]:
        flag = "  <- within budget" if cur <= BISTABILITY_BUDGET_PA else ""
        summary = ", ".join(
            f"{k.replace('_ps', '').replace('_mv', '')}={v:g}" for k, v in overrides.items()
        )
        print(f"{cur:8.2f}pA {med:10.1f}mV  {summary}{flag}")

    print(
        "\nA region satisfying the constraint is not a fitted parameter set. Two "
        "measured cells cannot determine nine parameters, and every injected current "
        "in the model was chosen against the old conductances -- see the module "
        "docstring before changing anything."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
