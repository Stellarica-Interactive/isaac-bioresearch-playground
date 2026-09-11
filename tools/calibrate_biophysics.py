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

§5F supplied one. RMD's resting potential is **−69.5 mV**, measured, and a real
RMD sits there *while embedded in a real animal*. So the assumed network has a
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

**It is not a fit.** One constraint cannot determine nine parameters, and this
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

#: The six RMD cells, and the resting potential Nicoletti et al. 2019 measured.
RMD_CELLS = ("RMDDL", "RMDDR", "RMDVL", "RMDVR", "RMDL", "RMDR")
MEASURED_REST_MV = -69.5

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
    idx = np.array([runtime.network.index(c) for c in RMD_CELLS])
    for _ in range(int(settle_ms)):
        runtime.state[0, idx] = MEASURED_REST_MV
        runtime.step()
    runtime.state[0, idx] = MEASURED_REST_MV

    net = runtime.network
    v, s = runtime.state[0], runtime.state[1]
    # The network's drive decomposes as constant - conductance * V; see
    # common.neural.hybrid.HybridRuntime.network_drive.
    constant = (net.g_gap @ v + net.g_syn_e_rev @ s)[idx]
    conductance = (net.gap_row_sum + net.g_syn @ s)[idx]
    current = constant - conductance * MEASURED_REST_MV

    others = [i for i in range(net.n) if i not in set(idx.tolist())]
    return float(np.median(v[others])), float(np.abs(current).max())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--scan", action="store_true", help="search the grid")
    parser.add_argument("--top", type=int, default=8)
    args = parser.parse_args(argv)

    median, current = evaluate()
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
        "\nA region satisfying the constraint is not a fitted parameter set. One "
        "measured cell cannot determine nine parameters, and every injected current "
        "in the model was chosen against the old conductances -- see the module "
        "docstring before changing anything."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
