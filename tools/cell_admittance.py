"""Can a measured cell's own channels inject energy at the gait frequency?

    .venv/Scripts/python tools/cell_admittance.py
    .venv/Scripts/python tools/cell_admittance.py --models vb6_nicoletti2024 --frequencies 0.1 0.5

model_assumptions 5AN.9 found that the forward gait mode's damping is not set by
any connection weight: the loop needs something that injects energy *at the gait
frequency*. 5AC asked whether a conductance-based VB6 oscillates on its own and
found it does not -- but a cell does not have to oscillate on its own to
destabilise a loop. It only has to stop being a resistor at the right frequency.

That is a property of the cell alone, and it is exactly computable. Hold the cell
at a voltage ``V0`` with whatever steady current that takes, linearise every
channel, gate and calcium pool around that point, and the cell's response to a
small sinusoidal current is its admittance ``Y(i 2 pi f)`` in nS. The real part
is the conductance the cell presents at that frequency:

* ``Re Y > 0`` everywhere: the cell dissipates at every frequency. Coupled into
  any loop it can only add damping. The graded leaky integrator is this, with
  ``Re Y = G_leak`` at every frequency.
* ``Re Y < 0`` at ``f = 0``: negative slope conductance -- regenerative, the
  ingredient of a latch or a plateau (RMD's bistability, 5F).
* ``Re Y < 0`` in a band away from zero: the cell *amplifies* in that band. This
  is the one that can destabilise an oscillatory mode, if the band covers the
  mode's frequency and the negative conductance outweighs what the network
  loads the cell with.

Every number comes from the published models as imported; nothing is fitted.
Holding at ``V0`` is a probe, not a claim about where the cell sits in the
animal: the loop holds B-type cells depolarised, by an amount the graded model
sets, so the scan covers the range.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, ".")
from common.neural.conductance import ConductanceModel  # noqa: E402

MODELS = (
    "vb6_nicoletti2024",
    "va5_nicoletti2024",
    "vd5_nicoletti2024",
    "aval_nicoletti2024",
    "avar_nicoletti2024",
    "rim_nicoletti2024",
    "aiy_nicoletti2024",
    "rmd_nicoletti2019",
    "awc_nicoletti2019",
)

#: The 2024 calcium pool's floor (cadiff.mod), mM. A pool sitting on it is
#: clamped, not dynamic, and is held fixed in the linearisation.
CA_FLOOR_MM = 1e-4


def clamped_steady_state(model: ConductanceModel, v0: float) -> tuple[np.ndarray, bool]:
    """Gates and calcium at steady state with the membrane held at ``v0``.

    A clamped run gets close; damped Newton on everything but ``V`` finishes.
    """
    state = model.resting_state(1).copy()
    for _ in range(40000):  # 2 s at the default 0.05 ms
        state[0] = v0
        state = model.step(state, 0.05)
    state[0] = v0
    rows = np.arange(1, model.n_state)
    for _ in range(30):
        r = model.derivatives(state)[1:, 0]
        if np.abs(r).max() < 1e-12:
            break
        jac = np.zeros((rows.size, rows.size))
        for k, row in enumerate(rows):
            h = 1e-7 * max(1e-3, abs(state[row, 0]))
            plus, minus = state.copy(), state.copy()
            plus[row] += h
            minus[row] -= h
            jac[:, k] = (model.derivatives(plus)[1:, 0] - model.derivatives(minus)[1:, 0]) / (2 * h)
        step = np.linalg.lstsq(jac, -r, rcond=None)[0]
        state[1:, 0] += step
        state[1:, 0] = np.clip(state[1:, 0], 0.0, None)
    ca_row = model._row.get("ca_intra1")
    at_floor = (
        ca_row is not None
        and model._uses_2024_calcium_pool
        and (state[ca_row, 0] <= CA_FLOOR_MM * (1 + 1e-9))
    )
    return state, at_floor


def linearise(model: ConductanceModel, state: np.ndarray, i_hold: float, frozen: set[int]):
    """Jacobian of the full field at ``state`` under the holding current, per ms."""
    n = model.n_state
    a = np.zeros((n, n))
    for k in range(n):
        if k in frozen:
            continue
        h = 1e-6 * max(1e-3, abs(state[k, 0]))
        plus, minus = state.copy(), state.copy()
        plus[k] += h
        minus[k] -= h
        a[:, k] = (
            model.derivatives(plus, i_hold)[:, 0] - model.derivatives(minus, i_hold)[:, 0]
        ) / (2 * h)
    # A state with no dynamics at all -- an all-zero row -- is a constant here,
    # and leaving it in makes ``s I - A`` singular at 0 Hz. AVAL has one.
    dead = {k for k in range(1, n) if not np.any(a[k])}
    keep = [k for k in range(n) if k not in frozen and k not in dead]
    return a[np.ix_(keep, keep)]


def admittance(a: np.ndarray, c_pf: float, freqs_hz: np.ndarray) -> np.ndarray:
    """``Y = I / V`` in nS for a current injected into the membrane, row 0."""
    out = np.zeros(freqs_hz.size, dtype=complex)
    e0 = np.zeros(a.shape[0])
    e0[0] = 1.0
    for k, f in enumerate(freqs_hz):
        s = 2j * np.pi * f / 1000.0  # per ms
        h = np.linalg.solve(s * np.eye(a.shape[0]) - a, e0 / c_pf)[0]  # mV per pA
        out[k] = 1.0 / h
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--models", nargs="+", default=list(MODELS))
    parser.add_argument("--holding", type=float, nargs="+", default=list(range(-80, -4, 5)))
    parser.add_argument(
        "--frequencies", type=float, nargs="+", default=[0.0, 0.05, 0.17, 0.5, 2.0, 10.0]
    )
    args = parser.parse_args()
    freqs = np.array(args.frequencies)

    for model_id in args.models:
        model = ConductanceModel.load(model_id)
        c = model.p["c"]
        print(f"\n=== {model_id}  (C = {c:g} pF, {len(model.channels)} channels)")
        head = (
            f"  {'V0':>6} {'I hold':>9} "
            + " ".join(f"{f'ReY {f:g}Hz':>11}" for f in freqs)
            + f" {'own growth':>11} {'own freq':>9}"
        )
        print(head)
        print("  " + "-" * (len(head) - 2))
        for v0 in args.holding:
            state, at_floor = clamped_steady_state(model, float(v0))
            i_hold = float(sum(model.currents(state).values())[0])
            frozen = {model._row["ca_intra1"]} if at_floor else set()
            a = linearise(model, state, i_hold, frozen)
            y = admittance(a, c, freqs)
            # The cell's own stability when held there by current alone.
            lam = np.linalg.eigvals(a) * 1000.0  # per s
            lead = lam[np.argmax(lam.real)]
            cells = " ".join(f"{v.real:+11.4f}" for v in y)
            flag = "  [Ca at floor]" if at_floor else ""
            print(
                f"  {v0:6.0f} {i_hold:+9.3f} {cells} {lead.real:+11.3f} "
                f"{abs(lead.imag) / (2 * np.pi):8.3f}Hz{flag}"
            )
    print()
    print("  ReY in nS: the conductance the cell presents at that frequency.")
    print("  Negative away from 0 Hz = the cell amplifies there; negative at 0 Hz =")
    print("  regenerative. 'own growth' > 0 = the cell is unstable held there by")
    print("  current alone, so the row describes an operating point it will leave.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
