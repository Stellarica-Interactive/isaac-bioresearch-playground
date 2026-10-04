"""Can this network oscillate at all? Ask the eigenvalues.

    .venv/Scripts/python tools/analyse_dynamics.py

`tools/analyse_modes.py` showed the loop's steady-state response has a single
dominant mode, and model_assumptions 5AM.4a sharpened what that means: a
travelling wave is ``sin(kx)cos(wt) - cos(kx)sin(wt)``, two temporal components
in quadrature, and this model only ever produces one. But a steady-state
transfer function contains no time, so it cannot say whether a second temporal
component is *possible* -- only that one is not present at rest.

The eigenvalues of the linearised dynamics can. Oscillation of any kind needs
**complex-conjugate** eigenvalues: a real eigenvalue is a pure relaxation toward
or away from the fixed point, and a system with only real eigenvalues can decay,
latch or run away but cannot ring. So:

* **no complex eigenvalues** -- the network cannot oscillate, under any input,
  and no amount of tuning the things twelve hypotheses tuned would ever produce
  a wave. The missing ingredient is structural.
* **complex eigenvalues, all heavily damped** -- it can ring but dies in a cycle
  or two, which is what §5P.3's hand-over looks like.
* **a lightly damped complex pair near the gait frequency** -- the network has a
  latent oscillator and something is failing to excite or sustain it.

The Jacobian is taken by finite differences on `GradedLeakyIntegrator.derivatives`
rather than derived by hand, so it matches the code by construction. The state is
``[V; s]`` for every cell, so the matrix is 2n by 2n.

Linearised at rest, which is a real limitation and is stated: a driven network
sits elsewhere. See model_assumptions 5AN.
"""

from __future__ import annotations

import argparse

import numpy as np

from common.data.schemas import CellCategory
from common.neural.synapses import UnknownSignPolicy
from worm.loader import load
from worm.neural.config import RUNTIME_OVERLAYS, build_runtime

#: Gait frequency used throughout, Hz -- the frequency a latent oscillator would
#: have to sit near to matter.
GAIT_HZ = 0.5


def jacobian(runtime, step: float = 1.0e-4) -> np.ndarray:
    """Finite-difference Jacobian of the full state derivative at the current state.

    Central differences, column by column. ``step`` is relative to each state
    variable's scale: voltages are tens of mV and activations are fractions, so a
    single absolute step would be either too coarse for one or lost in rounding
    for the other.
    """
    model, network = runtime.model, runtime.network
    state = runtime.state.copy()
    i_ext = np.zeros(network.n)
    flat = state.reshape(-1)
    size = flat.size
    out = np.zeros((size, size))
    for k in range(size):
        delta = step * max(1.0, abs(flat[k]))
        plus = flat.copy()
        minus = flat.copy()
        plus[k] += delta
        minus[k] -= delta
        f_plus = model.derivatives(plus.reshape(state.shape), i_ext, network).reshape(-1)
        f_minus = model.derivatives(minus.reshape(state.shape), i_ext, network).reshape(-1)
        out[:, k] = (f_plus - f_minus) / (2.0 * delta)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--dataset", default="cook_2019_herm")
    parser.add_argument("--unknown-sign", default="exclude")
    parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="NAME=VALUE",
        help="biophysical override, as in the runner",
    )
    args = parser.parse_args()

    overrides = {}
    for item in args.param:
        name, _, value = item.partition("=")
        overrides[name.strip()] = float(value)

    connectome, _ = load(args.dataset, annotations=RUNTIME_OVERLAYS)
    muscles = {c.id for c in connectome.cells if c.category is CellCategory.MUSCLE}
    cells = tuple(
        c.id for c in connectome.cells if c.category is CellCategory.NEURON or c.id in muscles
    )
    runtime, report = build_runtime(
        args.dataset,
        unknown_sign=UnknownSignPolicy(args.unknown_sign),
        cells=cells,
        connectome=connectome,
        parameter_overrides=overrides or None,
    )
    print(report.summary())
    n = runtime.network.n
    print(f"\n  state is [V; s] for {n} cells: a {2 * n} by {2 * n} Jacobian")

    jac = jacobian(runtime)
    eig = np.linalg.eigvals(jac)  # per ms

    real = eig.real
    imag = eig.imag
    complex_mask = np.abs(imag) > 1e-9 * np.max(np.abs(eig))
    n_complex_pairs = int(np.count_nonzero(complex_mask & (imag > 0)))
    unstable = int(np.count_nonzero(real > 0))

    print(f"\n  eigenvalues: {eig.size}")
    print(f"  complex-conjugate pairs: {n_complex_pairs}")
    print(f"  with positive real part (unstable): {unstable}")
    print(f"  slowest decay: {1.0 / -real[real < 0].max():.1f} ms" if np.any(real < 0) else "")

    if n_complex_pairs == 0:
        print()
        print("  -> NO complex eigenvalues. The linearised network can only relax:")
        print("     it cannot ring, under any input, at any gain. A travelling wave")
        print("     needs oscillation, and oscillation needs complex eigenvalues.")
        return 0

    # The pairs that matter: least damped, and nearest the gait frequency.
    pairs = eig[complex_mask & (imag > 0)]
    freq_hz = pairs.imag / (2.0 * np.pi) * 1000.0  # per ms -> Hz
    # Damping ratio zeta = -Re / |lambda|; 1 is critically damped, 0 undamped.
    zeta = -pairs.real / np.abs(pairs)
    cycles_to_decay = np.where(pairs.real < 0, (pairs.imag / (2.0 * np.pi)) / -pairs.real, np.inf)

    order = np.argsort(zeta)
    print()
    print("  least-damped complex pairs:")
    print(f"  {'freq':>10} {'damping ratio':>14} {'cycles to 1/e':>14}")
    print("  " + "-" * 42)
    for k in order[:8]:
        print(f"  {freq_hz[k]:8.3f}Hz {zeta[k]:14.4f} {cycles_to_decay[k]:14.3f}")

    near = np.argmin(np.abs(freq_hz - GAIT_HZ))
    print()
    print(
        f"  pair nearest the {GAIT_HZ} Hz gait: {freq_hz[near]:.3f} Hz, damping ratio "
        f"{zeta[near]:.4f}, decays in {cycles_to_decay[near]:.3f} cycles"
    )
    print()
    if np.min(zeta) > 0.7:
        print("  -> every oscillation is heavily overdamped: the network can ring")
        print("     in principle but dies within a fraction of a cycle. That is a")
        print("     relaxation system with complex arithmetic, not an oscillator.")
    else:
        print("  -> a lightly damped mode exists. Find what excites it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
