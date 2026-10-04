"""Why is the B-type output one signal? Ask the wiring, not the simulation.

    .venv/Scripts/python tools/analyse_modes.py

Twelve hypotheses have been tested against the model's inability to produce a
travelling wave, and every one left the B-type adjacent correlation at or near
+1.00 with 100 per cent of the population's variance in a single component
(`negative_result.md`, model_assumptions 5Q-5U, 5AK, 5AL). Each changed a
parameter or added a mechanism. None asked whether the *structure* permits
anything else.

This does. The graded leaky integrator's voltage equation is

    C dV/dt = -G_leak (V - E_leak) + (Ggap V - V . gap_row_sum)
              + sum_j Gsyn[i,j] s_j (E[i,j] - V_i) + I_ext

so with the synaptic activations held at their resting values, the steady
response to an injected current pattern is ``dV = A^-1 dI`` where

    A = diag(G_leak + gap_row_sum + s* . Gsyn_row_sum) - Ggap

which is the matrix `neuron_models` already solves for the resting potentials.
Proprioception injects into the B-type motor neurons and the body reads their
activations back, so **the block of ``A^-1`` from B-type inputs to B-type
voltages is the linear transfer function of the whole sensorimotor loop**.

Two things come out of it. Its singular values say how many independent things
that loop can do: if the first dominates absolutely, no gain, sign policy or
proprioceptive law could ever have produced a wave. And the spatial shape of
those modes says whether any of them *is* a wave -- a mode whose entries all
share a sign is a whole-body contraction and cannot propagate however strongly
it is driven.

Linear response is a real limitation and is stated rather than hidden: the model
is nonlinear, the sigmoid can gate cells off, and 5L showed the resting state
matters. What this bounds is what the *linearisation* can carry, which is where
every measured failure so far has lived. See model_assumptions 5AM.
"""

from __future__ import annotations

import argparse

import numpy as np

from common.data.schemas import CellCategory
from common.neural.synapses import UnknownSignPolicy
from worm.body.geometry import BodyPlan
from worm.body.neural_bridge import PROPRIOCEPTIVE_CLASSES, Proprioception
from worm.loader import load
from worm.neural.config import RUNTIME_OVERLAYS, build_runtime


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--dataset", default="cook_2019_herm")
    parser.add_argument(
        "--unknown-sign",
        default="exclude",
        help="what to do with the 1846 chemical connections of unknown sign",
    )
    parser.add_argument("--top", type=int, default=8, help="how many singular values to print")
    parser.add_argument(
        "--loop",
        choices=("open", "closed"),
        default="closed",
        help="open holds the synaptic activations at rest, which answers what "
        "the wiring permits. closed includes the chemical feedback through "
        "ds*/dV, which is the actual sensorimotor loop and the right object for "
        "asking why its output is one signal.",
    )
    args = parser.parse_args()

    # load returns (connectome, overlay report); the runtime wants the first.
    connectome, _ = load(args.dataset, annotations=RUNTIME_OVERLAYS)
    plan = BodyPlan()
    muscles = {c.id for c in connectome.cells if c.category is CellCategory.MUSCLE}
    cells = tuple(
        c.id for c in connectome.cells if c.category is CellCategory.NEURON or c.id in muscles
    )
    runtime, report = build_runtime(
        args.dataset,
        unknown_sign=UnknownSignPolicy(args.unknown_sign),
        cells=cells,
        connectome=connectome,
    )
    print(report.summary())

    network = runtime.network
    params = runtime.model.params
    # The same matrix neuron_models solves the resting potentials with: the
    # conductance the network presents to a voltage perturbation, with synaptic
    # activations held at rest.
    diag = params.g_leak_ns + network.gap_row_sum + params.s_at_rest * network.g_syn.sum(axis=1)
    a = np.diag(diag) - network.g_gap

    # The chemical feedback the open-loop version omits.
    #
    # The sigmoid sits `v_threshold_offset_mv` above each cell's OWN resting
    # potential, so at rest every cell's argument is -offset whatever its resting
    # potential is (see `s_at_rest`). Every cell therefore has the same feedback
    # gain, and the loop has no cell-to-cell variation in how strongly it can
    # recruit anybody -- which is itself a finding about the parameterisation.
    rest = runtime.state[0]
    thresholds = runtime.model._thresholds(network)
    z = params.beta_per_mv * (rest - thresholds)
    phi = 1.0 / (1.0 + np.exp(-z))
    dphi = params.beta_per_mv * phi * (1.0 - phi)
    denom = params.a_d_per_ms + params.a_r_per_ms * phi
    # s* = a_r phi / (a_d + a_r phi), so ds*/dphi = a_r a_d / denom^2.
    ds_dv = params.a_r_per_ms * params.a_d_per_ms / denom**2 * dphi

    feedback = network.g_syn * (network.e_rev - rest[:, None]) * ds_dv[None, :]
    matrix = a if args.loop == "open" else a - feedback
    inverse = np.linalg.inv(matrix)

    print()
    if args.loop == "closed":
        print("  transfer function: CLOSED loop, chemical feedback included")
        spread = float(ds_dv.max() - ds_dv.min())
        print(f"  ds*/dV at rest: {ds_dv.mean():.4e} per mV, spread {spread:.2e}")
        print(f"  feedback block: {int(np.count_nonzero(feedback))} non-zero entries")
    else:
        print("  transfer function: OPEN loop, activations held at rest")

    proprio = Proprioception.build(connectome, plan, classes=PROPRIOCEPTIVE_CLASSES)
    targets = [c for c in proprio.targets if c in network.cell_ids]
    index = np.array([network.index(c) for c in targets], dtype=np.int64)
    print()
    print(f"  {len(targets)} proprioceptive targets in the network: {targets[:6]} ...")

    block = inverse[np.ix_(index, index)]
    left, singular, _ = np.linalg.svd(block)
    total = float(np.sum(singular**2))

    print()
    print("  singular values of the B-type to B-type transfer block:")
    print(f"  {'i':>3} {'sigma':>12} {'sigma/sigma_1':>14} {'cumulative var':>15}")
    print("  " + "-" * 48)
    running = 0.0
    for i, value in enumerate(singular[: args.top]):
        running += float(value**2)
        print(
            f"  {i + 1:3d} {value:12.4e} {value / singular[0]:14.4f} {running / total * 100:14.2f}%"
        )

    first = float(singular[0] ** 2 / total)
    print()
    print(f"  the first mode carries {first * 100:.2f}% of the response")
    print(f"  condition number sigma_1/sigma_last: {singular[0] / singular[-1]:.1f}")
    print()
    if first > 0.99:
        print("  -> effectively rank one. Every input pattern produces the same")
        print("     output pattern scaled, so no proprioceptive law, gain or sign")
        print("     policy could have produced a wave through this block.")
    else:
        print("  -> the wiring admits several independent modes, so a single")
        print("     signal is not forced by the connectivity.")

    # A mode whose entries all share a sign is a whole-body contraction; one that
    # alternates along the body is a spatial wave. That distinction decides
    # whether any input could produce propagation.
    order = np.argsort([int(c[2:]) for c in targets])
    print()
    print("  spatial structure of the leading modes, head to tail:")
    print("  " + " ".join(f"{targets[i]:>5}" for i in order))
    for k in range(min(3, left.shape[1])):
        vector = left[:, k][order]
        vector = vector / np.abs(vector).max()
        flips = int(np.count_nonzero(np.diff(np.sign(vector))))
        row = " ".join(f"{v:+5.2f}" for v in vector)
        print(f"  {row}   mode {k + 1}, {flips} sign changes")
    print()
    print("  Zero sign changes is a whole-body contraction and cannot propagate.")
    print("  A travelling wave needs alternating modes AND an input exciting them.")

    # Does the input the body actually produces excite anything but mode 1?
    position = np.array(
        [proprio.sensed_segment[proprio.targets.index(c)] for c in targets], dtype=np.float64
    )
    sign = np.array([1.0 if c.startswith("DB") else -1.0 for c in targets])
    print()
    print("  how much of each input pattern lands in mode 1:")
    for name, pattern in (
        ("uniform bend, what a latch gives", sign * np.ones(len(targets))),
        ("one wavelength on the body", sign * np.sin(2.0 * np.pi * position / plan.n_joints)),
        ("1.5 wavelengths, the real gait", sign * np.sin(3.0 * np.pi * position / plan.n_joints)),
    ):
        unit = pattern / np.linalg.norm(pattern)
        print(f"    {name:34s} {float((unit @ left[:, 0]) ** 2) * 100:6.2f}%")

    print()
    print("  Linear response only: the model is nonlinear and the sigmoid can gate")
    print("  cells off. This bounds what the linearisation can carry, which is")
    print("  where every measured failure so far has lived. See 5AM.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
