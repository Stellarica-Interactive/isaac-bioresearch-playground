"""Track the loop's forward-crawling mode across a parameter sweep.

    .venv/Scripts/python tools/track_gait_mode.py --param-name g_gap_ps --values 100 50 20 5

model_assumptions 5AN.6 found that the connectome-driven loop already contains a
forward-travelling mode with the real gait's wavelength -- phase step -24.6 deg
per joint against the scripted gait's -23.1 -- and that it is damped (zeta 0.87).
The target is therefore a Hopf bifurcation of *that* mode: its damping driven to
zero, so the loop sustains the gait instead of letting it decay.

Reporting the least-damped mode overall is the wrong instrument for that, and
5AN.6 is the demonstration: at e_inh_mv = -90 a mode does go unstable, but it is a
standing wave at 0.072 Hz while the forward mode stays at 0.84. So this classifies
every oscillatory mode by its shape and follows the forward-travelling ones.

A mode counts as forward-travelling if its magnitude-weighted phase step along the
joints is within ``--band`` of the scripted gait's, with the sign calibrated on a
gait known to move head-first (5AN.6): negative is forward in the eigenvector
convention used here.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, ".")
from tools.analyse_loop import Loop, fixed_point, jacobian  # noqa: E402,F401

#: The scripted gait's phase step per joint in the eigenvector convention,
#: measured on a body known to move head-first (model_assumptions 5AN.6).
FORWARD_STEP_DEG = -23.09


def modes(loop: Loop, settle_s: float):
    """Find the loop's equilibrium, linearise there, return (rates, eigenvectors, bend)."""
    x0, jac = fixed_point(loop, settle_s=settle_s)
    mu, vec = np.linalg.eig(jac)
    lam = np.log(mu.astype(complex)) / loop.dt
    bend = float(np.degrees(np.abs(loop.joints(x0)).max()))
    return lam, vec, bend


def classify(loop: Loop, lam: np.ndarray, vec: np.ndarray):
    """Every oscillatory mode: (frequency, damping, phase step, growth rate, index)."""
    n, _, m, j = loop.sizes
    nyquist_hz = 0.5 / loop.dt
    out = []
    for k in np.flatnonzero(lam.imag > 1e-6):
        hz = lam.imag[k] / (2.0 * np.pi)
        if abs(hz - nyquist_hz) < 1.0:
            continue  # the integrator's alternating mode, not the model
        joints = vec[2 * n + m : 2 * n + m + j, k]
        mag = np.abs(joints)
        if mag.max() < 1e-12:
            continue
        step = np.diff(np.degrees(np.unwrap(np.angle(joints))))
        weight = np.minimum(mag[:-1], mag[1:])
        mean_step = float(np.sum(weight * step) / np.sum(weight))
        zeta = float(-lam.real[k] / abs(lam[k]))
        out.append((hz, zeta, mean_step, float(lam.real[k]), int(k)))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--param-name", required=True)
    parser.add_argument("--values", type=float, nargs="+", required=True)
    parser.add_argument("--settle-s", type=float, default=10.0)
    parser.add_argument(
        "--band",
        type=float,
        default=12.0,
        help="how far from the gait's phase step, in degrees, still counts as forward",
    )
    parser.add_argument(
        "--max-hz", type=float, default=5.0, help="highest frequency a gait mode can have"
    )
    parser.add_argument("--proprioceptive-mv", type=float, default=20.0)
    parser.add_argument(
        "--proprio-rate",
        type=float,
        default=0.0,
        help="fraction of proprioception that senses the rate of bending (5AK)",
    )
    parser.add_argument("--extra", action="append", default=[], metavar="NAME=VALUE")
    args = parser.parse_args()

    print(
        f"  tracking the forward-travelling mode: phase step within "
        f"{args.band:g} deg of {FORWARD_STEP_DEG:+.2f}"
    )
    print()
    header = (
        f"  {args.param_name:>12} {'bend':>7} {'fwd freq':>9} {'fwd damping':>12} "
        f"{'fwd step':>9} {'least-damped overall':>22}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))
    for value in args.values:
        ns = argparse.Namespace(
            unknown_sign="exclude",
            physics_hz=240.0,
            neural_dt_ms=1.0,
            torque_scale=3.0e-3,
            proprioceptive_mv=args.proprioceptive_mv,
            proprio_rate=args.proprio_rate,
            command_mv=20.0,
            settle_s=args.settle_s,
            open_loop=False,
            param=[f"{args.param_name}={value}", *args.extra],
        )
        loop = Loop(ns)
        try:
            lam, vec, bend = modes(loop, args.settle_s)
        except RuntimeError as err:
            print(f"  {value:12g} {err}")
            continue
        found = classify(loop, lam, vec)
        # A gait is a fraction of a hertz. With the rate term on, the 50 ms finite
        # difference amplifies fast loop modes (57 and 110 Hz) and some land in the
        # forward phase band by accident; they are not candidates for a gait.
        forward = [
            f for f in found if abs(f[2] - FORWARD_STEP_DEG) < args.band and f[0] <= args.max_hz
        ]
        overall = min(found, key=lambda f: f[1]) if found else None
        if forward:
            best = min(forward, key=lambda f: f[1])
            fwd = f"{best[0]:7.3f}Hz {best[1]:+12.4f} {best[2]:+8.2f}d"
        else:
            fwd = f"{'none':>9} {'-':>12} {'-':>9}"
        tail = f"{overall[0]:.3f}Hz z{overall[1]:+.3f} s{overall[2]:+.0f}" if overall else "-"
        print(f"  {value:12g} {bend:6.1f}d {fwd} {tail:>22}")
    print()
    print("  fwd damping crossing zero is a Hopf bifurcation of the gait itself:")
    print("  the loop would then sustain a forward wave rather than let it decay.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
