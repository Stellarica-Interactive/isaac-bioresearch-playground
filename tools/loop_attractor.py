"""What does the loop settle into, and which way does it crawl? Simulation, not linearisation.

    .venv/Scripts/python tools/loop_attractor.py --proprio-rate 1 --proprioceptive-mv 200
    .venv/Scripts/python tools/loop_attractor.py --proprio-rate 1 --proprioceptive-mv 200 \\
        --proprio-exclude VB8,VB9,VB10,VB11

model_assumptions 5AN.15a found the loop's equilibrium locally stable while the
loop, started from rest, ran to a large oscillation that crawled backward: an
attractor the linearisation could not see. So this runs the loop -- the
runner's own chain, without Isaac -- from rest and measures what it settles into
after a transient:

* ``period``: the dominant period of the joint angles.
* ``step``: the phase step from each joint to the next at that frequency,
  magnitude-weighted. Negative is a bend arriving later at the joint behind, a
  wave moving tailward -- forward crawling. The scripted gait is about -23 deg.
* ``axial``: the centroid's velocity along the body's own tail-to-head axis,
  in body lengths per second. Positive is head-first.
* ``share``: that axial displacement as a fraction of the path, so a body
  spinning in place reads near zero however far its centroid wanders.
* ``amp``: the standard deviation of the joint angles, degrees.
* ``turn``: how fast the body axis rotates, degrees per second. A forward
  crawl with a constant bend bias goes round in circles; this catches it.
* ``bias``: the time-averaged bend, mean over joints, degrees. Positive is
  dorsal (worm/body/muscles.py).
* ``clear``: the closest the body came to itself, mm. Negative means it passed
  through itself, which only ``--self-contact`` prevents.

Same body, drag and muscles as the runner's ``--quasistatic`` path, linear drag.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, ".")
from tools.analyse_loop import Loop  # noqa: E402


class _Without:
    """Proprioception with some targets' input removed -- a sensory lesion."""

    def __init__(self, inner, cells) -> None:
        self.inner, self.cells = inner, set(cells)

    def currents(self, angles, **kw):
        return {
            c: (0.0 if c in self.cells else v) for c, v in self.inner.currents(angles, **kw).items()
        }


def measure(loop: Loop, seconds: float, transient: float) -> dict[str, float]:
    every = max(1, int(round(0.05 / loop.dt)))  # 20 Hz samples
    angles, centroids, axes = [], [], []
    clear = np.inf
    for k in range(int(seconds / loop.dt)):
        loop.step()
        if k * loop.dt >= transient and k % every == 0:
            c = loop.body.segment_centres()
            angles.append(loop.body.joint_angles.copy())
            clear = min(clear, loop.body.min_self_distance_m())
            centroids.append(c.mean(axis=0))
            axis = c[0] - c[-1]
            axes.append(axis / np.linalg.norm(axis))
    a = np.array(angles)
    x = np.array(centroids)
    ax = np.array(axes)
    dt = every * loop.dt
    steps = np.diff(x, axis=0)
    axial = float(np.sum(np.einsum("ij,ij->i", steps, ax[:-1])))
    heading = np.unwrap(np.arctan2(ax[:, 1], ax[:, 0]))
    turn = float(np.degrees(heading[-1] - heading[0]) / (len(a) * dt))
    path = float(np.sum(np.linalg.norm(steps, axis=1)))
    length = float(loop.plan.total_length_m)

    centred = a - a.mean(axis=0)
    spectrum = np.fft.rfft(centred, axis=0)
    freqs = np.fft.rfftfreq(len(a), dt)
    power = (np.abs(spectrum) ** 2).sum(axis=1)
    k = int(np.argmax(power[1:]) + 1)
    x_f = spectrum[k]
    pair = x_f[1:] * np.conj(x_f[:-1])
    step = float(np.degrees(np.angle(pair.sum())))
    return {
        "period": 1.0 / freqs[k],
        "step": step,
        "axial": axial / length / (len(a) * dt),
        "share": axial / path if path > 0 else 0.0,
        "amp": float(np.degrees(centred.std(axis=0).max())),
        "turn": turn,
        "bias": float(np.degrees(a.mean())),
        "clear": 1000.0 * clear,
        # Bend amplitude joint by joint, head first: where along the body the
        # wave actually lives. A real crawl bends the whole length and starts at
        # the head.
        "profile": np.degrees(centred.std(axis=0)),
        # Each joint's phase at the dominant frequency, relative to joint 0,
        # unwrapped from head to tail. A wave starting at the head and running
        # tailward falls steadily; a head that is not part of the wave shows as a
        # jump between the head joints and the rest.
        "phase": np.degrees(np.unwrap(np.angle(x_f * np.conj(x_f[0])))),
        # Dominant period of the head (joints 0-5) and of the body (10-19)
        # separately: one shared rhythm, or two oscillators that do not lock.
        "period_head": _dominant_period(centred[:, :6], dt),
        "period_body": _dominant_period(centred[:, 10:20], dt),
    }


def _dominant_period(block: np.ndarray, dt: float) -> float:
    spectrum = np.fft.rfft(block, axis=0)
    power = (np.abs(spectrum) ** 2).sum(axis=1)
    freqs = np.fft.rfftfreq(block.shape[0], dt)
    k = int(np.argmax(power[1:]) + 1)
    return float(1.0 / freqs[k])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--proprioceptive-mv", type=float, default=20.0)
    parser.add_argument("--proprio-rate", type=float, default=0.0)
    parser.add_argument("--proprio-exclude", default="")
    parser.add_argument("--open-loop", action="store_true", help="cut proprioception: the control")
    parser.add_argument(
        "--muscle-tau-ms", type=float, default=None, help="muscle activation time constant"
    )
    parser.add_argument("--seconds", type=float, default=80.0)
    parser.add_argument("--transient", type=float, default=20.0)
    parser.add_argument("--param", action="append", default=[], metavar="NAME=VALUE")
    parser.add_argument(
        "--gap-scale",
        action="append",
        default=[],
        metavar="CLASS=FACTOR",
        help="scale one class of gap junction: neuron, muscle or neuron-muscle",
    )
    parser.add_argument(
        "--perturb-deg",
        type=float,
        default=0.0,
        help="random initial joint angles of this size, to probe for other attractors",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--muscle-gap-scale",
        type=float,
        default=None,
        help="body-wall muscle coupling; default committed, 1 as before 5AO.5",
    )
    parser.add_argument(
        "--head-proprio",
        action="store_true",
        help="SMD neurons proprioceptive, as Yeon et al. 2018 measured (5AO.7)",
    )
    parser.add_argument("--head-offset", type=float, default=None, help="segments anterior")
    parser.add_argument(
        "--head-classes",
        default="SMD",
        help="head proprioceptive classes; SMD is measured (Yeon et al. 2018), others are not",
    )
    parser.add_argument("--head-mv", type=float, default=None)
    parser.add_argument("--head-rate", type=float, default=None)
    parser.add_argument("--head-receptive", type=float, default=None, help="fraction of body")
    parser.add_argument(
        "--a-type-proprio",
        action="store_true",
        help="ASSUMED probe: VA/DA sense the body behind them, mirroring B-type",
    )
    parser.add_argument(
        "--command", default="AVBL,AVBR", help="command interneurons driven, as the runner's"
    )
    parser.add_argument(
        "--torque-scale",
        type=float,
        default=3.0e-3,
        help="peak muscle torque scale, as the runner's --torque-scale",
    )
    parser.add_argument(
        "--self-contact",
        action="store_true",
        help="the body cannot pass through itself, as the runner's --self-collision",
    )
    args = parser.parse_args()

    loop = Loop(
        argparse.Namespace(
            unknown_sign="exclude",
            physics_hz=240.0,
            neural_dt_ms=1.0,
            torque_scale=args.torque_scale,
            proprioceptive_mv=args.proprioceptive_mv,
            proprio_rate=args.proprio_rate,
            muscle_tau_ms=args.muscle_tau_ms,
            command_mv=20.0,
            settle_s=0.0,
            open_loop=args.open_loop,
            param=args.param,
            gap_scale=args.gap_scale,
            self_contact=args.self_contact,
            command_cells=args.command,
            a_type_proprio=args.a_type_proprio,
            head_proprio=args.head_proprio,
            head_offset=args.head_offset,
            head_classes=tuple(c.strip() for c in args.head_classes.split(",") if c.strip()),
            head_mv=args.head_mv,
            head_rate=args.head_rate,
            head_receptive=args.head_receptive,
            muscle_gap_scale=args.muscle_gap_scale,
        )
    )
    excluded = [c.strip() for c in args.proprio_exclude.split(",") if c.strip()]
    if excluded:
        loop.proprio = _Without(loop.proprio, excluded)
    if args.perturb_deg:
        rng = np.random.default_rng(args.seed)
        loop.body.joint_angles = np.radians(args.perturb_deg) * rng.standard_normal(loop.sizes[3])
    r = measure(loop, args.seconds, args.transient)
    print(
        f"  {'OPEN LOOP, ' if args.open_loop else ''}"
        f"rate {args.proprio_rate:g}, {args.proprioceptive_mv:g} mV/rad"
        f"{', no ' + '/'.join(excluded) if excluded else ''}"
        f"{f', muscle tau {args.muscle_tau_ms:g} ms' if args.muscle_tau_ms else ''}"
        f"{', ' + ', '.join(args.param) if args.param else ''}"
        f"{', gap ' + ', '.join(args.gap_scale) if args.gap_scale else ''}"
        f"{', self-contact' if args.self_contact else ''}"
        f"{f', command {args.command}' if args.command != 'AVBL,AVBR' else ''}"
        f"{', A-type proprio' if args.a_type_proprio else ''}"
        f"{f', head proprio {args.head_classes}' if args.head_proprio else ''}"
        f"{f' (offset {args.head_offset:g})' if args.head_offset else ''}"
        f"{f' ({args.head_mv:g} mV)' if args.head_mv else ''}"
        f"{f' (rate {args.head_rate:g})' if args.head_rate is not None else ''}"
        f"{f' (field {args.head_receptive:g})' if args.head_receptive else ''}"
        f"{f', torque {args.torque_scale:g}' if args.torque_scale != 3.0e-3 else ''}"
        f"{f', start {args.perturb_deg:g} deg seed {args.seed}' if args.perturb_deg else ''}: "
        f"period {r['period']:5.2f} s, step {r['step']:+6.1f} deg, "
        f"axial {r['axial']:+.4f} BL/s, share {r['share']:+.2f}, amp {r['amp']:5.1f} deg, "
        f"turn {r['turn']:+6.1f} deg/s, bias {r['bias']:+5.2f} deg, clear {r['clear']:+.2f} mm"
    )
    print(
        "    bend amplitude by joint, head to tail, deg: "
        + " ".join(f"{v:.1f}" for v in r["profile"])
    )
    print(
        f"    period head {r['period_head']:.2f} s, body {r['period_body']:.2f} s; "
        "phase by joint vs joint 0, every 2nd, deg: "
        + " ".join(f"{v:+.0f}" for v in r["phase"][::2])
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
