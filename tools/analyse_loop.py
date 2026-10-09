"""Can the closed sensorimotor loop oscillate? Eigenvalues of the real loop.

    .venv/Scripts/python tools/analyse_loop.py

`tools/analyse_dynamics.py` shows the neural network on its own is a pure
relaxation system: every complex eigenvalue it has is critically damped to four
decimal places, so it cannot ring. But a travelling wave is a property of the
*loop* -- neurons drive muscles, muscles bend the body, the body is sensed and fed
back -- and a chain of non-oscillating elements can oscillate once it is closed
on itself, provided the feedback is negative and the lags add up to enough phase.
The body supplies lag the network alone does not have.

So this linearises the closed loop exactly as the runner executes it, one physics
step at a time:

    1. read the joint angles
    2. clear inputs, inject the command drive and the proprioceptive currents
    3. step the neural network ``neural_substeps`` times
    4. step the muscles toward the neural drive
    5. step the quasi-static body under the resulting joint torques

The loop state is ``[V; s; muscle activation; joint angles]``. The body's origin
and heading are left out on purpose: proprioception reads joint angles only and
the drag law is invariant under rotating or translating the whole body, so they
cannot feed back.

The one-step map is linearised by central finite differences at the loop's
**equilibrium** with the command drive on, found by Newton's method on
``M(x) = x``, rather than at rest, since a driven loop does not sit at rest. Its
eigenvalues ``mu`` are converted to continuous rates ``lambda = log(mu) / dt``,
from which frequency and damping follow.

An earlier version linearised wherever the loop happened to be after twenty
seconds of running, on the assumption that it would have settled. It had not:
the committed loop never settles, it runs a 22-second relaxation cycle, and
twenty seconds is the top of a burst. Eigenvalues taken on a moving trajectory
are not stability exponents. :func:`fixed_point` refuses to return anything that
is not an equilibrium to ``tol``. See model_assumptions 5AN.8.

Linear drag throughout: the nonlinear law of §5AH has a fixed-point iteration
inside every step and would make the finite differences measure the iteration's
tolerance as much as the dynamics. See model_assumptions 5AN.
"""

from __future__ import annotations

import argparse
from dataclasses import replace

import numpy as np

from common.data.schemas import CellCategory
from common.neural.stimulus import scale_for_depolarisation, solve_for_depolarisation
from common.neural.synapses import UnknownSignPolicy
from worm.body.geometry import BodyPlan
from worm.body.muscles import MuscleModel, MuscleParameters
from worm.body.neural_bridge import (
    DEFAULT_SENSING_OFFSET,
    PROPRIOCEPTIVE_CLASSES,
    RATE_LAG_MS,
    MuscleDrive,
    Proprioception,
)
from worm.body.quasistatic import QuasiStaticBody
from worm.importers.naming import body_wall_muscle_ids
from worm.loader import load
from worm.neural.config import RUNTIME_OVERLAYS, build_runtime, rescale_gap_junctions

GAIT_HZ = 0.5


def parse_gap_scale(items: list[str]) -> dict[str, float]:
    """``["muscle=0.05", ...]`` to ``{"muscle": 0.05}``, as ``--gap-scale`` takes them."""
    out = {}
    for item in items:
        name, _, value = item.partition("=")
        if not value:
            raise SystemExit(f"--gap-scale expects CLASS=FACTOR, got {item!r}")
        out[name.strip()] = float(value)
    return out


class Loop:
    """The runner's closed loop, without Isaac, as a map on its state."""

    def __init__(self, args: argparse.Namespace) -> None:
        overrides = {}
        for item in args.param:
            name, _, value = item.partition("=")
            overrides[name.strip()] = float(value)

        connectome, _ = load("cook_2019_herm", annotations=RUNTIME_OVERLAYS)
        self.plan = BodyPlan()
        # The runner's cell set exactly: every neuron and the 95 body-wall muscles.
        # This used to take every muscle cell in the connectome -- 135, including
        # pharyngeal, vulval and anal muscles the runner leaves out -- so every
        # result before model_assumptions 5AN.18 was computed on 437 cells against
        # the runner's 397. See 5AN.18 for what that changed.
        muscles = body_wall_muscle_ids()
        cells = tuple(
            c.id for c in connectome.cells if c.category is CellCategory.NEURON or c.id in muscles
        )
        self.runtime, _ = build_runtime(
            "cook_2019_herm",
            unknown_sign=UnknownSignPolicy(args.unknown_sign),
            cells=cells,
            connectome=connectome,
            dt_ms=args.neural_dt_ms,
            parameter_overrides=overrides or None,
            # None is the committed, measurement-derived coupling between
            # body-wall muscle cells (parameters.toml [body_wall_muscle]); every
            # result before model_assumptions 5AO.5 was computed with 1.0.
            muscle_gap_scale=getattr(args, "muscle_gap_scale", None),
        )
        gap_scale = parse_gap_scale(getattr(args, "gap_scale", None) or [])
        if gap_scale:
            # Before anything is calibrated against the network, as --param is.
            rescale_gap_junctions(self.runtime, connectome, gap_scale)
        self.dt = 1.0 / args.physics_hz
        self.substeps = max(1, int(round(self.dt * 1000.0 / args.neural_dt_ms)))

        self.bridge = MuscleDrive.build(self.plan, self.runtime.network.cell_ids)
        proprio = Proprioception.build(connectome, self.plan, classes=PROPRIOCEPTIVE_CLASSES)
        b_type = [c for c in proprio.targets if c in self.runtime.network.cell_ids]
        self.proprio = replace(
            proprio,
            gain_pa_per_rad=scale_for_depolarisation(self.runtime, b_type, args.proprioceptive_mv),
        )
        self.proprioception = not args.open_loop

        # Rate-sensitive proprioception (model_assumptions 5AK) compares the body
        # with itself RATE_LAG_MS ago, exactly as the runner does. Those past
        # postures are state: a map that dropped them would linearise a different
        # law. Carried only when the rate term is on, so every result taken
        # without it keeps the state it was measured with.
        rate = float(getattr(args, "proprio_rate", 0.0) or 0.0)
        self.proprio = replace(self.proprio, rate_fraction=rate)

        # A probe, off by default: A-type motor neurons sensing the body *behind*
        # them, the mirror of B-type, as Gao et al. 2018 propose ("proprioception
        # may ... serve as feedbacks to regulate A-MN oscillation") and nobody has
        # measured. Same law, same target depolarisation, gain calibrated against
        # the A-type cells' own conductances. See model_assumptions 5AO.
        # Further proprioceptive populations, each with the B-type law and its
        # gain calibrated against its own cells' conductances.
        self.extra_proprio: list[Proprioception] = []

        def add_population(
            classes: tuple[str, ...],
            offset: float,
            *,
            mv: float | None = None,
            law: float | None = None,
            receptive: float | None = None,
        ) -> None:
            pop = Proprioception.build(
                connectome,
                self.plan,
                classes=classes,
                offset=offset,
                **({"receptive_fraction": receptive} if receptive is not None else {}),
            )
            cells = [c for c in pop.targets if c in self.runtime.network.cell_ids]
            self.extra_proprio.append(
                replace(
                    pop,
                    gain_pa_per_rad=scale_for_depolarisation(
                        self.runtime, cells, args.proprioceptive_mv if mv is None else mv
                    ),
                    rate_fraction=rate if law is None else law,
                )
            )

        if getattr(args, "a_type_proprio", False):
            add_population(("VA", "DA"), -DEFAULT_SENSING_OFFSET)
        # The head: SMD neurons are proprioceptive -- SMDD's calcium follows dorsal
        # head bending, SMDV's ventral, and forced bends excite each even in
        # paralysed worms (Yeon et al. 2018). Same law and receptive-field rule as
        # B-type, so the only new claim is the one that was measured. 5AO.7.
        if getattr(args, "head_proprio", False):
            add_population(
                tuple(getattr(args, "head_classes", None) or ("SMD",)),
                getattr(args, "head_offset", None) or DEFAULT_SENSING_OFFSET,
                mv=getattr(args, "head_mv", None),
                law=getattr(args, "head_rate", None),
                receptive=getattr(args, "head_receptive", None),
            )
        self.lag_steps = max(1, int(round(RATE_LAG_MS / (self.dt * 1000.0))))
        self.history: list[np.ndarray] = []

        tau = getattr(args, "muscle_tau_ms", None)
        self.muscle = MuscleModel(
            self.plan,
            MuscleParameters(
                peak_torque_scale=args.torque_scale,
                **({"activation_tau_ms": tau} if tau is not None else {}),
            ),
        )
        self.body = QuasiStaticBody(
            self.plan,
            muscle=self.muscle.params,
            self_contact=bool(getattr(args, "self_contact", False)),
        )

        self.runtime.run(2000.0)
        # The command interneurons the drive goes to, as the runner's --command.
        named = getattr(args, "command_cells", None) or "AVBL,AVBR"
        wanted = [c.strip() for c in named.split(",") if c.strip()]
        cells_in = [c for c in wanted if c in self.runtime.network.cell_ids]
        self.command = solve_for_depolarisation(self.runtime, cells_in, args.command_mv)

        n = self.runtime.network.n
        self.sizes = (n, n, self.muscle.activation.size, self.plan.n_joints)
        if rate > 0.0:
            self.history = [self.body.joint_angles.copy() for _ in range(self.lag_steps)]

        #: Currents *added* after every other input, pA. ``inject`` sets a cell's
        #: current rather than adding to it, and proprioception is injected after
        #: the command, so anything a probe puts into ``command`` for a B-type cell
        #: is silently overwritten. This is the path for probes.
        self.extra: dict[str, float] = {}

    # -- the state ---------------------------------------------------------

    def joints(self, x: np.ndarray) -> np.ndarray:
        """The joint angles within a state vector. Not ``x[-j:]``: with the rate
        term on, past postures follow them."""
        n, _, m, j = self.sizes
        return x[2 * n + m : 2 * n + m + j]

    def get(self) -> np.ndarray:
        return np.concatenate(
            [
                self.runtime.state[0],
                self.runtime.state[1],
                self.muscle.activation.reshape(-1),
                self.body.joint_angles,
                *self.history,
            ]
        )

    def put(self, x: np.ndarray) -> None:
        n, _, m, j = self.sizes
        self.runtime.state[0] = x[:n]
        self.runtime.state[1] = x[n : 2 * n]
        self.muscle.activation[:] = x[2 * n : 2 * n + m].reshape(self.muscle.activation.shape)
        self.body.joint_angles = x[2 * n + m : 2 * n + m + j].copy()
        start = 2 * n + m + j
        self.history = [
            x[start + k * j : start + (k + 1) * j].copy() for k in range(len(self.history))
        ]
        # Rotation- and translation-invariant: these cannot feed back.
        self.body.origin = np.zeros(2)
        self.body.heading = 0.0

    def step(self) -> None:
        """One physics step, in the runner's order."""
        angles = self.body.joint_angles.copy()
        self.runtime.clear_inputs()
        self.runtime.inject_many(self.command)
        if self.proprioception:
            if self.history:
                # Oldest first: history[0] is the posture lag_steps steps ago.
                currents = self.proprio.currents(
                    angles,
                    earlier_angles_rad=self.history[0],
                    elapsed_ms=self.lag_steps * self.dt * 1000.0,
                )
            else:
                currents = self.proprio.currents(angles)
            self.runtime.inject_many(currents)
            for pop in self.extra_proprio:
                if self.history:
                    self.runtime.inject_many(
                        pop.currents(
                            angles,
                            earlier_angles_rad=self.history[0],
                            elapsed_ms=self.lag_steps * self.dt * 1000.0,
                        )
                    )
                else:
                    self.runtime.inject_many(pop.currents(angles))
        if self.history:
            self.history = [*self.history[1:], angles]
        for cell, current in self.extra.items():
            self.runtime.i_ext_pa[self.runtime.network.index(cell)] += current
        for _ in range(self.substeps):
            self.runtime.step()
        self.muscle.step(self.bridge.drive(self.runtime.state[1]), dt_ms=self.dt * 1000.0)
        self.body.step(self.muscle.joint_torques(), dt_s=self.dt)

    def map(self, x: np.ndarray) -> np.ndarray:
        self.put(x)
        self.step()
        return self.get()


def jacobian(loop: Loop, x0: np.ndarray) -> np.ndarray:
    """Central-difference Jacobian of the loop's one-step map at ``x0``."""
    size = x0.size
    jac = np.zeros((size, size))
    for k in range(size):
        h = 1.0e-5 * max(1.0, abs(x0[k]))
        plus, minus = x0.copy(), x0.copy()
        plus[k] += h
        minus[k] -= h
        jac[:, k] = (loop.map(plus) - loop.map(minus)) / (2.0 * h)
    return jac


def fixed_point(
    loop: Loop,
    *,
    settle_s: float = 10.0,
    average_s: float = 30.0,
    tol: float = 1.0e-9,
    max_iter: int = 15,
    start: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """The loop's equilibrium and the Jacobian there, by damped Newton iteration.

    Starts from the trajectory's time average, which lies inside any cycle the
    loop is running and is a better guess than any single point on it. Each step
    is limited to 5 mV and 5 degrees so a poor Jacobian far from the solution
    cannot throw the iterate out of the region the model is valid in.

    The loop has more than one equilibrium -- the committed one is a C-shape,
    and scaling one pathway by 10% can land Newton on an S-shaped saddle
    instead. To follow *one* equilibrium as a parameter changes, pass the
    previous one as ``start``; the trajectory is then not used at all.

    Raises if no equilibrium is found: a linearisation anywhere else would be
    reported as if it meant something.
    """
    n, _, _, j = loop.sizes
    if start is None:
        for _ in range(int(settle_s / loop.dt)):
            loop.step()
        samples = []
        for _ in range(int(average_s / loop.dt)):
            loop.step()
            samples.append(loop.get())
        x = np.mean(samples, axis=0)
    else:
        x = start.copy()
    for _ in range(max_iter):
        jac = jacobian(loop, x)
        r = loop.map(x) - x
        if float(np.abs(r).max()) < tol:
            return x, jac
        dx = np.linalg.solve(jac - np.eye(x.size), -r)
        scale = min(
            1.0,
            5.0 / max(1e-12, float(np.abs(dx[:n]).max())),
            np.radians(5.0) / max(1e-12, float(np.abs(loop.joints(dx)).max())),
        )
        x = x + scale * dx
        x[n : 2 * n] = np.clip(x[n : 2 * n], 0.0, 1.0)
    residual = float(np.abs(loop.map(x) - x).max())
    raise RuntimeError(f"no equilibrium found: residual {residual:.2e} after {max_iter} steps")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--unknown-sign", default="exclude")
    parser.add_argument("--physics-hz", type=float, default=240.0)
    parser.add_argument("--neural-dt-ms", type=float, default=1.0)
    parser.add_argument("--torque-scale", type=float, default=3.0e-3)
    parser.add_argument("--proprioceptive-mv", type=float, default=20.0)
    parser.add_argument("--command-mv", type=float, default=20.0)
    parser.add_argument("--settle-s", type=float, default=10.0)
    parser.add_argument(
        "--open-loop",
        action="store_true",
        help="cut proprioception, so the body is driven but never sensed: the "
        "control that says what closing the loop adds",
    )
    parser.add_argument("--param", action="append", default=[], metavar="NAME=VALUE")
    args = parser.parse_args()

    loop = Loop(args)
    print(f"  loop state: {sum(loop.sizes)} variables {loop.sizes}")
    print(
        f"  {'OPEN' if args.open_loop else 'CLOSED'} loop, {loop.substeps} neural steps "
        f"per {loop.dt * 1000:.2f} ms physics step"
    )

    # The loop's equilibrium with the command on: not rest, and not wherever a
    # trajectory happens to be, which is what this used to report.
    x0, jac = fixed_point(loop, settle_s=args.settle_s)
    bend = float(np.degrees(np.abs(loop.joints(x0)).max()))
    residual = float(np.abs(loop.map(x0) - x0).max())
    print(f"  equilibrium: max joint angle {bend:.2f} deg, residual {residual:.1e}")

    mu = np.linalg.eigvals(jac)
    # Discrete multipliers to continuous rates, per second.
    lam = np.log(mu.astype(complex)) / loop.dt
    rate, omega = lam.real, lam.imag
    unstable = int(np.count_nonzero(np.abs(mu) > 1.0 + 1e-9))
    oscillating = (np.abs(omega) > 1e-6) & (omega > 0)
    print(f"\n  eigenvalues: {mu.size}, multipliers outside the unit circle: {unstable}")
    print(f"  complex pairs: {int(np.count_nonzero(oscillating))}")

    if not oscillating.any():
        print("\n  -> no oscillatory modes at all. The closed loop can only relax.")
        return 0

    freq = omega[oscillating] / (2.0 * np.pi)
    zeta = -rate[oscillating] / np.abs(lam[oscillating])
    order = np.argsort(zeta)
    print("\n  least-damped oscillatory modes of the closed loop:")
    print(f"  {'freq':>10} {'damping ratio':>14} {'growth /s':>11}")
    print("  " + "-" * 38)
    for k in order[:8]:
        print(f"  {freq[k]:8.3f}Hz {zeta[k]:14.4f} {rate[oscillating][k]:+11.4f}")

    near = int(np.argmin(np.abs(freq - GAIT_HZ)))
    print(f"\n  nearest the {GAIT_HZ} Hz gait: {freq[near]:.3f} Hz, damping {zeta[near]:.4f}")
    print()
    if zeta.min() > 0.7:
        print("  -> every oscillation is heavily overdamped, with the body in the loop.")
        print("     Closing the loop does not create an oscillator. See 5AN.")
    elif zeta.min() > 0.0:
        print("  -> a lightly damped mode exists: the loop can ring, and something")
        print("     is failing to excite or sustain it. See 5AN.")
    else:
        print("  -> an UNSTABLE oscillatory mode: the loop should self-oscillate.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
