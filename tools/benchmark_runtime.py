"""Benchmark the neural runtime, and check a result against its own assumptions.

    python tools/benchmark_runtime.py speed
    python tools/benchmark_runtime.py accuracy
    python tools/benchmark_runtime.py sweep --cell AVAL --stimulate ASHL

``speed`` measures throughput per integrator and timestep.

``accuracy`` compares each integrator against a fine reference, so a chosen
timestep can be justified rather than guessed.

``sweep`` is the important one. Nine of the ten biophysical parameters are
assumptions, and every unknown-sign policy changes the network. This varies them
one at a time and reports how much the answer moves. A result that survives the
sweep is worth something; one that does not is a statement about our parameter
choices, not about the worm.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass

import numpy as np

from common.neural.runtime import Integrator, NeuralRuntime
from common.neural.synapses import UnknownSignPolicy, WeightScaling
from worm.neural.config import build_runtime, parameter_provenance

DEFAULT_DATASET = "cook_2019_herm"


@dataclass(frozen=True)
class Timing:
    integrator: str
    dt_ms: float
    steps: int
    seconds: float

    @property
    def steps_per_second(self) -> float:
        return self.steps / self.seconds if self.seconds else float("inf")

    @property
    def seconds_per_simulated_second(self) -> float:
        return self.seconds / (self.steps * self.dt_ms / 1000.0)


def _run(
    dataset: str, integrator: Integrator, dt: float, duration_ms: float
) -> tuple[NeuralRuntime | None, float]:
    """Run one configuration. Returns ``(runtime, seconds)``, or ``(None, seconds)``
    if the integration diverged -- which is a result worth printing, not a crash."""
    runtime, _ = build_runtime(
        dataset, unknown_sign=UnknownSignPolicy.EXCLUDE, dt_ms=dt, integrator=integrator
    )
    runtime.inject_many({"ASHL": 20.0, "ASHR": 20.0})
    start = time.perf_counter()
    try:
        with np.errstate(over="ignore", invalid="ignore"):
            runtime.run(duration_ms)
    except FloatingPointError:
        return None, time.perf_counter() - start
    return runtime, time.perf_counter() - start


def cmd_speed(args: argparse.Namespace) -> int:
    rows: list[Timing] = []
    for integrator, dt in [
        (Integrator.EULER, 0.005),
        (Integrator.RK4, 0.005),
        (Integrator.RK4, 0.02),
        (Integrator.EXPONENTIAL, 0.02),
        (Integrator.EXPONENTIAL, 0.1),
        (Integrator.EXPONENTIAL, 0.5),
    ]:
        runtime, seconds = _run(args.dataset, integrator, dt, args.duration_ms)
        if runtime is None:
            print(f"{str(integrator):13s} dt={dt:<7.3f} DIVERGED", file=sys.stderr)
            continue
        rows.append(Timing(str(integrator), dt, runtime.step_count, seconds))

    if args.format == "json":
        print(json.dumps([r.__dict__ | {"steps_per_second": r.steps_per_second} for r in rows],
                         indent=2))
        return 0

    print(f"{args.dataset}, {args.duration_ms:.0f} ms simulated\n")
    print(f"{'integrator':13s} {'dt (ms)':>8s} {'steps':>8s} {'wall (s)':>9s} "
          f"{'steps/s':>10s} {'s per sim s':>12s}")
    print("-" * 66)
    for r in rows:
        print(f"{r.integrator:13s} {r.dt_ms:8.3f} {r.steps:8d} {r.seconds:9.2f} "
              f"{r.steps_per_second:10,.0f} {r.seconds_per_simulated_second:12.2f}")
    return 0


def cmd_accuracy(args: argparse.Namespace) -> int:
    """How wrong is each integrator at each timestep, against a fine RK4 reference."""
    reference, _ = _run(args.dataset, Integrator.RK4, 0.002, args.duration_ms)
    if reference is None:
        print("the reference integration itself diverged", file=sys.stderr)
        return 1
    ref_v = reference.state[0].copy()

    print(f"{args.dataset}, {args.duration_ms:.0f} ms, vs RK4 at dt = 0.002 ms\n")
    print(f"{'integrator':13s} {'dt (ms)':>8s} {'max |dV| (mV)':>15s} {'wall (s)':>9s}")
    print("-" * 50)
    for integrator, dt in [
        (Integrator.RK4, 0.005),
        (Integrator.RK4, 0.02),
        (Integrator.EXPONENTIAL, 0.02),
        (Integrator.EXPONENTIAL, 0.1),
        (Integrator.EXPONENTIAL, 0.5),
        (Integrator.EXPONENTIAL, 1.0),
    ]:
        runtime, seconds = _run(args.dataset, integrator, dt, args.duration_ms)
        if runtime is None:
            print(f"{str(integrator):13s} {dt:8.3f} {'DIVERGED':>15s} {seconds:9.2f}")
            continue
        err = float(np.max(np.abs(runtime.state[0] - ref_v)))
        print(f"{str(integrator):13s} {dt:8.3f} {err:15.6f} {seconds:9.2f}")
    print(
        "\nThe stability limit here is about 0.021 ms for RK4: the busiest cells have "
        "\nmembrane time constants near 0.008 ms because gap-junction conductance adds up."
    )
    return 0


def cmd_sweep(args: argparse.Namespace) -> int:
    """Vary each assumed parameter and report how far the answer moves."""
    stimulus = {cell: args.current_pa for cell in args.stimulate}

    def measure(
        *,
        unknown_sign: UnknownSignPolicy = UnknownSignPolicy.EXCLUDE,
        scaling: WeightScaling = WeightScaling.LINEAR,
        overrides: dict[str, float] | None = None,
    ) -> float:
        runtime, _ = build_runtime(
            args.dataset,
            unknown_sign=unknown_sign,
            weight_scaling=scaling,
            parameter_overrides=overrides,
        )
        rest = runtime.voltage(args.cell)
        runtime.inject_many(stimulus)
        runtime.run(args.duration_ms)
        return runtime.voltage(args.cell) - rest

    baseline = measure()
    print(
        f"{args.dataset}: response of {args.cell} to {args.current_pa} pA into "
        f"{', '.join(args.stimulate)} for {args.duration_ms:.0f} ms\n"
    )
    print(f"baseline (exclude unsigned, linear scaling): {baseline:+.4f} mV\n")
    print(f"{'variation':44s} {'response':>10s} {'vs base':>10s}")
    print("-" * 68)

    def report(label: str, value: float) -> None:
        delta = value - baseline
        rel = f"{delta / baseline:+.1%}" if abs(baseline) > 1e-9 else "n/a"
        print(f"{label:44s} {value:+10.4f} {rel:>10s}")

    for policy in UnknownSignPolicy:
        if policy is UnknownSignPolicy.EXCLUDE:
            continue
        report(f"unknown-sign policy = {policy}", measure(unknown_sign=policy))
    for scaling in (WeightScaling.SQRT, WeightScaling.BINARY):
        report(f"weight scaling = {scaling}", measure(scaling=scaling))

    for p in parameter_provenance():
        if not p.is_assumed:
            continue
        if abs(p.value) < 1e-12:
            # A multiplicative sweep does nothing to a parameter that is zero
            # (e_exc_mv is 0 mV), so perturb it additively instead.
            variants = [(f"{p.name} = {v:+g} {p.unit}", v) for v in (-10.0, 10.0)]
        else:
            variants = [
                (f"{p.name} x {f} (= {p.value * f:g} {p.unit})", p.value * f) for f in (0.5, 2.0)
            ]
        for label, value in variants:
            report(label, measure(overrides={p.name: value}))

    print(
        "\nA result that moves a lot here is a statement about our assumptions, "
        "not about the animal.\nSee worm/neural/parameters.toml and "
        "docs/model_assumptions.md."
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--dataset", default=DEFAULT_DATASET)
    p.add_argument("--format", default="table", choices=["table", "json"])
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("speed")
    s.add_argument("--duration-ms", type=float, default=100.0)
    s.set_defaults(func=cmd_speed)

    a = sub.add_parser("accuracy")
    a.add_argument("--duration-ms", type=float, default=50.0)
    a.set_defaults(func=cmd_accuracy)

    w = sub.add_parser("sweep")
    w.add_argument("--cell", default="AVAL", help="cell whose response is measured")
    w.add_argument("--stimulate", nargs="+", default=["ASHL", "ASHR"])
    w.add_argument("--current-pa", type=float, default=20.0)
    w.add_argument("--duration-ms", type=float, default=200.0)
    w.set_defaults(func=cmd_sweep)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
