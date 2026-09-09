"""The simulator: explicit state, one deterministic timestep, checkpoint and restore.

    runtime = NeuralRuntime.build(connectome, unknown_sign=UnknownSignPolicy.EXCLUDE)
    runtime.inject("ASHL", 5.0)          # pA into a sensory neuron
    runtime.run(duration_ms=500)
    runtime.voltage("AVAL")

What one timestep does, in words
-------------------------------

Each cell holds two numbers: membrane voltage ``V`` and synaptic activation ``s``.
A step computes four currents into every cell -- leak pulling it toward rest, gap
junctions pulling it toward its electrically coupled neighbours, chemical synapses
pulling it toward each incoming synapse's reversal potential in proportion to how
active the presynaptic cell is, and whatever current we inject -- divides the total
by capacitance to get ``dV/dt``, updates ``s`` toward a sigmoid of the cell's own
voltage, and advances both by ``dt`` using the chosen integrator.

There is no event, no threshold crossing, no spike. Everything is continuous.
``docs/neural_runtime.md`` works through the arithmetic on a two-neuron network.

Determinism
-----------

Given the same connectome, parameters, integrator, timestep and inputs, this
produces bit-identical results across runs and processes. No randomness is used
anywhere; if noise is added later it must come from a seeded generator stored in
the checkpoint. The test suite asserts this.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import ClassVar

import numpy as np

from common.data.schemas import Connectome
from common.neural.neuron_models import (
    GradedLeakyIntegrator,
    GradedLeakyIntegratorParameters,
    NeuronModel,
    prepare,
)
from common.neural.synapses import (
    NetworkMatrices,
    UnknownSignPolicy,
    WeightScaling,
    build_matrices,
)

CHECKPOINT_VERSION = "1"


class Integrator(StrEnum):
    """How the state is advanced.

    ``EULER``
        One forward Euler step. Simple enough to check by hand, which is why it is
        here; it needs a small ``dt`` and is the least accurate.

    ``RK4``
        Classical fourth-order Runge-Kutta. Four derivative evaluations per step.
        The default: accurate, and its error falls as ``dt^4``.

    ``EXPONENTIAL``
        Exponential Euler. Solves each cell's own linear leak-plus-conductance term
        exactly over the step and treats the coupling explicitly, so it stays stable
        at timesteps where Euler blows up. The right choice when a hub neuron's time
        constant is far shorter than the timescale of interest.
    """

    EULER = "euler"
    RK4 = "rk4"
    EXPONENTIAL = "exponential"


#: How far past a cell's membrane time constant each integrator stays bounded.
#: Forward Euler on dV/dt = -(G/C)V is stable for dt < 2 tau; RK4 for about 2.785
#: tau; the exponential integrator solves that term exactly and has no such limit.
_STABILITY_MARGIN: dict[Integrator, float] = {
    Integrator.EULER: 2.0,
    Integrator.RK4: 2.785,
    Integrator.EXPONENTIAL: float("inf"),
}


@dataclass(frozen=True, slots=True)
class StabilityReport:
    """Whether the chosen timestep can be trusted for this network."""

    fastest_time_constant_ms: float
    dt_ms: float
    integrator: Integrator
    stability_limit_ms: float
    stable: bool
    recommended_dt_ms: float

    @property
    def dt_over_tau(self) -> float:
        return self.dt_ms / self.fastest_time_constant_ms

    def to_dict(self) -> dict[str, float | str | bool]:
        return {
            "fastest_time_constant_ms": self.fastest_time_constant_ms,
            "dt_ms": self.dt_ms,
            "integrator": str(self.integrator),
            "stability_limit_ms": self.stability_limit_ms,
            "dt_over_tau": self.dt_over_tau,
            "stable": self.stable,
            "recommended_dt_ms": self.recommended_dt_ms,
        }

    def summary(self) -> str:
        verdict = "stable" if self.stable else "UNSTABLE"
        return (
            f"dt = {self.dt_ms} ms, {self.integrator}, fastest tau = "
            f"{self.fastest_time_constant_ms:.4f} ms "
            f"(dt/tau = {self.dt_over_tau:.2f}, {verdict})"
        )


@dataclass
class ActivityRecorder:
    """Records state over time. Off by default; recording every step of a large
    network for a long run is a lot of memory."""

    cells: tuple[str, ...]
    every: int = 1
    times_ms: list[float] = field(default_factory=list)
    states: list[np.ndarray] = field(default_factory=list)

    def maybe_record(self, step: int, t_ms: float, state: np.ndarray, idx: np.ndarray) -> None:
        if step % self.every == 0:
            self.times_ms.append(t_ms)
            self.states.append(state[:, idx].copy())

    def as_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        """``(times, states)`` with states shaped ``(n_samples, n_vars, n_cells)``."""
        if not self.states:
            return np.zeros(0), np.zeros((0, 0, 0))
        return np.asarray(self.times_ms), np.stack(self.states)


@dataclass
class NeuralRuntime:
    """A network, its state, and the means to advance it."""

    network: NetworkMatrices
    model: NeuronModel
    dt_ms: float = 0.1
    integrator: Integrator = Integrator.RK4

    state: np.ndarray = field(init=False)
    i_ext_pa: np.ndarray = field(init=False)
    t_ms: float = field(default=0.0, init=False)
    step_count: int = field(default=0, init=False)
    recorder: ActivityRecorder | None = field(default=None, init=False)

    def __post_init__(self) -> None:
        self.state = self.model.initial_state(self.network)
        self.i_ext_pa = np.zeros(self.network.n, dtype=np.float64)

    # -- construction ------------------------------------------------------

    @classmethod
    def build(
        cls,
        connectome: Connectome,
        *,
        unknown_sign: UnknownSignPolicy,
        params: GradedLeakyIntegratorParameters | None = None,
        weight_scaling: WeightScaling = WeightScaling.LINEAR,
        cells: Sequence[str] | None = None,
        dt_ms: float = 0.1,
        integrator: Integrator = Integrator.RK4,
    ) -> NeuralRuntime:
        """Build a runtime from a connectome.

        ``unknown_sign`` is required. Around half of chemical connections have no
        predicted polarity and every neuromuscular junction is unsigned, so there is
        no honest default -- see :class:`~common.neural.synapses.UnknownSignPolicy`.
        """
        p = params or GradedLeakyIntegratorParameters()
        network = build_matrices(
            connectome,
            unknown_sign=unknown_sign,
            g_syn_ps=p.g_syn_ps,
            g_gap_ps=p.g_gap_ps,
            e_exc_mv=p.e_exc_mv,
            e_inh_mv=p.e_inh_mv,
            e_leak_mv=p.e_leak_mv,
            weight_scaling=weight_scaling,
            cells=tuple(cells) if cells is not None else None,
        )
        model = prepare(GradedLeakyIntegrator(params=p), network)
        return cls(network=network, model=model, dt_ms=dt_ms, integrator=integrator)

    # -- inputs and outputs ------------------------------------------------

    def inject(self, cell_id: str, current_pa: float) -> None:
        """Hold a constant current into one cell. This is the sensory input path."""
        self.i_ext_pa[self.network.index(cell_id)] = current_pa

    def inject_many(self, currents: Mapping[str, float]) -> None:
        for cell_id, value in currents.items():
            self.inject(cell_id, value)

    def clear_inputs(self) -> None:
        self.i_ext_pa[:] = 0.0

    def voltage(self, cell_id: str) -> float:
        return float(self.state[0, self.network.index(cell_id)])

    def activation(self, cell_id: str) -> float:
        """Synaptic activation: how strongly this cell is driving its targets."""
        return float(self.state[1, self.network.index(cell_id)])

    def voltages(self) -> dict[str, float]:
        return dict(zip(self.network.cell_ids, self.state[0].tolist(), strict=True))

    def activations(self, cells: Sequence[str] | None = None) -> dict[str, float]:
        ids = tuple(cells) if cells is not None else self.network.cell_ids
        return {c: self.activation(c) for c in ids}

    def current_breakdown(self, cell_id: str) -> dict[str, float]:
        """The four current terms for one cell, in pA. For understanding *why* it moved."""
        terms = self.model.currents(  # type: ignore[attr-defined]
            self.state[0], self.state[1], self.i_ext_pa, self.network
        )
        i = self.network.index(cell_id)
        return {k: float(v[i]) for k, v in terms.items()}

    def record(self, cells: Sequence[str] | None = None, every: int = 1) -> ActivityRecorder:
        ids = tuple(cells) if cells is not None else self.network.cell_ids
        self.recorder = ActivityRecorder(cells=ids, every=every)
        return self.recorder

    # -- integration -------------------------------------------------------

    def _derivs(self, state: np.ndarray) -> np.ndarray:
        return self.model.derivatives(state, self.i_ext_pa, self.network)

    def _step_euler(self, dt: float) -> np.ndarray:
        return self.state + dt * self._derivs(self.state)

    def _step_rk4(self, dt: float) -> np.ndarray:
        k1 = self._derivs(self.state)
        k2 = self._derivs(self.state + 0.5 * dt * k1)
        k3 = self._derivs(self.state + 0.5 * dt * k2)
        k4 = self._derivs(self.state + dt * k3)
        return self.state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    def _step_exponential(self, dt: float) -> np.ndarray:
        """Solve each cell's own linear term exactly; treat coupling explicitly.

        For ``dV/dt = (I_total - G_total V) / C`` with the neighbours held fixed
        across the step, the exact solution is a relaxation toward ``I/G`` with time
        constant ``C/G``. That removes the stiffness the leak-plus-conductance term
        would otherwise impose, so large ``dt`` stays bounded instead of diverging.
        """
        model = self.model
        v, s = self.state[0], self.state[1]
        g_total = model.total_conductance(s, self.network)  # type: ignore[attr-defined]
        terms = model.currents(v, s, self.i_ext_pa, self.network)  # type: ignore[attr-defined]
        # Total current with the V-dependence of this cell's own terms removed.
        i_drive = (
            terms["leak"] + terms["gap"] + terms["chemical"] + terms["external"]
        ) + g_total * v
        v_inf = i_drive / g_total
        decay = np.exp(-dt * g_total / model.params.c_m_pf)  # type: ignore[attr-defined]
        v_next = v_inf + (v - v_inf) * decay

        ds = self._derivs(self.state)[1]
        s_next = np.clip(s + dt * ds, 0.0, 1.0)
        return np.vstack([v_next, s_next])

    # ClassVar, not a field: a dataclass would otherwise treat this as per-instance state.
    _STEPPERS: ClassVar[dict[Integrator, str]] = {
        Integrator.EULER: "_step_euler",
        Integrator.RK4: "_step_rk4",
        Integrator.EXPONENTIAL: "_step_exponential",
    }

    def step(self, dt_ms: float | None = None) -> None:
        """Advance one timestep."""
        dt = self.dt_ms if dt_ms is None else dt_ms
        stepper: Callable[[float], np.ndarray] = getattr(self, self._STEPPERS[self.integrator])
        nxt = stepper(dt)
        if not np.all(np.isfinite(nxt)):
            raise FloatingPointError(
                f"state became non-finite at t={self.t_ms:.3f} ms (step {self.step_count}). "
                f"dt={dt} ms is too large for this network: the fastest membrane time "
                f"constant is {self.fastest_time_constant_ms():.4f} ms. Reduce dt, or use "
                f"Integrator.EXPONENTIAL."
            )
        self.state = nxt
        self.t_ms += dt
        self.step_count += 1
        if self.recorder is not None:
            idx = np.array([self.network.index(c) for c in self.recorder.cells])
            self.recorder.maybe_record(self.step_count, self.t_ms, self.state, idx)

    def run(self, duration_ms: float, dt_ms: float | None = None) -> None:
        dt = self.dt_ms if dt_ms is None else dt_ms
        n_steps = int(round(duration_ms / dt))
        for _ in range(n_steps):
            self.step(dt)

    # -- stability ---------------------------------------------------------

    def fastest_time_constant_ms(self) -> float:
        return self.model.fastest_time_constant_ms(self.network)  # type: ignore[attr-defined]

    def stability_report(self) -> StabilityReport:
        """Whether the timestep is small enough, and why.

        Forward Euler on ``dV/dt = -(G/C) V`` is stable only for ``dt < 2 tau``, and
        needs to be far smaller than that to be accurate. RK4 tolerates roughly
        ``dt < 2.8 tau``. The exponential integrator has no such limit for the leak
        term, though gap-junction coupling is still explicit.
        """
        tau = self.fastest_time_constant_ms()
        limit = _STABILITY_MARGIN[self.integrator] * tau
        return StabilityReport(
            fastest_time_constant_ms=tau,
            dt_ms=self.dt_ms,
            integrator=self.integrator,
            stability_limit_ms=limit,
            stable=bool(self.dt_ms < limit),
            recommended_dt_ms=min(self.dt_ms, tau / 10.0),
        )

    # -- checkpointing -----------------------------------------------------

    def network_fingerprint(self) -> str:
        """Hash of the wiring and parameters, so a checkpoint cannot be restored
        into a different network and quietly produce nonsense."""
        h = hashlib.sha256()
        h.update("|".join(self.network.cell_ids).encode())
        h.update(np.ascontiguousarray(self.network.g_gap).tobytes())
        h.update(np.ascontiguousarray(self.network.g_syn).tobytes())
        h.update(np.ascontiguousarray(self.network.e_rev).tobytes())
        h.update(str(self.model.params).encode())  # type: ignore[attr-defined]
        return h.hexdigest()

    def save_checkpoint(self, path: Path) -> None:
        """Write the complete simulated state.

        Complete means complete: every state variable of every cell, the clock, the
        step count, and the standing external input, plus a fingerprint of the
        network they belong to. Restoring gives a numerically identical continuation.

        This is continuity of a *numerical simulation*. It is not preservation of a
        mind, and the connectome it runs on contains no memories to preserve.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        meta = {
            "checkpoint_version": CHECKPOINT_VERSION,
            "network_fingerprint": self.network_fingerprint(),
            "cell_ids": list(self.network.cell_ids),
            "state_names": list(self.model.state_names),
            "t_ms": self.t_ms,
            "step_count": self.step_count,
            "dt_ms": self.dt_ms,
            "integrator": str(self.integrator),
            "unknown_sign_policy": str(self.network.unknown_sign_policy),
            "weight_scaling": str(self.network.weight_scaling),
        }
        np.savez_compressed(
            path,
            state=self.state,
            i_ext_pa=self.i_ext_pa,
            v_threshold_mv=self.model.v_threshold_mv,  # type: ignore[attr-defined]
            meta=json.dumps(meta, sort_keys=True),
        )

    def load_checkpoint(self, path: Path, *, allow_different_network: bool = False) -> None:
        data = np.load(Path(path), allow_pickle=False)
        meta = json.loads(str(data["meta"]))

        if meta["checkpoint_version"] != CHECKPOINT_VERSION:
            raise ValueError(
                f"checkpoint version {meta['checkpoint_version']!r}, this build expects "
                f"{CHECKPOINT_VERSION!r}"
            )
        if list(self.network.cell_ids) != meta["cell_ids"]:
            raise ValueError("checkpoint was taken on a different set of cells")
        if not allow_different_network:
            got = self.network_fingerprint()
            if got != meta["network_fingerprint"]:
                raise ValueError(
                    "checkpoint was taken on a different network (wiring or parameters "
                    f"differ): fingerprint {meta['network_fingerprint'][:12]}... vs "
                    f"{got[:12]}...\nPass allow_different_network=True only if you are "
                    "deliberately restoring a neural state into an altered network, such "
                    "as a lesion experiment."
                )

        self.state = np.asarray(data["state"], dtype=np.float64)
        self.i_ext_pa = np.asarray(data["i_ext_pa"], dtype=np.float64)
        self.t_ms = float(meta["t_ms"])
        self.step_count = int(meta["step_count"])

    # -- description -------------------------------------------------------

    def describe(self) -> str:
        return (
            f"{self.network.describe()}\n"
            f"  {self.stability_report().summary()}\n"
            f"  t = {self.t_ms:.1f} ms after {self.step_count} steps"
        )
