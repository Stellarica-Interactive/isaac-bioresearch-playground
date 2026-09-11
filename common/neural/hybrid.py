"""A network where a few named cells use their own measured biophysics.

    hybrid = HybridRuntime.build(runtime, cells=("RMDDL", "RMDDR", "RMDVL", "RMDVR"))
    hybrid.step()

Almost every cell in this model is a graded leaky integrator, because for almost
every cell that is all anybody has measured. A handful are different: RMD and
AWC^on were characterised channel by channel by Nicoletti et al. 2019, and six
more cells in 2024. This runs those cells on their own published dynamics and
leaves the rest alone.

Why bother
----------

``docs/model_assumptions.md`` §5C.9 established that the graded model has no limit
cycle **anywhere** -- every input drives it to a fixed point, and sustained
variation after the transient measured 4e-5 mV, which is floating-point residue.
A network of relaxation elements cannot produce a gait no matter how it is wired,
so the connectome never got a fair test.

RMD is the cell to start with. It innervates body wall muscle directly (§5C.8),
and it has **plateau potentials** -- it latches on, and stays on until something
switches it off. That bistability is precisely the nonlinearity the graded model
lacks, and two bistable cells wired to inhibit each other is the classic way a
nervous system makes a rhythm.

RMD alone is not a rhythm generator: held at any current from 2 to 150 pA it
latches and sits there (measured; a light switch, not a metronome). The rhythm, if
there is one, has to come from the coupling -- which is exactly the question this
module exists to ask, and which cannot be asked while RMD is a relaxation element.

How the two models are joined
-----------------------------

The conductance model owns the membrane potential of its cells. The network owns
everything else about them: they receive gap-junction and chemical current from
the graded cells like anyone else, and they drive their targets through the same
synaptic activation variable ``s``, which follows their voltage.

Each network step:

1. The synaptic, gap-junction and injected current arriving at each intrinsic cell
   is read from the network, with the rest of the network held fixed.
2. The conductance model is substepped with that as its external current -- it
   needs ``dt <= 0.05 ms`` where the network runs at 1 ms, so it takes many small
   steps inside one large one.
3. The graded runtime steps normally, then the intrinsic cells' voltages are
   overwritten with the conductance model's answer.

That is an explicit operator split: during the substeps the intrinsic cell sees a
frozen network, and during the network step the network sees a frozen intrinsic
cell. It is first-order accurate in the coupling, which is the same order as the
gap-junction coupling in the exponential integrator already, and it avoids
integrating the whole 397-cell network twenty times per millisecond.

Two things that must not be double-counted
------------------------------------------

**Leak.** Nicoletti's RMD has ``g_leak = 0.4 nS`` as one of fourteen ionic
currents. The graded model applies its own ``10 pS`` to every cell. Both would
otherwise act on the same membrane -- forty times too much leak, silently shifting
the resting potential. The graded leak is removed for these cells by construction:
their voltage comes from the conductance model, which never sees it.

**Capacitance.** 1.2 pF measured against 1.0 pF assumed. Same resolution: the
conductance model integrates with its own value.

What is assumed here
--------------------

The activation threshold. Every graded cell's output sigmoid is centred on that
cell's own resting potential, so at rest it sits at its midpoint. An intrinsic
cell rests at -69.5 mV rather than the graded model's operating point, so its
sigmoid is re-centred there to keep the same convention. Without it RMD would sit
at the bottom of its output curve and be silent. This is a modelling choice, and
it decides how strongly a plateau is transmitted downstream.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from common.neural.conductance import DEFAULT_DT_MS, ConductanceModel
from common.neural.runtime import NeuralRuntime

#: Cells with a published conductance-based model, by model id. Only these may be
#: run intrinsically; see the module docstring on why the list is not longer.
INTRINSIC_MODELS: dict[str, tuple[str, ...]] = {
    "rmd_nicoletti2019": ("RMDDL", "RMDDR", "RMDVL", "RMDVR", "RMDL", "RMDR"),
}


@dataclass
class HybridRuntime:
    """A :class:`NeuralRuntime` in which some cells run a conductance model."""

    runtime: NeuralRuntime
    model: ConductanceModel
    cells: tuple[str, ...]
    indices: np.ndarray
    state: np.ndarray = field(init=False)
    substep_ms: float = DEFAULT_DT_MS

    def __post_init__(self) -> None:
        self.state = self.model.resting_state(len(self.cells))
        # Hand the network the intrinsic cells' true resting potential before
        # anything runs, so the first step does not begin with a spurious jump.
        self.runtime.state[0, self.indices] = self.state[0]

    # -- construction ------------------------------------------------------

    @classmethod
    def build(
        cls,
        runtime: NeuralRuntime,
        *,
        model_id: str = "rmd_nicoletti2019",
        cells: tuple[str, ...] | None = None,
        substep_ms: float = DEFAULT_DT_MS,
        recentre_threshold: bool = True,
    ) -> HybridRuntime:
        """Promote ``cells`` to the conductance model, leaving every other cell alone.

        Cells absent from the network are skipped rather than raising, so a body
        model that excludes the head still works.
        """
        wanted = cells if cells is not None else INTRINSIC_MODELS[model_id]
        present = tuple(c for c in wanted if c in runtime.network.cell_ids)
        if not present:
            raise ValueError(
                f"none of {wanted} are in this network; nothing to promote. "
                "Build the runtime with those cells, or pass a different set."
            )
        indices = np.array([runtime.network.index(c) for c in present], dtype=np.int64)
        model = ConductanceModel.load(model_id)

        if recentre_threshold:
            # See the module docstring: keep the convention that a cell sits at the
            # midpoint of its output sigmoid at rest.
            rest = float(model.resting_state()[0, 0])
            thresholds = np.array(
                runtime.model.v_threshold_mv,  # type: ignore[attr-defined]
                dtype=np.float64,
                copy=True,
            )
            thresholds[indices] = rest
            runtime.model = replace(  # type: ignore[type-var]
                runtime.model,
                v_threshold_mv=thresholds,
            )

        return cls(
            runtime=runtime,
            model=model,
            cells=present,
            indices=indices,
            substep_ms=substep_ms,
        )

    # -- the loop ----------------------------------------------------------

    def network_drive(self) -> tuple[np.ndarray, np.ndarray]:
        """What the network presents to the intrinsic cells: ``(current, conductance)``.

        Split deliberately rather than returned as a single current. Both the gap
        and chemical terms have the form ``constant - g * V``:

            i_gap = Ggap @ v            - v * rowsum(Ggap)
            i_syn = (Gsyn o E) @ s      - v * (Gsyn @ s)

        With the rest of the network held fixed across a step, the first part of
        each is a constant and the second is a conductance the intrinsic cell sees
        on its own membrane. Handing over only the total current at the start of the
        step discards the conductance, and with it the restoring force: a gap
        junction pulls a cell back toward its neighbours, and a frozen current
        cannot. RMD carries up to 3.8 nS of gap conductance on 1.2 pF, so the frozen
        version diverged within nine milliseconds.

        The graded leak is excluded, because these cells have their own (see the
        module docstring).
        """
        net = self.runtime.network
        v, s = self.runtime.state[0], self.runtime.state[1]
        constant = (net.g_gap @ v + net.g_syn_e_rev @ s + self.runtime.i_ext_pa)[self.indices]
        conductance = (net.gap_row_sum + net.g_syn @ s)[self.indices]
        return np.asarray(constant), np.asarray(conductance)

    def step(self, dt_ms: float | None = None) -> None:
        """Advance the whole network one step."""
        dt = self.runtime.dt_ms if dt_ms is None else dt_ms
        current, conductance = self.network_drive()

        # The intrinsic cells, substepped against a frozen network.
        steps = max(1, int(round(dt / self.substep_ms)))
        inner = dt / steps
        for _ in range(steps):
            self.state = self.model.step(self.state, inner, current, conductance)

        self.runtime.step(dt)
        # The conductance model owns these voltages; discard what the graded model
        # computed for them. Their synaptic activation `s` is left to the graded
        # runtime, which drives it from this voltage, so they signal downstream
        # through the same path as every other cell.
        self.runtime.state[0, self.indices] = self.state[0]

    def run(self, duration_ms: float, dt_ms: float | None = None) -> None:
        dt = self.runtime.dt_ms if dt_ms is None else dt_ms
        for _ in range(int(round(duration_ms / dt))):
            self.step(dt)

    # -- inspection --------------------------------------------------------

    def voltage(self, cell_id: str) -> float:
        return self.runtime.voltage(cell_id)

    def voltages(self) -> dict[str, float]:
        return dict(zip(self.cells, self.state[0].tolist(), strict=True))

    def describe(self) -> str:
        return (
            f"{len(self.cells)} cells on {self.model.model_id}: "
            f"{', '.join(self.cells)}\n"
            f"  substep {self.substep_ms} ms inside a {self.runtime.dt_ms} ms network step\n"
            f"  resting: " + ", ".join(f"{c} {v:.1f}mV" for c, v in self.voltages().items())
        )
