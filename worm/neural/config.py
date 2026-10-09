"""Load the *C. elegans* neural model configuration.

Parameters live in ``parameters.toml`` rather than in code so that every one of
them carries a visible confidence level and citation, and so that a run can say
exactly which values it used. See that file before trusting any number.

    from worm.neural.config import build_runtime
    runtime, report = build_runtime("cook_2019_herm", unknown_sign="exclude")
"""

from __future__ import annotations

import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, replace
from functools import lru_cache
from importlib.resources import files

import numpy as np

from common.data.schemas import CellCategory, Connectome
from common.neural.neuron_models import (
    GradedLeakyIntegrator,
    GradedLeakyIntegratorParameters,
    prepare,
)
from common.neural.runtime import Integrator, NeuralRuntime, StabilityReport
from common.neural.synapses import UnknownSignPolicy, WeightScaling
from worm.loader import load

CONFIG_FILE = "parameters.toml"

#: Overlays needed for a meaningful simulation. The two sign overlays are included
#: here and nowhere else by default: the runtime is the one place synaptic sign is
#: actually used, so this is where the opt-in belongs and where it gets reported.
#:
#: They are not equivalent. ``polarity`` is *predicted* from gene expression and
#: covers interneuronal connections; ``nmj`` is *measured* physiology and covers the
#: body wall neuromuscular junctions that drive movement. The reports keep them
#: separate so a result can say which it leaned on.
RUNTIME_OVERLAYS = ("classes", "sim", "nt", "polarity", "nmj")


@dataclass(frozen=True)
class ParameterProvenance:
    """One parameter, with where it came from and how much to trust it."""

    name: str
    value: float
    unit: str
    confidence: str
    source: str
    note: str = ""
    improvement: str = ""

    @property
    def is_assumed(self) -> bool:
        return self.confidence == "assumed"


@lru_cache(maxsize=1)
def _raw_config() -> dict[str, object]:
    text = (files("worm.neural") / CONFIG_FILE).read_text(encoding="utf-8")
    return tomllib.loads(text)


def parameter_provenance() -> list[ParameterProvenance]:
    """Every parameter with its confidence and citation, for reporting."""
    params = _raw_config()["parameters"]
    assert isinstance(params, dict)
    return [
        ParameterProvenance(
            name=name,
            value=float(entry["value"]),
            unit=str(entry.get("unit", "")),
            confidence=str(entry["confidence"]),
            source=str(entry.get("source", "")),
            note=str(entry.get("note", "")).strip(),
            improvement=str(entry.get("improvement", "")).strip(),
        )
        for name, entry in sorted(params.items())
    ]


def model_parameters(
    overrides: Mapping[str, float] | None = None,
) -> GradedLeakyIntegratorParameters:
    """Parameters from the config file, with optional overrides for a sweep."""
    values = {p.name: p.value for p in parameter_provenance()}
    if overrides:
        unknown = set(overrides) - set(values)
        if unknown:
            raise KeyError(f"not parameters of this model: {sorted(unknown)}")
        values.update(overrides)
    return GradedLeakyIntegratorParameters(**values)


def body_wall_muscle_gap_scale() -> float:
    """The committed factor on gap junctions between body-wall muscle cells.

    From the measured coupling ratio, not assumed; see the ``[body_wall_muscle]``
    section of ``parameters.toml``.
    """
    section = _raw_config().get("body_wall_muscle", {})
    assert isinstance(section, dict)
    return float(section.get("gap_scale", 1.0))


def integration_defaults() -> tuple[float, Integrator]:
    section = _raw_config()["integration"]
    assert isinstance(section, dict)
    return float(section["dt_ms"]), Integrator(str(section["integrator"]))


def assumed_parameters() -> list[str]:
    """Names of every parameter with no direct biological backing.

    This is most of them, which is the point of surfacing the list.
    """
    return [p.name for p in parameter_provenance() if p.is_assumed]


def provenance_report() -> str:
    lines = [
        "Neural model parameters",
        "",
        f"{'parameter':16s} {'value':>9s} {'unit':6s} confidence",
        "-" * 62,
    ]
    for p in parameter_provenance():
        lines.append(f"{p.name:16s} {p.value:9.4g} {p.unit:6s} {p.confidence}")
    assumed = assumed_parameters()
    lines += [
        "",
        f"{len(assumed)} of {len(parameter_provenance())} parameters are ASSUMED: "
        + ", ".join(assumed),
        "",
        "In this model's lineage the connectivity is the only measured quantity.",
        "See worm/neural/parameters.toml for the provenance of each value.",
    ]
    return "\n".join(lines)


@dataclass(frozen=True)
class RuntimeReport:
    """What a built runtime is actually standing on."""

    dataset_id: str
    n_cells: int
    unknown_sign_policy: UnknownSignPolicy
    weight_scaling: WeightScaling
    sign_summary: str
    polarity_coverage: str
    assumed_parameters: tuple[str, ...]
    stability: StabilityReport

    def summary(self) -> str:
        return (
            f"{self.dataset_id}: {self.n_cells} cells\n"
            f"  sign: {self.sign_summary}\n"
            f"  polarity overlay: {self.polarity_coverage}\n"
            f"  unknown-sign policy: {self.unknown_sign_policy}, "
            f"weight scaling: {self.weight_scaling}\n"
            f"  {self.stability.summary()}\n"
            f"  {len(self.assumed_parameters)} assumed parameters: "
            f"{', '.join(self.assumed_parameters)}"
        )


def build_runtime(
    dataset_id: str = "cook_2019_herm",
    *,
    unknown_sign: UnknownSignPolicy | str,
    cells: str | tuple[str, ...] = "neurons",
    weight_scaling: WeightScaling | str = WeightScaling.LINEAR,
    dt_ms: float | None = None,
    integrator: Integrator | str | None = None,
    connectome: Connectome | None = None,
    parameter_overrides: Mapping[str, float] | None = None,
    muscle_gap_scale: float | None = None,
) -> tuple[NeuralRuntime, RuntimeReport]:
    """Build a runtime from a named dataset, reporting what it rests on.

    ``unknown_sign`` is required and has no default; see
    :class:`~common.neural.synapses.UnknownSignPolicy`.

    ``cells`` may be ``"neurons"`` (the default -- see the ``[network]`` note in
    ``parameters.toml`` for why muscles are excluded for now), ``"all"``, or an
    explicit tuple of cell ids.

    ``parameter_overrides`` is a mapping rather than keyword arguments so that
    biophysical parameters cannot be confused with the build options above -- which
    matters when sweeping, since almost every parameter is an assumption.

    When body-wall muscle cells are in the network, the gap junctions between
    them are scaled by ``muscle_gap_scale`` -- by default the committed,
    measurement-derived :func:`body_wall_muscle_gap_scale`; pass 1.0 for the
    uniform coupling every result before model_assumptions 5AO.5 used.
    """
    policy = UnknownSignPolicy(unknown_sign)
    scaling = WeightScaling(weight_scaling)

    if connectome is None:
        connectome, reports = load(dataset_id, annotations=RUNTIME_OVERLAYS)
        # Both sign overlays are reported, and separately: 'polarity' is predicted
        # from gene expression, 'nmj' is measured physiology. Collapsing them into
        # one number would hide which kind of evidence a result leaned on.
        by_id = {r.overlay_id: r.summary() for r in reports}
        polarity = "; ".join(
            by_id.get(name, f"{name}: not applied") for name in ("polarity", "nmj")
        )
    else:
        polarity = "supplied by caller"

    if cells == "neurons":
        selected: tuple[str, ...] | None = tuple(
            c.id for c in connectome.cells if c.category is CellCategory.NEURON
        )
    elif cells == "all":
        selected = None
    else:
        selected = tuple(cells)

    default_dt, default_integrator = integration_defaults()
    runtime = NeuralRuntime.build(
        connectome,
        unknown_sign=policy,
        params=model_parameters(parameter_overrides),
        weight_scaling=scaling,
        cells=selected,
        dt_ms=default_dt if dt_ms is None else dt_ms,
        integrator=default_integrator if integrator is None else Integrator(integrator),
    )
    scale = body_wall_muscle_gap_scale() if muscle_gap_scale is None else float(muscle_gap_scale)
    in_network = set(runtime.network.cell_ids)
    has_muscle = any(
        c.category is CellCategory.MUSCLE for c in connectome.cells if c.id in in_network
    )
    if has_muscle and scale != 1.0:
        rescale_gap_junctions(runtime, connectome, {"muscle": scale})
    report = RuntimeReport(
        dataset_id=dataset_id,
        n_cells=runtime.network.n,
        unknown_sign_policy=policy,
        weight_scaling=scaling,
        sign_summary=runtime.network.census.summary(),
        polarity_coverage=polarity,
        assumed_parameters=tuple(assumed_parameters()),
        stability=runtime.stability_report(),
    )
    return runtime, report


#: The kinds of gap junction :func:`rescale_gap_junctions` can scale separately,
#: named by the two cells a junction joins.
GAP_JUNCTION_CLASSES = ("neuron", "muscle", "neuron-muscle")


def rescale_gap_junctions(
    runtime: NeuralRuntime, connectome: Connectome, factors: Mapping[str, float]
) -> None:
    """Scale gap junctions by the kinds of cell they join, then rebuild as a build would.

    ``g_gap_ps`` is one number for every gap junction in the network: between
    neurons, between body-wall muscle cells, and between the two. Nothing measured
    says they share a conductance, and model_assumptions 5AN.20 needed all of them
    twenty times weaker for the loop to crawl forward. This lets the three be set
    apart, so a result can say which of them it actually needed.

    ``factors`` maps a name in :data:`GAP_JUNCTION_CLASSES` to a multiplier;
    classes not named are left alone. Every cell's sigmoid midpoint is then
    re-solved against the new network and the state reset to its resting state --
    exactly what building the runtime with a different ``g_gap_ps`` does -- so
    scaling all three classes by one factor is the same model as changing
    ``g_gap_ps`` by it, which a test holds.
    """
    unknown = set(factors) - set(GAP_JUNCTION_CLASSES)
    if unknown:
        raise ValueError(
            f"unknown gap-junction classes {sorted(unknown)}; use {GAP_JUNCTION_CLASSES}"
        )
    category = {c.id: c.category for c in connectome.cells}
    muscle = np.array([category.get(c) is CellCategory.MUSCLE for c in runtime.network.cell_ids])
    both_muscle = muscle[:, None] & muscle[None, :]
    neither = ~muscle[:, None] & ~muscle[None, :]
    masks = {"muscle": both_muscle, "neuron": neither, "neuron-muscle": ~both_muscle & ~neither}
    g = runtime.network.g_gap.copy()
    for name, factor in factors.items():
        g[masks[name]] *= float(factor)
    network = replace(runtime.network, g_gap=g)
    model = prepare(GradedLeakyIntegrator(params=runtime.model.params), network)
    runtime.network = network
    runtime.model = model
    runtime.state = model.initial_state(network)
