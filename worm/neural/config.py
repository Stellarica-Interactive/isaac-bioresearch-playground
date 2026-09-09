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
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files

from common.data.schemas import CellCategory, Connectome
from common.neural.neuron_models import GradedLeakyIntegratorParameters
from common.neural.runtime import Integrator, NeuralRuntime, StabilityReport
from common.neural.synapses import UnknownSignPolicy, WeightScaling
from worm.loader import load

CONFIG_FILE = "parameters.toml"

#: Overlays needed for a meaningful simulation. ``polarity`` is included here and
#: nowhere else by default: the runtime is the one place synaptic sign is actually
#: used, so this is where the opt-in belongs and where it must be reported.
RUNTIME_OVERLAYS = ("classes", "sim", "nt", "polarity")


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
    """
    policy = UnknownSignPolicy(unknown_sign)
    scaling = WeightScaling(weight_scaling)

    if connectome is None:
        connectome, reports = load(dataset_id, annotations=RUNTIME_OVERLAYS)
        polarity = next(
            (r.summary() for r in reports if r.overlay_id == "polarity"), "not applied"
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
