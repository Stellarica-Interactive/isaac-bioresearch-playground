"""Load the *C. elegans* neural model configuration.

Parameters live in ``parameters.toml`` rather than in code so that every one of
them carries a visible confidence level and citation, and so that a run can say
exactly which values it used. See that file before trusting any number.

    from worm.neural.config import build_runtime
    runtime, report = build_runtime("cook_2019_herm", unknown_sign="exclude")
"""

from __future__ import annotations

import tomllib
from collections.abc import Mapping, Sequence
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


@dataclass(frozen=True)
class ExtrasynapticLink:
    """A coupling between two cells that is not a synapse, so no connectome has it.

    Neuromodulators and transmitter released outside synapses reach receptors on
    cells the releasing cell never contacts. Shen et al. 2016 measured two such
    links in the head (model_assumptions 5AO.11), and nothing in an electron
    micrograph shows either. ``sign`` is "excitatory" or "inhibitory";
    ``g_ns`` is the conductance, which for these links is never measured.
    """

    pre: str
    post: str
    g_ns: float
    sign: str


def add_extrasynaptic(runtime: NeuralRuntime, links: Sequence[ExtrasynapticLink]) -> list[str]:
    """Add couplings the connectome cannot show, then rebuild as a build would.

    Each link is a graded chemical connection from ``pre`` onto ``post``, with the
    network's own synaptic dynamics and the reversal potential of its sign. Where
    the pair is already connected the two combine exactly: both are gated by the
    same presynaptic activation, so the conductances add and the reversal
    potential is their conductance-weighted mean. Links naming a cell that is not
    in the network -- removed by a lesion, say -- are skipped. Every cell's
    sigmoid midpoint is then re-solved and the state reset, as
    :func:`rescale_gap_junctions` does, so the network rests where a network built
    with these connections would. Returns the links added, as "PRE->POST".
    """
    net = runtime.network
    params = runtime.model.params
    reversal = {"excitatory": params.e_exc_mv, "inhibitory": params.e_inh_mv}
    g_syn = net.g_syn.copy()
    e_rev = net.e_rev.copy()
    added = []
    for link in links:
        if link.sign not in reversal:
            raise ValueError(f"sign must be excitatory or inhibitory, not {link.sign!r}")
        if link.g_ns < 0.0:
            raise ValueError(f"negative conductance for {link.pre}->{link.post}")
        if link.pre not in net.cell_ids or link.post not in net.cell_ids:
            continue
        i, j = net.index(link.post), net.index(link.pre)
        g_old = g_syn[i, j]
        g_new = g_old + link.g_ns
        if g_new > 0.0:
            e_rev[i, j] = (g_old * e_rev[i, j] + link.g_ns * reversal[link.sign]) / g_new
        g_syn[i, j] = g_new
        added.append(f"{link.pre}->{link.post}")
    network = replace(net, g_syn=g_syn, e_rev=e_rev)
    model = prepare(GradedLeakyIntegrator(params=params), network)
    runtime.network = network
    runtime.model = model
    runtime.state = model.initial_state(network)
    return added


#: The pairs Shen et al. 2016 recorded active together: SMDD and RMED during dorsal
#: head bends, SMDV and RMEV during ventral.
_SMD_RME_PAIRS = (
    ("SMDDL", "RMED"),
    ("SMDDR", "RMED"),
    ("SMDVL", "RMEV"),
    ("SMDVR", "RMEV"),
)


def smd_rme_loop(
    excite_ns: float, inhibit_ns: float, *, pairing: str = "same"
) -> list[ExtrasynapticLink]:
    """The head's measured gain control, as links (Shen et al. 2016; 5AO.11).

    SMD excites RME extrasynaptically, through the muscarinic receptor GAR-2 on
    RME, and RME inhibits SMD extrasynaptically, through the GABA-B receptor
    GBB-1/GBB-2 on SMD. Measured: both signs, that RME's activity follows SMD's
    and the head's bending (SMDD with RMED, SMDV with RMEV), and the effect --
    without RME the head bends further, with RME driven it bends less. Not
    measured: either conductance, and which SMD each RME inhibits. ``pairing``
    "same" returns each RME's inhibition to the SMD cells that drive it, a gain
    control per side; "all" lets each RME inhibit all four.
    """
    if pairing not in ("same", "all"):
        raise ValueError(f"pairing must be 'same' or 'all', not {pairing!r}")
    links = [ExtrasynapticLink(smd, rme, excite_ns, "excitatory") for smd, rme in _SMD_RME_PAIRS]
    if pairing == "same":
        targets = [(rme, smd) for smd, rme in _SMD_RME_PAIRS]
    else:
        targets = [(rme, smd) for rme in ("RMED", "RMEV") for smd, _ in _SMD_RME_PAIRS]
    links += [ExtrasynapticLink(rme, smd, inhibit_ns, "inhibitory") for rme, smd in targets]
    return links
