"""Annotation overlays for *C. elegans*.

Each overlay attaches one kind of published annotation to a connectome that was
parsed from electron-microscopy data alone. See :mod:`common.data.overlay` for
why these are kept separate from the anatomy.

Three overlays exist today:

``classes``
    Anatomical neuron class (``AVAL`` -> ``AVA``), from Wang et al. 2024.

``sim``
    Coarse sensory / interneuron / motor role, from WormAtlas cell listings and
    the Cook et al. 2019 groupings.

``nt``
    Neurotransmitter usage, from the Wang et al. 2024 reporter-allele atlas.

A polarity overlay (which chemical synapses are excitatory or inhibitory) is
**not** included yet. That value is not measured anywhere; it can only be
predicted from receptor gene expression, and letting a prediction reach a
dynamics model without a deliberate opt-in is the most likely way this project
would end up fooling itself. See ``docs/model_assumptions.md``.
"""

from __future__ import annotations

import csv
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, replace
from functools import lru_cache
from importlib.resources import files

from common.data.overlay import AnnotationOverlay, OverlayReport
from common.data.schemas import (
    Cell,
    CellCategory,
    Confidence,
    Connection,
    Connectome,
    Neurotransmitter,
    Provenance,
    Sign,
    SIMRole,
    SynapseType,
)
from worm.importers.sources import source_file


@lru_cache(maxsize=8)
def _rows(name: str) -> tuple[Mapping[str, str], ...]:
    text = (files("worm.data") / "annotations" / name).read_text(encoding="utf-8")
    return tuple(csv.DictReader(text.splitlines()))


def _provenance(source_id: str, note: str) -> Provenance:
    spec = source_file(source_id)
    return Provenance(
        source_id=source_id,
        kind="annotation",
        citation=spec.citation,
        url=spec.url,
        license=spec.license,
        doi=spec.doi,
        file_name=spec.filename or None,
        sha256=spec.sha256 or None,
        notes=note,
    )


def _report(
    overlay_id: str, annotated: set[str], connectome: Connectome, table_ids: set[str]
) -> OverlayReport:
    eligible = {c.id for c in connectome.cells if c.category is CellCategory.NEURON}
    return OverlayReport(
        overlay_id=overlay_id,
        matched=len(annotated),
        eligible=len(eligible),
        unmatched=tuple(sorted(eligible - annotated)),
        unknown_in_source=tuple(sorted(table_ids - {c.id for c in connectome.cells})),
    )


def _rebuild(connectome: Connectome, updated: Mapping[str, Cell]) -> Connectome:
    return replace(
        connectome, cells=tuple(updated.get(c.id, c) for c in connectome.cells)
    )


# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NeuronClassOverlay:
    """Anatomical neuron class: the group a bilateral or radial set belongs to.

    ``AVAL`` and ``AVAR`` are both class ``AVA``; ``RMDDL``, ``RMDDR``, ``RMDVL``,
    ``RMDVR``, ``RMDL`` and ``RMDR`` are all class ``RMD``. This grouping is
    anatomical, published, and *not* derivable by stripping letters from the name,
    which is why it is read from a table rather than computed.
    """

    overlay_id: str = "classes"
    source_id: str = "wang_2024_nt_atlas"

    def provenance(self) -> Provenance:
        return _provenance(self.source_id, "Neuron class column of Supplementary File 2.")

    def apply(self, c: Connectome) -> tuple[Connectome, OverlayReport]:
        table = {r["cell_id"]: r["class_name"] for r in _rows("neuron_classes.csv")}
        updated: dict[str, Cell] = {}
        for cell in c.cells:
            klass = table.get(cell.id)
            if klass is None:
                continue
            updated[cell.id] = replace(
                cell,
                class_name=klass,
                field_sources={**cell.field_sources, "class_name": self.source_id},
            )
        return _rebuild(c, updated), _report(self.overlay_id, set(updated), c, set(table))


@dataclass(frozen=True)
class SIMRoleOverlay:
    """Coarse sensory / interneuron / motor classification.

    Not a measurement made by any connectome. It is a curated summary of decades
    of ablation, imaging and anatomy work. A cell may carry several roles, and
    seven cells (CANL/R, MI, MCL/R, NSML/R) carry
    :attr:`~common.data.schemas.SIMRole.UNKNOWN` because the published label for
    them does not resolve to a single role and we decline to invent one.
    """

    overlay_id: str = "sim"
    source_id: str = "wormatlas_cook_2019_via_cect"

    def provenance(self) -> Provenance:
        return _provenance(self.source_id, "Coarse functional role per cell.")

    def apply(self, c: Connectome) -> tuple[Connectome, OverlayReport]:
        table: dict[str, list[SIMRole]] = {}
        for r in _rows("sim_roles.csv"):
            table.setdefault(r["cell_id"], []).append(SIMRole(r["role"]))
        updated: dict[str, Cell] = {}
        for cell in c.cells:
            roles = table.get(cell.id)
            if not roles:
                continue
            updated[cell.id] = replace(
                cell,
                roles=tuple(roles),
                field_sources={**cell.field_sources, "roles": self.source_id},
            )
        return _rebuild(c, updated), _report(self.overlay_id, set(updated), c, set(table))


@dataclass(frozen=True)
class NeurotransmitterOverlay:
    """Neurotransmitter usage from CRISPR reporter-allele expression.

    The source distinguishes three things that a naive import would flatten, and
    which are preserved in ``worm/data/annotations/neurotransmitters.csv``:

    * clean reporter expression;
    * *dim and variable* expression — real but weaker evidence, marked ``*`` in
      the source;
    * **uptake rather than synthesis** — the cell does not make the transmitter
      but takes it up from neighbours (AVFL/R for GABA; CANL/R and ASIL/R for
      betaine). Such a cell can still release the transmitter, but the genetics
      and the manipulability differ.

    Sixteen neurons express no known transmitter pathway gene at all. That is a
    published *result*, not a missing measurement, and it is recorded as
    :attr:`~common.data.schemas.Neurotransmitter.UNKNOWN` with the evidence string
    ``orphan_no_pathway_gene_detected`` rather than left blank.
    """

    overlay_id: str = "nt"
    source_id: str = "wang_2024_nt_atlas"

    def provenance(self) -> Provenance:
        return _provenance(
            self.source_id,
            "Reporter-allele neurotransmitter assignments, Supplementary File 2 "
            "(hermaphrodite).",
        )

    def apply(self, c: Connectome) -> tuple[Connectome, OverlayReport]:
        table: dict[str, list[Neurotransmitter]] = {}
        for r in _rows("neurotransmitters.csv"):
            table.setdefault(r["cell_id"], []).append(Neurotransmitter(r["neurotransmitter"]))
        updated: dict[str, Cell] = {}
        for cell in c.cells:
            nts = table.get(cell.id)
            if not nts:
                continue
            updated[cell.id] = replace(
                cell,
                neurotransmitters=tuple(nts),
                field_sources={**cell.field_sources, "neurotransmitters": self.source_id},
            )
        return _rebuild(c, updated), _report(self.overlay_id, set(updated), c, set(table))


@dataclass(frozen=True)
class PolarityOverlay:
    """Predicted excitatory/inhibitory sign for chemical synapses. **Opt-in.**

    Not in :data:`DEFAULT_OVERLAYS`, and that is deliberate. No experiment measures
    synaptic sign across a connectome; these values combine the presynaptic
    transmitter with postsynaptic receptor gene expression, so every one is
    :attr:`~common.data.schemas.Confidence.PREDICTED`. A neural simulation is
    extremely sensitive to whether a connection is + or -, so a prediction arriving
    unnoticed is the fastest route to confident nonsense.

    Three limits worth knowing before switching it on:

    * **Roughly half of connections get a definite sign** (1327 excitatory, 425
      inhibitory of 3638). 471 are ``mixed`` -- the target expresses both excitatory
      and inhibitory receptors for that transmitter, so the net effect is genuinely
      undetermined rather than merely unmeasured -- and 1415 have no receptor match
      at all.
    * **No neuromuscular junctions.** The source covers interneuronal connections
      only, so it supplies no sign for the synapses that actually drive muscle. See
      ``docs/model_assumptions.md``.
    * Predictions are made against the WormWiring/Cook reconstruction, so coverage of
      a different dataset's edges will be partial. The overlay reports what it matched.
    """

    overlay_id: str = "polarity"
    source_id: str = "fenyves_2020_polarity"

    def provenance(self) -> Provenance:
        return _provenance(
            self.source_id,
            "PREDICTED synaptic polarity from neurotransmitter and receptor gene "
            "expression (S1 Data, NT+R method). Never measured.",
        )

    def apply(self, c: Connectome) -> tuple[Connectome, OverlayReport]:
        table: dict[tuple[str, str], Sign] = {}
        for r in _rows("polarity_fenyves2020.csv"):
            table[(r["pre"], r["post"])] = Sign(r["sign"])

        updated: list[Connection] = []
        matched = 0
        eligible = 0
        for e in c.connections:
            if e.synapse_type is not SynapseType.CHEMICAL:
                updated.append(e)
                continue
            eligible += 1
            sign = table.get((e.pre, e.post))
            if sign is None:
                updated.append(e)
                continue
            matched += 1
            updated.append(
                replace(
                    e,
                    sign=sign,
                    sign_confidence=Confidence.PREDICTED,
                    field_sources={**e.field_sources, "sign": self.source_id},
                )
            )

        chemical_pairs = {(e.pre, e.post) for e in c.chemical()}
        report = OverlayReport(
            overlay_id=self.overlay_id,
            matched=matched,
            eligible=eligible,
            unmatched=tuple(
                sorted(f"{e.pre}->{e.post}" for e in c.chemical() if (e.pre, e.post) not in table)
            ),
            unknown_in_source=tuple(
                sorted(f"{a}->{b}" for a, b in table if (a, b) not in chemical_pairs)
            ),
        )
        return replace(c, connections=tuple(updated)), report


#: Sign of a fast chemical synapse onto **body wall muscle**, by presynaptic
#: transmitter, with the paper that establishes each. Deliberately short: it covers
#: only the two transmitters whose effect on body muscle has been measured directly.
#:
#: Glutamate, dopamine and unknown-transmitter synapses onto muscle are left
#: unsigned. Dopamine in particular acts through G-protein-coupled receptors on a
#: slow, modulatory timescale, so giving it a reversal potential would not be a
#: cautious guess -- it would be the wrong kind of model.
NMJ_SIGN_BY_TRANSMITTER: dict[Neurotransmitter, tuple[Sign, str]] = {
    Neurotransmitter.ACETYLCHOLINE: (Sign.EXCITATORY, "richmond_1999_nmj_receptors"),
    Neurotransmitter.GABA: (Sign.INHIBITORY, "mcintire_1993_gaba_inhibitory"),
}


@dataclass(frozen=True)
class NeuromuscularPolarityOverlay:
    """Sign for the synapses that actually drive muscle. **Opt-in.**

    Why this exists separately from :class:`PolarityOverlay`: the Fenyves
    predictions cover interneuronal connections only, so every one of the ~1000
    neuron-to-muscle synapses is unsigned. Those are precisely the synapses a
    locomotion model depends on, so without this there is no motor output at all.

    Why it is better evidence than the interneuron predictions, rather than more
    of the same guessing: acetylcholine and GABA at the *C. elegans* body wall
    neuromuscular junction were established by direct patch-clamp recording and
    mutant analysis, not by inference from gene expression. Muscle expresses two
    nicotinic acetylcholine receptors and one GABA receptor, and ``unc-49`` is
    required postsynaptically for GABA's inhibitory effect on body muscle. So
    these carry :attr:`~common.data.schemas.Confidence.PUBLISHED_ANNOTATION`
    rather than ``PREDICTED``.

    Scope, deliberately narrow:

    * **Body wall muscle only.** Pharyngeal muscle has different pharmacology --
      glutamate is inhibitory there, through a glutamate-gated chloride channel --
      so applying a body-wall rule to it would be wrong. Vulval, uterine, anal and
      intestinal muscle are likewise left alone.
    * **Acetylcholine and GABA only**, covering 93% of body wall neuromuscular
      junctions in Cook 2019. The remaining 7% (glutamate from IL1 and RIM,
      dopamine from ADE and CEP, and a handful with no known transmitter) stay
      unsigned and are reported.

    Requires the ``nt`` overlay to have run first, since it reads the presynaptic
    transmitter. :func:`get_overlays` orders them correctly.

    A satisfying check that falls out of this: the cholinergic cells innervating
    body wall muscle turn out to be exactly the AS, DA, DB, VA, VB and VC classes,
    and the GABAergic ones exactly DD and VD -- which is the textbook division into
    excitatory and inhibitory motor neurons, arrived at from two independent
    datasets rather than assumed.
    """

    overlay_id: str = "nmj"
    source_id: str = "richmond_1999_nmj_receptors"

    def provenance(self) -> Provenance:
        return _provenance(
            self.source_id,
            "Sign of body wall neuromuscular junctions from presynaptic transmitter "
            "identity. Measured physiology, not a gene-expression prediction.",
        )

    def apply(self, c: Connectome) -> tuple[Connectome, OverlayReport]:
        from worm.importers.naming import body_wall_muscle_ids

        muscles = body_wall_muscle_ids()
        transmitters = {cell.id: cell.neurotransmitters for cell in c.cells}

        updated: list[Connection] = []
        matched = 0
        eligible = 0
        unmatched: list[str] = []

        for e in c.connections:
            if e.synapse_type is not SynapseType.CHEMICAL or e.post not in muscles:
                updated.append(e)
                continue
            eligible += 1

            pre_nt = transmitters.get(e.pre, ())
            assignment = next(
                (NMJ_SIGN_BY_TRANSMITTER[n] for n in pre_nt if n in NMJ_SIGN_BY_TRANSMITTER),
                None,
            )
            if assignment is None:
                unmatched.append(f"{e.pre}->{e.post}")
                updated.append(e)
                continue

            sign, source = assignment
            matched += 1
            updated.append(
                replace(
                    e,
                    sign=sign,
                    sign_confidence=Confidence.PUBLISHED_ANNOTATION,
                    field_sources={**e.field_sources, "sign": source},
                )
            )

        report = OverlayReport(
            overlay_id=self.overlay_id,
            matched=matched,
            eligible=eligible,
            unmatched=tuple(sorted(unmatched)),
            unknown_in_source=(),
        )
        return replace(c, connections=tuple(updated)), report


def polarity_evidence() -> Iterator[Mapping[str, str]]:
    """Full polarity table including the basis string and source row."""
    yield from _rows("polarity_fenyves2020.csv")


def neurotransmitter_evidence() -> Iterator[Mapping[str, str]]:
    """Full neurotransmitter table including evidence strings and source rows."""
    yield from _rows("neurotransmitters.csv")


OVERLAYS: dict[str, AnnotationOverlay] = {
    "classes": NeuronClassOverlay(),
    "sim": SIMRoleOverlay(),
    "nt": NeurotransmitterOverlay(),
    # Both opt-in. Never add either to DEFAULT_OVERLAYS: a simulation is extremely
    # sensitive to synaptic sign, so it should never arrive without being asked for.
    "polarity": PolarityOverlay(),
    "nmj": NeuromuscularPolarityOverlay(),
}

#: Overlays that must run after others. ``nmj`` reads the presynaptic transmitter,
#: so the ``nt`` overlay has to have populated it first.
_OVERLAY_ORDER = ("classes", "sim", "nt", "polarity", "nmj")

DEFAULT_OVERLAYS = ("classes", "sim", "nt")


def get_overlays(names: tuple[str, ...] = DEFAULT_OVERLAYS) -> list[AnnotationOverlay]:
    """Resolve overlay names, applying dependency order.

    Order matters: ``nmj`` derives a synapse's sign from its presynaptic
    transmitter, so ``nt`` must have run first. Sorting here rather than trusting
    the caller means a plausible-looking argument order cannot silently produce an
    overlay that matches nothing.
    """
    unknown = [n for n in names if n not in OVERLAYS]
    if unknown:
        raise KeyError(f"unknown overlay(s) {unknown}; known: {sorted(OVERLAYS)}")
    if "nmj" in names and "nt" not in names:
        raise ValueError(
            "the 'nmj' overlay reads each synapse's presynaptic transmitter, so it "
            "needs the 'nt' overlay as well; requesting it alone would silently "
            "annotate nothing"
        )
    return [OVERLAYS[n] for n in sorted(names, key=_OVERLAY_ORDER.index)]
