"""Generate ``docs/neurons.md`` — a reference for every cell in a dataset.

    python tools/generate_neuron_reference.py --dataset witvliet_2021_7

The table is **generated**, never hand-written, so it cannot drift away from the
data. Every column comes from a cited source:

* name expansion, lineage, description  -> WormAtlas cell listings
* functional role                        -> WormAtlas + Cook et al. 2019 groupings
* neurotransmitter                       -> Wang et al. 2024 reporter-allele atlas
* neuron class                           -> Wang et al. 2024
* partner counts                         -> the connectome itself (measured)

The hand-written parts of ``docs/neurons.md`` are confined to a clearly marked
"curated notes" section defined in :data:`CURATED`, each entry carrying its own
citation. Nothing in this file paraphrases biology from memory.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

from common.data.schemas import Cell, CellCategory, Connectome, SIMRole
from common.neural import graph as gr
from worm.importers.naming import cell_registry
from worm.loader import load

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT = REPO_ROOT / "docs" / "neurons.md"


@dataclass(frozen=True)
class Curated:
    """A hand-written note. Every one carries a citation; none is written from memory."""

    summary: str
    source: str


#: Short notes on the cells a reader is most likely to meet first — the ones that
#: stand out in the figures, and the ones the planned experiments depend on.
#: Keyed by *neuron class*, so ``AVA`` covers both AVAL and AVAR.
CURATED: dict[str, Curated] = {
    "ALM": Curated(
        "Anterior gentle-touch receptor. One of six touch receptor neurons; "
        "ablating them abolishes the response to gentle body touch. Anterior touch "
        "is associated with reversal.",
        "Chalfie M, et al. J Neurosci 5:956-964 (1985).",
    ),
    "AVM": Curated(
        "Anterior ventral gentle-touch receptor, born post-embryonically. Works "
        "with ALM for anterior touch.",
        "Chalfie M, et al. J Neurosci 5:956-964 (1985).",
    ),
    "PLM": Curated(
        "Posterior gentle-touch receptor. Posterior touch is associated with "
        "forward acceleration rather than reversal.",
        "Chalfie M, et al. J Neurosci 5:956-964 (1985).",
    ),
    "AVA": Curated(
        "Premotor command interneuron associated with backward locomotion. "
        "Among the most heavily connected cells in the animal.",
        "Chalfie M, et al. J Neurosci 5:956-964 (1985); Kawano T, et al. Neuron 72:572-586 (2011).",
    ),
    "AVB": Curated(
        "Premotor command interneuron associated with forward locomotion.",
        "Chalfie M, et al. J Neurosci 5:956-964 (1985); Kawano T, et al. Neuron 72:572-586 (2011).",
    ),
    "AVD": Curated(
        "Premotor interneuron in the backward-locomotion pathway, downstream of "
        "the anterior touch receptors.",
        "Chalfie M, et al. J Neurosci 5:956-964 (1985).",
    ),
    "AVE": Curated(
        "Premotor interneuron associated with backward locomotion, acting on the "
        "anterior body.",
        "Kawano T, et al. Neuron 72:572-586 (2011).",
    ),
    "PVC": Curated(
        "Premotor command interneuron associated with forward locomotion, "
        "downstream of the posterior touch receptors.",
        "Chalfie M, et al. J Neurosci 5:956-964 (1985).",
    ),
    "ASH": Curated(
        "Polymodal nociceptor: responds to noxious chemicals, high osmolarity and "
        "nose touch, and drives avoidance. A good example of why a single "
        "sensory/inter/motor label is a simplification.",
        "Kaplan JM, Horvitz HR. PNAS 90:2227-2231 (1993).",
    ),
    "ASE": Curated(
        "Principal salt-sensing chemosensory pair, and functionally left/right "
        "asymmetric: ASEL responds to salt increases and ASER to decreases. The "
        "two members of one class are not interchangeable.",
        "Suzuki H, et al. Nature 454:114-117 (2008).",
    ),
    "AWA": Curated(
        "Olfactory neuron for attractive volatile odours. Notable for the model: "
        "AWA fires genuine calcium-mediated all-or-none action potentials, which "
        "the common 'C. elegans neurons are graded, not spiking' summary omits.",
        "Bargmann CI, et al. Cell 74:515-527 (1993); Liu Q, et al. Cell 175:57-70 (2018).",
    ),
    "AWB": Curated(
        "Olfactory neuron mediating avoidance of repulsive volatile odours.",
        "Troemel ER, et al. Cell 91:161-169 (1997).",
    ),
    "AWC": Curated(
        "Olfactory neuron for attractive volatile odours, and asymmetric between "
        "left and right in the odorants it detects.",
        "Bargmann CI, et al. Cell 74:515-527 (1993); "
        "Wes PD, Bargmann CI. Nature 410:698-701 (2001).",
    ),
    "AIY": Curated(
        "First-layer amphid interneuron; a major integration point downstream of "
        "chemosensory and thermosensory input, implicated in the turn/run decision "
        "underlying chemotaxis.",
        "Gray JM, Hill JJ, Bargmann CI. PNAS 102:3184-3191 (2005).",
    ),
    "AIB": Curated(
        "First-layer amphid interneuron acting largely in opposition to AIY in "
        "the chemotaxis turning decision.",
        "Gray JM, Hill JJ, Bargmann CI. PNAS 102:3184-3191 (2005).",
    ),
    "AIA": Curated(
        "First-layer amphid interneuron receiving convergent chemosensory input.",
        "White JG, et al. Phil Trans R Soc Lond B 314:1-340 (1986).",
    ),
    "AIZ": Curated(
        "First-layer amphid interneuron in the chemotaxis pathway.",
        "Gray JM, Hill JJ, Bargmann CI. PNAS 102:3184-3191 (2005).",
    ),
    "RIM": Curated(
        "Motor/interneuron coupled to the backward-locomotion circuit; its "
        "published label records disagreement between studies about whether it is "
        "an interneuron or a motor neuron.",
        "Kawano T, et al. Neuron 72:572-586 (2011).",
    ),
    "RIA": Curated(
        "Highly connected ring interneuron integrating sensory input with head "
        "motor output; compartmentalized calcium dynamics within one cell.",
        "Hendricks M, et al. Nature 487:99-103 (2012).",
    ),
    "RMG": Curated(
        "Hub of a gap-junction 'hub-and-spoke' circuit: several sensory neurons "
        "are electrically coupled to RMG, which aggregates them.",
        "Macosko EZ, et al. Nature 458:1171-1175 (2009).",
    ),
    "URX": Curated(
        "Oxygen-sensing neuron, electrically coupled into the RMG hub circuit.",
        "Macosko EZ, et al. Nature 458:1171-1175 (2009).",
    ),
    "BAG": Curated(
        "Sensory neuron responding to carbon dioxide and to falling oxygen.",
        "Hallem EA, Sternberg PW. PNAS 105:8038-8043 (2008).",
    ),
    "AQR": Curated(
        "Oxygen-sensing neuron with a process exposed to the body cavity, part of "
        "the RMG hub circuit.",
        "Macosko EZ, et al. Nature 458:1171-1175 (2009).",
    ),
    "AFD": Curated(
        "Principal thermosensory neuron, setting the temperature the animal "
        "prefers based on prior experience.",
        "Mori I, Ohshima Y. Nature 376:344-348 (1995).",
    ),
    "ADE": Curated(
        "Dopaminergic mechanosensory neuron of the anterior deirid, involved in "
        "sensing bacterial lawn texture and the resulting slowing response.",
        "Sawin ER, Ranganathan R, Horvitz HR. Neuron 26:619-631 (2000).",
    ),
    "CEP": Curated(
        "Dopaminergic mechanosensory neuron of the head, acting with ADE and PDE "
        "in the food-induced slowing response.",
        "Sawin ER, Ranganathan R, Horvitz HR. Neuron 26:619-631 (2000).",
    ),
    "DVA": Curated(
        "Interneuron with a stretch-sensitive, proprioceptive role: it reports the "
        "body's own bending back to the nervous system, so it is both an "
        "interneuron by anatomy and a sensor by function.",
        "Li W, Feng Z, Sternberg PW, Xu XZS. Nature 440:684-687 (2006).",
    ),
    "CAN": Curated(
        "Makes no chemical synapses at all, yet the animal dies without it. Its "
        "function remains unclear, which is why this repository records its role "
        "as unknown rather than guessing.",
        "Forrester WC, Garriga G. Development 124:1831-1843 (1997).",
    ),
    "NSM": Curated(
        "Serotonergic neurosecretory-motor neuron of the pharynx; detects "
        "bacterial food and slows the animal. Combines sensory, motor and "
        "neurosecretory function, so a single SIM role does not apply.",
        "Sawin ER, Ranganathan R, Horvitz HR. Neuron 26:619-631 (2000); "
        "Rhoades JL, et al. Cell 176:85-97 (2019).",
    ),
    "RIP": Curated(
        "The only direct connection between the pharyngeal nervous system and the "
        "rest of the animal.",
        "Albertson DG, Thomson JN. Phil Trans R Soc Lond B 275:299-325 (1976).",
    ),
    "HSN": Curated(
        "Hermaphrodite-specific serotonergic motor neuron driving egg laying.",
        "Desai C, Garriga G, McIntire SL, Horvitz HR. Nature 336:638-646 (1988).",
    ),
    "PVD": Curated(
        "Highly branched nociceptor responding to harsh touch and cold.",
        "Way JC, Chalfie M. Genes Dev 3:1823-1833 (1989); Chatzigeorgiou M, et al. "
        "Nat Neurosci 13:861-868 (2010).",
    ),
    "ALA": Curated(
        "Interneuron mediating a sleep-like quiescent state after cellular stress.",
        "Hill AJ, et al. Curr Biol 24:2399-2405 (2014).",
    ),
    "RIS": Curated(
        "Interneuron that induces developmentally timed sleep.",
        "Turek M, Lewandrowski I, Bringmann H. Curr Biol 23:2215-2223 (2013).",
    ),
}


def _descriptions() -> dict[str, dict[str, str]]:
    text = (files("worm.data") / "annotations" / "cell_descriptions.csv").read_text(
        encoding="utf-8"
    )
    return {r["cell_id"]: r for r in csv.DictReader(text.splitlines())}


def _nt_evidence() -> dict[str, dict[str, str]]:
    text = (files("worm.data") / "annotations" / "neurotransmitters.csv").read_text(
        encoding="utf-8"
    )
    return {r["cell_id"]: r for r in csv.DictReader(text.splitlines())}


NT_SHORT = {
    "acetylcholine": "ACh",
    "glutamate": "Glu",
    "gaba": "GABA",
    "dopamine": "DA",
    "serotonin": "5-HT",
    "octopamine": "OA",
    "tyramine": "TA",
    "betaine": "betaine",
    "unknown": "—",
    "none": "none",
}
ROLE_SHORT = {
    SIMRole.SENSORY: "sensory",
    SIMRole.INTER: "inter",
    SIMRole.MOTOR: "motor",
    SIMRole.MODULATORY: "modulatory",
    SIMRole.UNKNOWN: "unknown",
}
EVIDENCE_MARK = {
    "reporter_expression": "",
    "reporter_expression_dim_variable": " *",
    "uptake_not_synthesis": " †",
    "orphan_no_pathway_gene_detected": "",
}


def render(c: Connectome, dataset_id: str) -> str:
    desc = _descriptions()
    nt_ev = _nt_evidence()
    registry = {k: v.type_label for k, v in cell_registry().items()}
    deg = {d.cell_id: d for d in gr.degrees(c)}
    neurons = [x for x in c.cells if x.category is CellCategory.NEURON]

    def partners(cell_id: str) -> int:
        d = deg[cell_id]
        return d.in_degree + d.out_degree + d.gap_partners

    hubs = sorted(neurons, key=lambda x: -partners(x.id))[:25]

    lines: list[str] = []
    add = lines.append

    add("# Neuron reference")
    add("")
    add(
        f"Every cell in **`{dataset_id}`** ({c.organism}, {c.sex}, {c.stage}), with what its "
        "name stands for, what it does, and how heavily it is wired."
    )
    add("")
    add(
        "> **This file is generated.** Run `python tools/generate_neuron_reference.py` "
        "to rebuild it. Edit the generator, not this file."
    )
    add("")
    add("## How to read a name")
    add("")
    add(
        "*C. elegans* neuron names are acronyms, and the same cell has the same name in "
        "every animal — the worm's cell lineage is invariant, so `AVAL` is one specific "
        "cell, not a cell type. `ADEL` is the **A**nterior **DE**irid neuron, **L**eft; "
        "`ADER` is its right-hand partner. Both belong to class `ADE`."
    )
    add("")
    add("The trailing letters are usually position:")
    add("")
    add("| Suffix | Meaning |")
    add("| --- | --- |")
    for suffix, meaning in [
        ("`L` / `R`", "left / right of the bilateral pair"),
        ("`D` / `V`", "dorsal / ventral"),
        ("`DL`, `DR`, `VL`, `VR`", "one of four radial quadrants"),
        ("a number", "position along the body, head to tail (`VB1` … `VB11`)"),
    ]:
        add(f"| {suffix} | {meaning} |")
    add("")
    add(
        "The **lineage** column is the cell-division path from the fertilized egg. It is "
        "the reason named cells exist at all: division is invariant, so the same cell "
        "appears in the same place in every hermaphrodite."
    )
    add("")

    # ---- curated notes -------------------------------------------------
    add("## Notable cells")
    add("")
    add(
        "Hand-written notes, each with its own citation, for the cells you are most "
        "likely to meet first — the outliers in the figures and the ones the planned "
        "experiments depend on. Everything below this section is generated from data."
    )
    add("")
    present = {x.class_name for x in neurons if x.class_name}
    for klass in sorted(CURATED):
        note = CURATED[klass]
        members = sorted(x.id for x in neurons if x.class_name == klass)
        here = ", ".join(f"`{m}`" for m in members) if members else "_not in this dataset_"
        add(f"### {klass}")
        add("")
        add(note.summary)
        add("")
        add(f"- **In `{dataset_id}`:** {here}")
        add(f"- **Source:** {note.source}")
        add("")
    absent = sorted(set(CURATED) - present)
    if absent:
        add(
            f"> Not present in `{dataset_id}`: {', '.join(f'`{k}`' for k in absent)}. "
            f"This is a {c.scope} reconstruction — see [datasets.md](datasets.md)."
        )
        add("")

    # ---- hubs ----------------------------------------------------------
    add("## Most heavily connected neurons")
    add("")
    add(
        "Ranked by total number of partners in this dataset. Wiring counts are "
        "**measured**; role and neurotransmitter are annotations from other work."
    )
    add("")
    add("| Neuron | Class | Role | NT | Chem in | Chem out | Gap partners | Description |")
    add("| --- | --- | --- | --- | ---: | ---: | ---: | --- |")
    for cell in hubs:
        dr = deg[cell.id]
        add(
            f"| `{cell.id}` | {cell.class_name or '—'} | "
            f"{'/'.join(ROLE_SHORT[r] for r in cell.roles) or '—'} | "
            f"{_nt_cell(cell, nt_ev)} | {dr.in_degree} | {dr.out_degree} | {dr.gap_partners} | "
            f"{desc.get(cell.id, {}).get('classification', '')} |"
        )
    add("")

    # ---- full table ----------------------------------------------------
    add("## All neurons in this dataset")
    add("")
    add(f"{len(neurons)} neurons. Sorted by name.")
    add("")
    add(
        "**Where the columns disagree, that disagreement is real and is left visible.** "
        "`DVA` is listed with role *sensory* but described as a ring interneuron: it is "
        "an interneuron by anatomy that acts as a stretch receptor, reporting the body's "
        "own bending. `RIM` carries both *inter* and *motor* because the published label "
        "records a difference between studies. `RMG` is a major hub with no known "
        "transmitter at all. Flattening any of these would hide something true."
    )
    add("")
    add("Neurotransmitter marks: `*` dim and variable reporter expression; ")
    add("`†` taken up from neighbouring cells rather than synthesized; ")
    add("`—` no known transmitter pathway gene detected (a published result, not a gap).")
    add("")
    add("| Neuron | Class | Name stands for | Role | NT | Partners | Description | Lineage |")
    add("| --- | --- | --- | --- | --- | ---: | --- | --- |")
    for cell in sorted(neurons, key=lambda x: x.id):
        info = desc.get(cell.id, {})
        add(
            f"| `{cell.id}` | {cell.class_name or '—'} | {info.get('name_expansion', '')} | "
            f"{'/'.join(ROLE_SHORT[r] for r in cell.roles) or '—'} | "
            f"{_nt_cell(cell, nt_ev)} | {partners(cell.id)} | "
            f"{info.get('classification', '')} | `{info.get('lineage', '')}` |"
        )
    add("")

    # ---- non-neurons ---------------------------------------------------
    others = [x for x in c.cells if x.category is not CellCategory.NEURON]
    if others:
        add("## Non-neuronal cells in this dataset")
        add("")
        add(
            f"{len(others)} cells that are not neurons. Muscles are the nervous system's "
            "actual output, and glia are increasingly understood to shape neural "
            "function, so neither is filtered out."
        )
        add("")
        add("| Cell | Category | Partners | Description |")
        add("| --- | --- | ---: | --- |")
        for cell in sorted(others, key=lambda x: (str(x.category), x.id)):
            info = desc.get(cell.id, {})
            # WormAtlas records no free-text description for most muscles and glia,
            # so fall back to the cell registry's type label rather than a blank cell.
            described = (
                info.get("classification", "")
                or info.get("name_expansion", "")
                or registry.get(cell.id, "")
            )
            add(f"| `{cell.id}` | {cell.category} | {partners(cell.id)} | {described} |")
        add("")

    add("## Sources")
    add("")
    add(
        "| Column | Source | Confidence |\n"
        "| --- | --- | --- |\n"
        "| name expansion, lineage, description | WormAtlas cell listings |"
        " published annotation |\n"
        "| role | WormAtlas + Cook et al. 2019 groupings | published annotation |\n"
        "| class, neurotransmitter | Wang et al. 2024, eLife 12:RP95402 | published annotation |\n"
        f"| partner counts | {c.provenance.citation} | **measured** |"
    )
    add("")
    add(
        "Full citations and licences: [references.md](references.md). What each source "
        "does and does not establish: [model_assumptions.md](model_assumptions.md)."
    )
    add("")
    return "\n".join(lines)


def _nt_cell(cell: Cell, nt_ev: dict[str, dict[str, str]]) -> str:
    if not cell.neurotransmitters:
        return "—"
    mark = EVIDENCE_MARK.get(nt_ev.get(cell.id, {}).get("evidence", ""), "")
    return "/".join(NT_SHORT.get(str(n), str(n)) for n in cell.neurotransmitters) + mark


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--dataset", default="witvliet_2021_7")
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--check", action="store_true", help="fail if the file is out of date")
    args = p.parse_args(argv)

    c, _ = load(args.dataset)
    text = render(c, args.dataset)

    if args.check:
        current = args.out.read_text(encoding="utf-8") if args.out.exists() else ""
        if current != text:
            print(
                f"{args.out} is out of date; run python tools/generate_neuron_reference.py",
                file=sys.stderr,
            )
            return 1
        print(f"{args.out} is up to date")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text, encoding="utf-8", newline="")
    print(f"wrote {args.out.relative_to(REPO_ROOT)} ({len(text.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
