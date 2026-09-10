"""Developer tool: derive our committed annotation tables from published sources.

This script is **not** imported by the package and is not needed to use it. It
exists so that every annotation table in ``worm/data/`` is *reproducible* and
*reviewable* rather than hand-copied: run it, and ``git diff`` shows exactly what
changed and why.

Outputs
-------

``worm/data/cells.csv``
    Canonical cell registry: which cell identifiers exist, and whether each is a
    neuron, a muscle, a glial cell or something else. Importers need this to
    classify graph nodes, so it is *structural* data rather than an annotation
    overlay.

``worm/data/annotations/sim_roles.csv``
    Coarse functional role (sensory / interneuron / motor) per cell. This is an
    annotation overlay: it is not measured by any connectome reconstruction.

``worm/data/annotations/neurotransmitters.csv``
``worm/data/annotations/neuron_classes.csv``
    Neurotransmitter usage and anatomical neuron class, from the Wang et al. 2024
    neurotransmitter atlas.

Sources
-------

``all_cell_info.csv`` from the MIT-licensed OpenWorm *C. elegans* Connectome
Toolbox, which compiles cell listings from WormAtlas together with the cell
groupings used by Cook et al. 2019.

``elife-95402-supp2-v1.xlsx``, Supplementary File 2 of Wang C, Vidal B, Sural S,
et al., "A neurotransmitter atlas of *C. elegans* males and hermaphrodites",
eLife 12:RP95402 (2024) — CRISPR knock-in reporter expression for every
neurotransmitter pathway gene, scored across all 302 hermaphrodite neurons.

Honesty note
------------

The source gives each cell a free-text ``Type`` label such as
``"Layer 1 interneuron"`` or ``"Amphid, nociceptive"``. Collapsing those 61
labels onto our four-value :class:`~common.data.schemas.SIMRole` enum is *our*
interpretation, so the mapping is written out explicitly in :data:`LABEL_MAP`
below rather than being inferred by substring matching. Every label must appear
there or the script fails; a new upstream label can never be silently guessed at.

Where a label genuinely does not resolve to a Sensory/Inter/Motor assignment —
``"Pharyngeal polymodal neuron"`` and ``"Canal neuron"`` — we record
:attr:`SIMRole.UNKNOWN` and keep the original label, rather than inventing a
classification. See ``docs/model_assumptions.md``.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from worm.importers.naming import canonical_cell_id

REPO_ROOT = Path(__file__).resolve().parent.parent
WORM_DATA = REPO_ROOT / "worm" / "data"

ALL_CELL_INFO_URL = (
    "https://raw.githubusercontent.com/openworm/ConnectomeToolbox/main/cect/data/all_cell_info.csv"
)
ALL_CELL_INFO_SOURCE_ID = "openworm_cect_all_cell_info"
SOURCE_REF = "wormatlas_cook_2019_via_cect"


# Category values, matching common.data.schemas.CellCategory
NEURON, MUSCLE, GLIA, OTHER = "neuron", "muscle", "glia", "other"
# Role values, matching common.data.schemas.SIMRole
SENSORY, INTER, MOTOR, MODULATORY, UNKNOWN = (
    "sensory",
    "inter",
    "motor",
    "modulatory",
    "unknown",
)


@dataclass(frozen=True)
class LabelMapping:
    category: str
    roles: tuple[str, ...]
    note: str = ""


def _n(*roles: str) -> LabelMapping:
    return LabelMapping(NEURON, roles)


#: Explicit mapping from the upstream free-text ``Type`` label to our category and
#: role enums. Exhaustive by construction: an unmapped label aborts the run.
LABEL_MAP: dict[str, LabelMapping] = {
    # --- sensory neurons -------------------------------------------------
    "Amphid": _n(SENSORY),
    "Amphid, nociceptive": _n(SENSORY),
    "Cephalic": _n(SENSORY),
    "Mechanosensory": _n(SENSORY),
    "Touch": _n(SENSORY),
    "Phasmid": _n(SENSORY),
    "O2, CO2, social signals, touch": _n(SENSORY),
    "Male sensory neuron": _n(SENSORY),
    "Male head sensory neuron": _n(SENSORY),
    # --- interneurons ----------------------------------------------------
    "Layer 1 interneuron": _n(INTER),
    "Layer 2 interneuron": _n(INTER),
    "Layer 3 interneuron": _n(INTER),
    "Category 4 interneuron": _n(INTER),
    "Pharyngeal interneuron": _n(INTER),
    "Male interneuron": _n(INTER),
    "Male head interneuron": _n(INTER),
    "Linker to pharynx": LabelMapping(
        NEURON, (INTER,), "RIPL/RIPR: the only cells linking somatic and pharyngeal systems"
    ),
    # --- motor neurons ---------------------------------------------------
    "Ventral cord motor neuron": _n(MOTOR),
    "Head motor neuron": _n(MOTOR),
    "Sublateral motor neuron": _n(MOTOR),
    "Pharyngeal motor neuron": _n(MOTOR),
    "Hermaphrodite specific motor neuron": _n(MOTOR),
    # --- neurons whose published label asserts two roles -----------------
    "Sublateral motor neuron; interneuron in White et al., 1986": LabelMapping(
        NEURON, (MOTOR, INTER), "the source label itself records disagreement between studies"
    ),
    "Layer 1 interneuron; motorneuron in White et al., 1986": LabelMapping(
        NEURON, (INTER, MOTOR), "the source label itself records disagreement between studies"
    ),
    # --- neurons we decline to classify ----------------------------------
    "Pharyngeal polymodal neuron": LabelMapping(
        NEURON,
        (UNKNOWN,),
        "MI, NSML/R, MCL/R combine sensory, motor and neurosecretory function; "
        "'polymodal' does not resolve to a single SIM role and we do not invent one",
    ),
    "Canal neuron": LabelMapping(
        NEURON,
        (UNKNOWN,),
        "CANL/CANR make no chemical synapses; essential for survival, function unclear",
    ),
    # --- muscle ----------------------------------------------------------
    "Main body muscle": LabelMapping(MUSCLE, ()),
    "Head muscle": LabelMapping(MUSCLE, ()),
    "Unspecified body wall muscle": LabelMapping(MUSCLE, ()),
    "Pharyngeal muscle": LabelMapping(MUSCLE, ()),
    "Vulval muscle": LabelMapping(MUSCLE, ()),
    "Uterine muscle": LabelMapping(MUSCLE, ()),
    "Anal/sphincter muscle": LabelMapping(MUSCLE, ()),
    "Intestinal muscles": LabelMapping(MUSCLE, ()),
    "Diagonal muscle (male specific)": LabelMapping(MUSCLE, ()),
    "Anterior oblique (male specific)": LabelMapping(MUSCLE, ()),
    "Posterior oblique (male specific)": LabelMapping(MUSCLE, ()),
    "Gubernacular erector (male specific)": LabelMapping(MUSCLE, ()),
    "Gubernacular retractor (male specific)": LabelMapping(MUSCLE, ()),
    "Caudal longitudinal muscle (male specific)": LabelMapping(MUSCLE, ()),
    "Anterior inner longitudinal muscle (male specific)": LabelMapping(MUSCLE, ()),
    "Posterior inner longitudinal muscle (male specific)": LabelMapping(MUSCLE, ()),
    "Posterior outer longitudinal muscle (male specific)": LabelMapping(MUSCLE, ()),
    "Dorsal spicule protractor (male specific)": LabelMapping(MUSCLE, ()),
    "Ventral spicule protractor (male specific)": LabelMapping(MUSCLE, ()),
    "Dorsal spicule retractor (male specific)": LabelMapping(MUSCLE, ()),
    "Ventral spicule retractor (male specific)": LabelMapping(MUSCLE, ()),
    # --- glia ------------------------------------------------------------
    "Sheath cell other than amphid sheath and phasmid": LabelMapping(GLIA, ()),
    "Pharyngeal glial cell": LabelMapping(GLIA, ()),
    # --- everything else -------------------------------------------------
    "GLR cell": LabelMapping(
        OTHER, (), "mesodermal cells gap-junctioned to head muscle and RME motor neurons"
    ),
    "Marginal cell of the pharynx": LabelMapping(OTHER, ()),
    "Pharyngeal epithelium": LabelMapping(OTHER, ()),
    "Pharyngeal basement membrane": LabelMapping(OTHER, ()),
    "Male ray structural cell": LabelMapping(OTHER, ()),
    "Excretory cell": LabelMapping(OTHER, ()),
    "Excretory gland": LabelMapping(OTHER, ()),
    "Head mesodermal cell": LabelMapping(OTHER, ()),
    "Hypodermis": LabelMapping(OTHER, ()),
    "Intestine": LabelMapping(OTHER, ()),
    "Proctodeum (male specific)": LabelMapping(OTHER, ()),
    "Gonad (male specific)": LabelMapping(OTHER, ()),
}

CELLS_COLUMNS = ("cell_id", "category", "type_label", "sex_specific", "source_ref")
SIM_ROLES_COLUMNS = ("cell_id", "role", "type_label", "source_ref")
DESCRIPTIONS_COLUMNS = ("cell_id", "name_expansion", "lineage", "classification", "source_ref")
POLARITY_COLUMNS = ("pre", "post", "sign", "basis", "primary_nt", "source_row", "source_ref")

#: WormAtlas leaves unfilled fields as this literal string. It means "not recorded",
#: which is different from an empty description, so it is dropped rather than stored.
_PLACEHOLDER = "To be added"
NT_COLUMNS = ("cell_id", "neurotransmitter", "evidence", "source_label", "source_row", "source_ref")
CLASS_COLUMNS = ("cell_id", "class_name", "source_row", "source_ref")


# ---------------------------------------------------------------------------
# Wang et al. 2024 neurotransmitter atlas
# ---------------------------------------------------------------------------

NT_ATLAS_URL = (
    "https://raw.githubusercontent.com/openworm/ConnectomeToolbox/main/cect/data/"
    "elife-95402-supp2-v1.xlsx"
)
NT_SOURCE_REF = "wang_2024_nt_atlas"

# Column positions in Supplementary File 2 (header on row 4, data from row 5).
_NT_COL_CLASS, _NT_COL_NEURON, _NT_COL_ASSIGNMENT, _NT_COL_COMMENT = 1, 2, 20, 23
_NT_HEADER_ROW = 4

# Neurotransmitter abbreviations, per the legend at the foot of the sheet.
_NT_TOKEN = {
    "ach": "acetylcholine",
    "glu": "glutamate",
    "gaba": "gaba",
    "da": "dopamine",
    "5-ht": "serotonin",
    "octopamine": "octopamine",
    "tyramine": "tyramine",
    "betaine": "betaine",
}

#: Wang et al. cannot reliably distinguish DB1 from DB3 and label those rows
#: "DB1/3" and "DB3/1". Both are assigned the same transmitter, so resolving the
#: pair in a fixed order cannot produce a wrong assignment for either cell.
_NT_CELL_ALIASES = {"DB1/3": "DB1", "DB3/1": "DB3"}


def _classify_nt(raw: str) -> tuple[str | None, str, str]:
    """Interpret one assignment cell.

    Returns ``(neurotransmitter, evidence, note)``. ``None`` means the source
    reports no transmitter for this neuron, which we record as an explicit gap
    rather than omitting the row.

    Three distinctions in the source are preserved because they are biologically
    real and would otherwise be flattened away:

    * a leading ``*`` marks *dim and variable* reporter expression, i.e. genuine
      but weaker evidence than a clean positive;
    * ``(uptake)`` marks a cell that does not synthesize the transmitter but takes
      it up from its neighbours — it can still release it, but the mechanism and
      the genetics differ;
    * ``unknown (orphan…)`` marks a neuron for which no known transmitter pathway
      gene is expressed. That is a *result*, not a missing measurement.
    """
    text = raw.strip()
    low = text.lower()

    if low.startswith("unknown") or low.startswith("*unknown"):
        return None, "orphan_no_pathway_gene_detected", text

    dim = text.startswith("*")
    uptake = "(uptake)" in low

    body = low.lstrip("*").replace("- new", "").replace("(uptake)", "").strip()
    token = _NT_TOKEN.get(body)
    if token is None:
        raise SystemExit(
            f"unrecognised neurotransmitter assignment {raw!r}. Add it to _NT_TOKEN "
            "explicitly rather than guessing; a wrong transmitter changes the sign of "
            "every synapse the cell makes."
        )

    if uptake:
        evidence = "uptake_not_synthesis"
    elif dim:
        evidence = "reporter_expression_dim_variable"
    else:
        evidence = "reporter_expression"
    return token, evidence, text


def parse_nt_atlas(path: Path) -> tuple[list[tuple[str, ...]], list[tuple[str, ...]]]:
    from openpyxl import load_workbook

    wb = load_workbook(path, read_only=True, data_only=True)
    rows = list(wb.worksheets[0].iter_rows(values_only=True))
    wb.close()

    nt_rows: list[tuple[str, ...]] = []
    class_rows: list[tuple[str, ...]] = []
    current_class: str | None = None

    for row_no, row in enumerate(rows[_NT_HEADER_ROW:], start=_NT_HEADER_ROW + 1):
        neuron = row[_NT_COL_NEURON]
        if not neuron:
            continue
        cell_id = str(neuron).strip()
        cell_id = _NT_CELL_ALIASES.get(cell_id, cell_id)

        klass = row[_NT_COL_CLASS]
        if klass:
            current_class = str(klass).strip()
        if current_class:
            class_rows.append((cell_id, current_class, str(row_no), NT_SOURCE_REF))

        assignment = row[_NT_COL_ASSIGNMENT]
        if not assignment:
            continue
        nt, evidence, label = _classify_nt(str(assignment))
        nt_rows.append(
            (cell_id, nt or UNKNOWN, evidence, label, str(row_no), NT_SOURCE_REF)
        )

    nt_rows.sort(key=lambda t: (t[0], t[1]))
    class_rows.sort(key=lambda t: t[0])
    return nt_rows, class_rows


# ---------------------------------------------------------------------------
# Fenyves et al. 2020 synaptic polarity predictions
# ---------------------------------------------------------------------------

FENYVES_URL = (
    "https://journals.plos.org/ploscompbiol/article/file?"
    "id=10.1371/journal.pcbi.1007974.s003&type=supplementary"
)
FENYVES_FILENAME = "pcbi.1007974.s003.xlsx"
FENYVES_SOURCE_REF = "fenyves_2020_polarity"
FENYVES_SHEET = "5. Sign prediction"

#: Sheet layout: two header rows, data from row 3 (1-based).
_FEN_HEADER_ROWS = 2
_FEN_COL_PRE, _FEN_COL_NT1, _FEN_COL_POST = 0, 1, 3
_FEN_COL_TYPE, _FEN_COL_POLARITY = 5, 16

#: How the source's polarity strings map onto our Sign enum.
#:
#: ``complex`` means the postsynaptic cell expresses BOTH excitatory and inhibitory
#: receptors for the presynaptic transmitter, so the net effect is genuinely
#: undetermined by this method — not merely unmeasured. ``no pred`` means no
#: receptor match was found at all. The two are different findings and are kept apart.
_FEN_SIGN = {
    "+": ("excitatory", "receptor_match"),
    "-": ("inhibitory", "receptor_match"),
    "complex": ("mixed", "both_excitatory_and_inhibitory_receptors"),
    "no pred": ("unknown", "no_receptor_match_found"),
}


def parse_fenyves_polarity(path: Path) -> list[tuple[str, ...]]:
    """Per-connection predicted synaptic polarity.

    **These are predictions, not measurements.** They combine the presynaptic
    neurotransmitter with postsynaptic ionotropic receptor gene expression: if the
    target expresses an excitatory receptor for that transmitter and no inhibitory
    one, the synapse is predicted excitatory, and vice versa.

    Every row of the source is emitted, including the ones with no prediction, so
    that "Fenyves considered this connection and could not call it" stays
    distinguishable from "this connection is not in Fenyves at all".

    Scope limit worth knowing before relying on this: the source covers
    **interneuronal connections only**. It contains no neuromuscular junctions, so
    it supplies no polarity for the synapses that actually drive muscle.
    """
    from openpyxl import load_workbook

    wb = load_workbook(path, read_only=True, data_only=True)
    if FENYVES_SHEET not in wb.sheetnames:
        raise SystemExit(f"{path.name}: expected a {FENYVES_SHEET!r} sheet")
    rows = list(wb[FENYVES_SHEET].iter_rows(values_only=True))
    wb.close()

    out: list[tuple[str, ...]] = []
    for row_no, row in enumerate(rows[_FEN_HEADER_ROWS:], start=_FEN_HEADER_ROWS + 1):
        pre, post = row[_FEN_COL_PRE], row[_FEN_COL_POST]
        if not pre or not post:
            continue
        edge_type = str(row[_FEN_COL_TYPE] or "").strip().lower()
        if edge_type != "chemical":
            continue  # polarity is meaningless for a gap junction
        raw = str(row[_FEN_COL_POLARITY] or "").strip()
        if raw not in _FEN_SIGN:
            raise SystemExit(
                f"{path.name} row {row_no}: unrecognised polarity {raw!r}. Add it to "
                "_FEN_SIGN explicitly; guessing a synapse's sign is exactly what this "
                "project must not do."
            )
        sign, basis = _FEN_SIGN[raw]
        # Canonicalise, exactly as the connectome importers do.
        #
        # Fenyves zero-pads the ventral cord motor neurons -- DB01, VA07, AS03 --
        # where Cook and Witvliet write DB1, VA7, AS3. Emitting the raw labels made
        # 634 of 3638 rows unmatchable, and 582 of those carried an actual sign, so
        # a quarter of every polarity prediction we have was being discarded. The
        # loss fell almost entirely on the ventral nerve cord: every AS, DA, DB, VA,
        # VB, DD and VD cell, which is to say the entire locomotor circuit.
        #
        # It failed quietly because the overlay reported it as "739 table entries
        # not in this dataset", which is indistinguishable from the legitimate case
        # of a table covering a different animal. See test_annotations.py, which now
        # asserts the match rate rather than trusting the summary line.
        out.append(
            (
                canonical_cell_id(str(pre).strip()),
                canonical_cell_id(str(post).strip()),
                sign,
                basis,
                str(row[_FEN_COL_NT1] or "").strip(),
                str(row_no),
                FENYVES_SOURCE_REF,
            )
        )
    out.sort(key=lambda t: (t[0], t[1]))
    return out


def fetch_text(url: str, cache: Path | None) -> str:
    if cache is not None and cache.exists():
        return cache.read_text(encoding="utf-8-sig")
    with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310 - fixed https URL
        text = r.read().decode("utf-8-sig")
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(text, encoding="utf-8", newline="")
    return text


def _clean(value: str | None) -> str:
    """Normalize a WormAtlas free-text field, dropping its 'not recorded' placeholder."""
    text = (value or "").strip().strip("-").strip()
    return "" if not text or _PLACEHOLDER in text else " ".join(text.split())


def parse_all_cell_info(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in csv.DictReader(io.StringIO(text)):
        cell_id = (row.get("Cell name") or "").strip()
        label = (row.get("Type") or "").strip()
        if not cell_id or not label:
            continue
        if cell_id in seen:
            raise ValueError(f"duplicate cell id in source: {cell_id!r}")
        seen.add(cell_id)
        rows.append(
            {
                "cell_id": cell_id,
                "type_label": label,
                "name_expansion": _clean(row.get("Name details")),
                "lineage": _clean(row.get("Lineage")),
                "classification": _clean(row.get("Classification")),
            }
        )
    return rows


def build_descriptions(rows: list[dict[str, str]]) -> list[tuple[str, ...]]:
    """Name etymology, embryonic lineage and functional description, per cell.

    *C. elegans* neuron names are acronyms — ``ADEL`` is "Anterior DEirid neuron
    Left" — so the expansion is genuinely informative rather than decorative.
    The lineage string is the invariant cell-division path from the zygote, which
    is why every hermaphrodite has the same named cells at all.

    Rows with no recorded text are omitted rather than stored blank: a missing
    description is a gap in WormAtlas, and should read as one.
    """
    out: list[tuple[str, ...]] = []
    for r in sorted(rows, key=lambda r: r["cell_id"]):
        if not any((r["name_expansion"], r["lineage"], r["classification"])):
            continue
        out.append(
            (
                r["cell_id"],
                r["name_expansion"],
                r["lineage"],
                r["classification"],
                SOURCE_REF,
            )
        )
    return out


def build_tables(rows: list[dict[str, str]]) -> tuple[list[tuple[str, ...]], list[tuple[str, ...]]]:
    unmapped = sorted({r["type_label"] for r in rows} - set(LABEL_MAP))
    if unmapped:
        raise SystemExit(
            "Upstream introduced Type labels this script does not know how to map:\n  "
            + "\n  ".join(repr(u) for u in unmapped)
            + "\n\nAdd them to LABEL_MAP explicitly. Do not guess by substring: a wrong "
            "sensory/motor assignment silently corrupts every downstream circuit."
        )

    cells: list[tuple[str, ...]] = []
    sim_roles: list[tuple[str, ...]] = []
    for r in sorted(rows, key=lambda r: r["cell_id"]):
        cell_id, label = r["cell_id"], r["type_label"]
        m = LABEL_MAP[label]
        sex = "male" if "male specific" in label.lower() or label.startswith("Male") else ""
        cells.append((cell_id, m.category, label, sex, SOURCE_REF))
        for role in m.roles:
            sim_roles.append((cell_id, role, label, SOURCE_REF))
    sim_roles.sort(key=lambda t: (t[0], t[1]))
    return cells, sim_roles


def write_csv(path: Path, columns: tuple[str, ...], rows: list[tuple[str, ...]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO(newline="")
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(columns)
    w.writerows(rows)
    path.write_text(buf.getvalue(), encoding="utf-8", newline="")
    print(f"wrote {path.relative_to(REPO_ROOT)} ({len(rows)} rows)")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--cache", type=Path, default=WORM_DATA / "raw" / "all_cell_info.csv")
    p.add_argument("--no-cache", action="store_true", help="always re-download")
    args = p.parse_args(argv)

    text = fetch_text(ALL_CELL_INFO_URL, None if args.no_cache else args.cache)
    rows = parse_all_cell_info(text)
    cells, sim_roles = build_tables(rows)

    write_csv(WORM_DATA / "cells.csv", CELLS_COLUMNS, cells)
    write_csv(WORM_DATA / "annotations" / "sim_roles.csv", SIM_ROLES_COLUMNS, sim_roles)
    write_csv(
        WORM_DATA / "annotations" / "cell_descriptions.csv",
        DESCRIPTIONS_COLUMNS,
        build_descriptions(rows),
    )

    nt_path = WORM_DATA / "raw" / "elife-95402-supp2-v1.xlsx"
    if not nt_path.exists():
        print(f"downloading {NT_ATLAS_URL}")
        nt_path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(NT_ATLAS_URL, timeout=120) as r:  # noqa: S310
            nt_path.write_bytes(r.read())
    nt_rows, class_rows = parse_nt_atlas(nt_path)
    write_csv(WORM_DATA / "annotations" / "neurotransmitters.csv", NT_COLUMNS, nt_rows)
    write_csv(WORM_DATA / "annotations" / "neuron_classes.csv", CLASS_COLUMNS, class_rows)

    fen_path = WORM_DATA / "raw" / FENYVES_FILENAME
    if not fen_path.exists():
        print(f"downloading {FENYVES_URL}")
        fen_path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(FENYVES_URL, timeout=180) as r:  # noqa: S310
            fen_path.write_bytes(r.read())
    polarity_rows = parse_fenyves_polarity(fen_path)
    write_csv(
        WORM_DATA / "annotations" / "polarity_fenyves2020.csv",
        POLARITY_COLUMNS,
        polarity_rows,
    )

    known = {c[0] for c in cells}
    stray = sorted({r[0] for r in nt_rows} - known)
    if stray:
        raise SystemExit(
            f"neurotransmitter atlas names not in the cell registry: {stray}. "
            "Add a canonicalization rule or an alias; do not drop them."
        )

    by_cat: dict[str, int] = {}
    for _, cat, *_ in cells:
        by_cat[cat] = by_cat.get(cat, 0) + 1
    print("\ncells by category: " + ", ".join(f"{k}={v}" for k, v in sorted(by_cat.items())))
    unresolved = sorted({c for c, r, *_ in sim_roles if r == UNKNOWN})
    print(f"cells with an unresolved SIM role ({len(unresolved)}): {unresolved}")

    nt_counts: dict[str, int] = {}
    for _, nt, *_ in nt_rows:
        nt_counts[nt] = nt_counts.get(nt, 0) + 1
    print("neurotransmitters: " + ", ".join(f"{k}={v}" for k, v in sorted(nt_counts.items())))
    orphans = sorted({r[0] for r in nt_rows if r[1] == UNKNOWN})
    print(f"neurons with no detected transmitter pathway gene ({len(orphans)}): {orphans}")

    pol: dict[str, int] = {}
    for r in polarity_rows:
        pol[r[2]] = pol.get(r[2], 0) + 1
    print("predicted synapse polarity: " + ", ".join(f"{k}={v}" for k, v in sorted(pol.items())))
    called = pol.get("excitatory", 0) + pol.get("inhibitory", 0)
    print(
        f"  a definite sign for {called}/{len(polarity_rows)} connections "
        f"({called / max(len(polarity_rows), 1):.1%}); "
        f"E:I = {pol.get('excitatory', 0) / max(pol.get('inhibitory', 1), 1):.2f}:1"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
