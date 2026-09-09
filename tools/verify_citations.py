"""Check every citation in the repository against CrossRef.

    python tools/verify_citations.py check          # verify what has a DOI
    python tools/verify_citations.py find           # look up citations with no DOI
    python tools/verify_citations.py check --json

Why this exists
---------------

This project's rule is that no biological claim gets stated without a source. A
citation written from memory breaks that rule just as thoroughly as an invented
parameter does, and it is harder to notice -- a plausible-looking author list
reads exactly like a correct one.

So citations are treated like every other kind of data here: resolved against an
authoritative source and asserted, rather than trusted. CrossRef is that source.
Given a DOI it returns the registered title, author list, journal, volume, pages
and year; this compares ours against it and reports every mismatch.

``find`` handles the other direction -- a free-text citation with no DOI is
searched for, so the DOI can be filled in and the entry becomes checkable.

Responses are cached under the scratch directory so repeated runs and CI do not
hammer a free public service.

Note on the CrossRef "polite pool": CrossRef gives faster service to requests
that include a contact address. This tool does **not** embed one, and will not
send any address unless you explicitly set ``CROSSREF_MAILTO`` yourself.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import tomllib
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CACHE_DIR = REPO_ROOT / ".citation-cache"
API = "https://api.crossref.org/works"
TIMEOUT = 45

#: Files scanned for `doi:[...](https://doi.org/...)` or `doi = "..."` entries.
SOURCE_FILES = [
    REPO_ROOT / "worm" / "data" / "sources.toml",
    REPO_ROOT / "docs" / "references.md",
    REPO_ROOT / "DATA_LICENSES.md",
    REPO_ROOT / "CITATION.cff",
    REPO_ROOT / "worm" / "neural" / "parameters.toml",
    REPO_ROOT / "worm" / "data" / "circuits" / "circuits.json",
]

_DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]*[A-Za-z0-9]")
_YEAR_RE = re.compile(r"\((\d{4})\)")

#: Places where CrossRef's own registered metadata is wrong and ours is right.
#:
#: These are not excuses -- each was checked against the article itself before being
#: listed, and the reason is recorded so a future reader can re-check rather than
#: take it on trust. Keyed by DOI, mapping the field to why we differ.
KNOWN_CROSSREF_ERRATA: dict[str, dict[str, str]] = {
    "10.1098/rstb.1976.0085": {
        "authors missing from our citation": (
            "CrossRef registers the second author as 'Thompson J. N.'. The correct "
            "spelling is Thomson -- the same person is registered as 'Thomson' on "
            "White et al. (10.1098/rstb.1976.0086) in the same journal volume, and "
            "on White et al. 1986. We keep Thomson."
        ),
        "names in our citation that CrossRef does not list": (
            "The same erratum seen from the other side: because CrossRef lists "
            "Thompson, our correct spelling Thomson looks like a name CrossRef does "
            "not know. Verified against the article and against the same author's "
            "other 1976 and 1986 papers, which spell it Thomson."
        ),
    },
    # The same OCR artefact appears on both McIntire papers.
    "10.1038/336638a0": {
        "authors missing from our citation": (
            "CrossRef registers 'Mclntire' with a lowercase L for the capital I in "
            "McIntire -- an OCR artefact. We keep McIntire."
        ),
    },
    "10.1038/364337a0": {
        "authors missing from our citation": (
            "CrossRef registers 'Mclntire' with a lowercase L for the capital I in "
            "McIntire -- an OCR artefact. We keep McIntire."
        ),
    },
    "10.1242/dev.124.9.1831": {
        "year": (
            "CrossRef dates the whole Development back-catalogue to 1995-02-01 as a "
            "bulk-registration artefact. The DOI itself encodes volume 124 issue 9, "
            "which is 1997. We keep 1997."
        ),
    },
    "10.1038/nature04538": {
        "authors missing from our citation": (
            "CrossRef records the last author's family name as 'Shawn Xu'; the "
            "article bylines him X. Z. Shawn Xu, which we render 'Xu XZS'."
        ),
    },
}


@dataclass
class Citation:
    """One citation found in the repository."""

    where: str
    text: str
    doi: str | None = None


@dataclass
class Mismatch:
    citation: Citation
    field: str
    ours: str
    crossref: str

    def __str__(self) -> str:
        return (
            f"{self.citation.where}\n"
            f"  {self.field}:\n"
            f"    ours     : {self.ours}\n"
            f"    crossref : {self.crossref}"
        )


@dataclass
class Result:
    checked: int = 0
    unresolvable: list[Citation] = field(default_factory=list)
    mismatches: list[Mismatch] = field(default_factory=list)


# ---------------------------------------------------------------------------
# CrossRef
# ---------------------------------------------------------------------------


def _fetch(url: str, cache_key: str) -> dict | None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cached = CACHE_DIR / f"{re.sub(r'[^A-Za-z0-9]+', '_', cache_key)[:120]}.json"
    if cached.exists():
        return json.loads(cached.read_text(encoding="utf-8"))

    import os

    agent = "isaac-bioresearch/0.1 (https://github.com/Stellarica-Interactive/isaac-bioresearch-playground)"
    # Opt-in only: no contact address is sent unless the user sets this themselves.
    if mailto := os.environ.get("CROSSREF_MAILTO"):
        agent += f" mailto:{mailto}"
    req = urllib.request.Request(url, headers={"User-Agent": agent})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:  # noqa: S310 - fixed https host
            payload = json.load(r)
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
        return None
    cached.write_text(json.dumps(payload), encoding="utf-8", newline="")
    return payload


def by_doi(doi: str) -> dict | None:
    payload = _fetch(f"{API}/{urllib.parse.quote(doi)}", f"doi_{doi}")
    return payload.get("message") if payload else None


def search(text: str) -> dict | None:
    query = urllib.parse.quote(text[:400])
    payload = _fetch(f"{API}?query.bibliographic={query}&rows=1", f"q_{text[:100]}")
    items = (payload or {}).get("message", {}).get("items", [])
    return items[0] if items else None


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------


def _fold(text: str) -> str:
    """Strip accents and case so `Sőti` matches `Soti`, which our ASCII TOML uses."""
    decomposed = unicodedata.normalize("NFKD", text)
    return "".join(c for c in decomposed if not unicodedata.combining(c)).lower()


def surnames(record: dict) -> list[str]:
    return [a.get("family", "") for a in record.get("author", []) if a.get("family")]


def _year(record: dict) -> str | None:
    parts = record.get("issued", {}).get("date-parts", [[None]])
    return str(parts[0][0]) if parts and parts[0] and parts[0][0] else None


def compare(citation: Citation, record: dict) -> list[Mismatch]:
    """Compare a free-text citation against the registered record.

    Deliberately loose on formatting and strict on facts: it checks that every
    registered surname appears in our text and that no *extra* surname-looking
    token has crept in, plus the year. Formatting styles vary; a fabricated
    co-author does not.
    """
    out: list[Mismatch] = []
    ours = _fold(citation.text)

    registered = surnames(record)
    abbreviated = "et al" in ours

    if abbreviated:
        # "Witvliet D, et al." is a legitimate short form; only the first author is
        # being claimed, so only the first author is checked.
        if registered and _fold(registered[0]) not in ours:
            out.append(
                Mismatch(citation, "first author", citation.text[:60], registered[0])
            )
    else:
        missing = [s for s in registered if _fold(s) not in ours]
        if missing:
            out.append(
                Mismatch(
                    citation,
                    "authors missing from our citation",
                    ", ".join(missing),
                    "; ".join(registered),
                )
            )

    # Catch invented co-authors: a capitalised surname-like token in our text that
    # is in no registered author and is not an ordinary word of the title/journal.
    title = " ".join(record.get("title", []) + record.get("container-title", []))
    allowed = {_fold(w) for w in re.findall(r"[A-Za-zÀ-ž']+", title)}
    for name in registered:
        # A hyphenated surname such as Ripoll-Sánchez must match both as a whole
        # and as its parts, since our text may be tokenised either way.
        allowed.add(_fold(name))
        allowed.update(_fold(part) for part in re.split(r"[-\s]+", name) if part)
    allowed |= _STOPWORDS
    author_part = citation.text.split(".")[0] if "." in citation.text else citation.text
    suspicious = [
        w
        for w in re.findall(r"\b[A-ZÀ-Ž][a-zà-ž]{3,}\b", author_part)
        if _fold(w) not in allowed
    ]
    if suspicious:
        out.append(
            Mismatch(
                citation,
                "names in our citation that CrossRef does not list",
                ", ".join(sorted(set(suspicious))),
                "; ".join(registered),
            )
        )

    year = _year(record)
    ours_year = _YEAR_RE.search(citation.text)
    if year and ours_year and ours_year.group(1) != year:
        out.append(Mismatch(citation, "year", ours_year.group(1), year))

    return out


_STOPWORDS = {
    "the", "and", "for", "with", "from", "using", "data", "gene", "expression",
    "nature", "neuron", "cell", "elife", "science", "supplementary", "information",
    "supp", "file", "table", "reader", "connectome", "connectomes", "prediction",
    "predicts", "sign", "balance", "polarity", "synaptic", "chemical", "synapse",
    "neuronal", "network", "whole", "animal", "both", "sexes", "across",
    "development", "reveal", "principles", "brain", "maturation", "atlas",
    "males", "hermaphrodites", "neurotransmitter", "computational", "biology",
    "physical", "review", "journal", "neuroscience", "proceedings", "national",
    "academy", "sciences", "current", "frontiers", "communications", "plos",
    "one", "dynamic", "simulation", "nematode", "withdrawal",
    "circuit", "predictions", "concerning", "function", "behavioral", "criteria",

    "active", "currents", "regulate", "sensitivity", "range", "neurons",
    "olfactory", "fire", "calcium", "mediated", "none", "action", "potentials",
    "quest", "hits", "plateau", "enteric", "motor", "synchronized", "underlying",
    "defecation", "program", "structure", "nervous", "system", "neural",
    "circuits", "touch", "dual", "mechanosensory", "chemosensory",
    "odorant", "selective", "genes", "mediate", "olfaction", "fundamental",
    "role", "pirouettes", "chemotaxis", "functional", "asymmetry", "taste",
    "hub", "spoke", "drives", "pheromone", "attraction",
    "social", "behaviour", "imbalancing", "junctions", "reduce", "backward",
    "activity", "bias", "forward", "locomotion", "signal", "propagation",
    "neuropeptidergic", "biophysical", "modeling", "single", "ion", "dynamics",
    "low", "dimensional", "functionality", "complex", "neurosensory",
    "integration", "multistability", "long", "timescale", "transients",
    "encoded", "model", "gaba", "two", "acetylcholine", "receptors", "junction",
    "gabaergic", "required", "acid", "decarboxylase", "transmission",
    "stretch", "receptor", "proprioceptive", "essential",
    "survival", "serotonergic", "sleep", "quiescent", "state", "cellular",
    "stress", "harsh", "cold", "nociceptor", "branched", "egg", "laying",
    "specific", "direct", "connection", "pharyngeal", "pharynx", "somatic",
    "thermosensory", "temperature", "preferred", "experience", "dopaminergic",
    "bacterial", "lawn", "texture", "slowing", "response", "food", "induced",
    "oxygen", "carbon", "dioxide", "sensing", "body", "cavity", "exposed",
    "compartmentalized", "within", "highly", "connected", "ring",
    "interneuron", "integrating", "sensory", "input", "head", "output",
    "asymmetric", "left", "right", "odorants", "detects", "attractive",
    "volatile", "odours", "avoidance", "repulsive", "premotor", "command",
    "associated", "anterior", "posterior", "gentle", "ablating",
    "abolishes", "reversal", "acceleration", "born", "embryonically", "works",
    "polymodal", "responds", "noxious", "chemicals", "high", "osmolarity",
    "nose", "example", "why", "label", "simplification",
    "principal", "salt", "pair", "increases", "decreases", "members", "class",
    "interchangeable", "notable", "genuine", "all", "or", "common", "summary",
    "omits", "first", "layer", "amphid", "major", "point", "downstream",
    "implicated", "turn", "run", "decision", "acting", "largely", "opposition",
    "turning", "receiving", "convergent", "pathway", "coupled", "records",
    "difference", "between", "studies", "whether", "several", "electrically",
    "aggregates", "them", "falling", "process", "innervate", "innervates",
    "muscles", "muscle", "vent",
}


# ---------------------------------------------------------------------------
# Collecting citations
# ---------------------------------------------------------------------------


def from_sources_toml() -> list[Citation]:
    path = REPO_ROOT / "worm" / "data" / "sources.toml"
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    out = []
    for key, entry in data.items():
        if not isinstance(entry, dict) or "citation" not in entry:
            continue
        out.append(
            Citation(
                where=f"worm/data/sources.toml [{key}]",
                text=str(entry["citation"]),
                doi=str(entry["doi"]) if entry.get("doi") else None,
            )
        )
    return out


def from_curated_notes() -> list[Citation]:
    """The hand-written neuron notes, which are the largest body of free text."""
    from tools.generate_neuron_reference import CURATED

    out = []
    for klass, note in sorted(CURATED.items()):
        parts = [p.strip() for p in note.source.split(";") if len(p.strip()) > 20]
        for i, part in enumerate(parts):
            # A note may cite two papers; the recorded DOI belongs to the last one,
            # which is the primary source for the claim. Earlier ones fall back to
            # a bibliographic search.
            doi = note.doi if (note.doi and i == len(parts) - 1) else None
            out.append(Citation(where=f"CURATED[{klass!r}]", text=part, doi=doi))
    return out


def from_circuits() -> list[Citation]:
    path = REPO_ROOT / "worm" / "data" / "circuits" / "circuits.json"
    data = json.loads(path.read_text(encoding="utf-8"))["circuits"]
    out = []
    for name, spec in sorted(data.items()):
        for part in str(spec.get("source", "")).split(";"):
            part = part.strip()
            if len(part) > 20:
                out.append(Citation(where=f"circuits.json [{name}]", text=part))
    return out


def from_references_md() -> list[Citation]:
    """Paragraphs of references.md that contain a DOI link."""
    text = (REPO_ROOT / "docs" / "references.md").read_text(encoding="utf-8")
    out = []
    for block in re.split(r"\n\s*\n", text):
        stripped = block.lstrip()
        if stripped.startswith("|"):
            continue  # a table is not a citation
        if stripped.startswith(("-", "*")):
            # A bullet list holds one citation per item. Flattening the whole block
            # would produce a soup of unrelated author names and false mismatches.
            chunks = [c for c in re.split(r"\n\s*[-*]\s+", "\n" + block) if c.strip()]
        else:
            chunks = [block]
        for chunk in chunks:
            doi = _DOI_RE.search(chunk)
            if not doi:
                continue
            flat = " ".join(chunk.split())
            out.append(
                Citation(where=f"docs/references.md: {flat[:60]}...", text=flat, doi=doi.group(0))
            )
    return out


def collect() -> list[Citation]:
    return (
        from_sources_toml() + from_references_md() + from_curated_notes() + from_circuits()
    )


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


def cmd_check(args: argparse.Namespace) -> int:
    result = Result()
    excused = 0
    for citation in collect():
        record = by_doi(citation.doi) if citation.doi else search(citation.text)
        if record is None:
            result.unresolvable.append(citation)
            continue
        result.checked += 1
        errata = KNOWN_CROSSREF_ERRATA.get((record.get("DOI") or "").lower(), {})
        for mismatch in compare(citation, record):
            if mismatch.field in errata:
                excused += 1
                continue
            result.mismatches.append(mismatch)

    if args.format == "json":
        print(
            json.dumps(
                {
                    "checked": result.checked,
                    "mismatches": [
                        {
                            "where": m.citation.where,
                            "field": m.field,
                            "ours": m.ours,
                            "crossref": m.crossref,
                        }
                        for m in result.mismatches
                    ],
                    "unresolvable": [c.where for c in result.unresolvable],
                },
                indent=2,
            )
        )
        return 1 if result.mismatches else 0

    print(f"checked {result.checked} citations against CrossRef")
    if excused:
        print(
            f"  ({excused} difference(s) excused as errors in CrossRef's own metadata; "
            "see KNOWN_CROSSREF_ERRATA)"
        )
    print()
    for m in result.mismatches:
        print(m)
        print()
    if result.unresolvable:
        print(f"could not resolve {len(result.unresolvable)}:")
        for c in result.unresolvable:
            print(f"  {c.where}")
        print()
    if result.mismatches:
        print(f"{len(result.mismatches)} MISMATCH(ES). A citation is data; fix it or drop it.")
        return 1
    print("all citations agree with CrossRef.")
    return 0


def cmd_find(args: argparse.Namespace) -> int:
    """Look up citations that carry no DOI, so one can be recorded."""
    for citation in collect():
        if citation.doi:
            continue
        record = search(citation.text)
        if record is None:
            print(f"{citation.where}\n  NO MATCH for: {citation.text[:90]}\n")
            continue
        names = "; ".join(surnames(record))
        print(f"{citation.where}")
        print(f"  ours    : {citation.text[:100]}")
        print(f"  crossref: {names} ({_year(record)}) {record.get('title', ['?'])[0][:70]}")
        print(f"  doi     : {record.get('DOI')}\n")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--format", default="table", choices=["table", "json"])
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("check").set_defaults(func=cmd_check)
    sub.add_parser("find").set_defaults(func=cmd_find)
    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
