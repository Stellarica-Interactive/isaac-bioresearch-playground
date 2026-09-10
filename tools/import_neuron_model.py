"""Parse an XPPAUT ``.ode`` model into our own parameter table.

    python tools/import_neuron_model.py --source nicoletti_2019_rmd_awc \\
        --out worm/neural/models/rmd_nicoletti2019.toml

Why a parser rather than typing the numbers in
----------------------------------------------

The RMD model has over a hundred numeric parameters across ten ion channels.
Transcribing those by hand is precisely the failure this repository has already
been bitten by twice -- a fabricated author list, and a DOI written from memory
that turned out to belong to a paper about lake bacteria. A wrong author name is
embarrassing and findable. A wrong Boltzmann midpoint is neither: it produces a
model that runs, looks plausible, and is quietly not the published one.

So the numbers are machine-extracted from the authors' own source, every one
carrying the line it came from, and the committed output is a reviewable diff.

**The code is deliberately the source of truth, not the paper.** Nicoletti et al.
2019 carries a published Correction (doi:10.1371/journal.pone.0256930): a plus
sign is missing from a Boltzmann denominator in the Methods, and twelve
supplementary equations are wrong as printed. The XPPAUT file is the code the
authors actually ran and does not contain those errors. Anyone transcribing from
the original PDF would import the mistakes.

What this does and does not do
------------------------------

It extracts **parameters, initial values, expressions and ODEs** verbatim, as
text, with line numbers. It does **not** evaluate them, and it is not a
translator: the channel equations are implemented by hand in
:mod:`common.neural.conductance` so that they are readable and reviewable, while
the error-prone bulk -- the constants -- is never retyped.
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from dataclasses import dataclass
from datetime import UTC, datetime
from hashlib import sha256
from importlib.resources import files
from pathlib import Path

#: ``par a=1,b=2 c=3`` -- XPPAUT accepts commas and whitespace interchangeably as
#: separators, and the Nicoletti files use both, sometimes on the same line.
_ASSIGNMENT = re.compile(r"([A-Za-z_][A-Za-z_0-9]*)\s*=\s*([^,\s]+)")


@dataclass(frozen=True)
class Entry:
    """One extracted definition, and where in the source file it came from."""

    name: str
    value: str
    line: int


@dataclass(frozen=True)
class OdeModel:
    """A parsed XPPAUT model."""

    source_file: str
    sha256: str
    parameters: tuple[Entry, ...]
    initial: tuple[Entry, ...]
    expressions: tuple[Entry, ...]
    odes: tuple[Entry, ...]
    settings: tuple[str, ...]

    def summary(self) -> str:
        return (
            f"{self.source_file}: {len(self.parameters)} parameters, "
            f"{len(self.odes)} state variables, {len(self.expressions)} expressions"
        )


def parse_ode(text: str, *, source_file: str, digest: str) -> OdeModel:
    """Extract every definition from XPPAUT source.

    Only the active model is read: XPPAUT stops at ``done``, and the Nicoletti
    files keep an alternative voltage-clamp protocol commented out below that
    marker. Reading past it would silently mix two simulation setups.
    """
    parameters: list[Entry] = []
    initial: list[Entry] = []
    expressions: list[Entry] = []
    odes: list[Entry] = []
    settings: list[str] = []

    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.lower() == "done":
            break
        if line.startswith("@"):
            settings.append(line.lstrip("@ ").strip())
            continue

        lowered = line.lower()
        if lowered.startswith(("par ", "param ", "p ")):
            body = line.split(None, 1)[1]
            parameters += [
                Entry(m.group(1), m.group(2), number) for m in _ASSIGNMENT.finditer(body)
            ]
            continue
        if lowered.startswith("init "):
            body = line.split(None, 1)[1]
            initial += [Entry(m.group(1), m.group(2), number) for m in _ASSIGNMENT.finditer(body)]
            continue
        if lowered.startswith("aux "):
            # Auxiliary quantities are XPPAUT plotting outputs, not state. Skipping
            # them is not a loss: each one duplicates an expression already parsed.
            continue

        # `name'=expr` is a differential equation; plain `name=expr` an expression.
        if "=" not in line:
            continue
        lhs, rhs = line.split("=", 1)
        lhs, rhs = lhs.strip(), rhs.strip()
        if lhs.endswith("'"):
            odes.append(Entry(lhs[:-1], rhs, number))
        elif re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", lhs):
            expressions.append(Entry(lhs, rhs, number))

    # XPPAUT lets a constant be written either as `par t_sk=6.3` or as a plain
    # assignment `t_sk=6.3`, and the Nicoletti files use both. Whichever form the
    # author happened to pick, it is a numeric constant, so it belongs with the
    # other constants rather than being unreachable in the expression table.
    literal = [e for e in expressions if _is_number(e.value)]
    parameters += literal
    expressions = [e for e in expressions if not _is_number(e.value)]

    return OdeModel(
        source_file=source_file,
        sha256=digest,
        parameters=tuple(parameters),
        initial=tuple(initial),
        expressions=tuple(expressions),
        odes=tuple(odes),
        settings=tuple(settings),
    )


def _is_number(text: str) -> bool:
    try:
        float(text)
    except ValueError:
        return False
    return True


def _source_entry(source_id: str) -> dict[str, object]:
    manifest = tomllib.loads((files("worm.data") / "sources.toml").read_text(encoding="utf-8"))
    if source_id not in manifest:
        raise SystemExit(f"no such source id in worm/data/sources.toml: {source_id!r}")
    entry = manifest[source_id]
    assert isinstance(entry, dict)
    return entry


def _toml_escape(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def render(model: OdeModel, entry: dict[str, object]) -> str:
    """Emit a deterministic TOML table. Sorted, so a re-import diffs cleanly."""
    lines = [
        "# GENERATED by tools/import_neuron_model.py -- do not edit by hand.",
        "#",
        "# Extracted verbatim from the authors' own XPPAUT source. Every entry",
        "# records the line it came from, so any value can be checked against the",
        "# original in seconds. Regenerate rather than editing.",
        "#",
        f"# {entry.get('citation', '')}",
        f"# doi: {entry.get('doi', '')}",
    ]
    if entry.get("correction_doi"):
        lines += [
            f"# CORRECTION: {entry['correction_doi']} -- the paper's printed equations",
            "# contain errors that this source file does not. The code is authoritative.",
        ]
    lines += [
        "",
        "[source]",
        f'source_id = "{_toml_escape(str(entry.get("source_id", "")))}"',
        f'file = "{_toml_escape(model.source_file)}"',
        f'sha256 = "{model.sha256}"',
        f'url = "{_toml_escape(str(entry.get("url", "")))}"',
        f'citation = "{_toml_escape(str(entry.get("citation", "")))}"',
        f'doi = "{_toml_escape(str(entry.get("doi", "")))}"',
        f'correction_doi = "{_toml_escape(str(entry.get("correction_doi", "")))}"',
        f'license = "{_toml_escape(str(entry.get("license", "")))}"',
        f'imported = "{datetime.now(UTC).date().isoformat()}"',
        'confidence = "published_model"',
        "",
        "# Numeric constants. These are the values that must never be retyped.",
        "[parameters]",
    ]
    for e in sorted(model.parameters, key=lambda x: x.name.lower()):
        lines.append(f"{e.name} = {e.value}  # line {e.line}")

    lines += ["", "[initial]"]
    for e in sorted(model.initial, key=lambda x: x.name.lower()):
        lines.append(f"{e.name} = {e.value}  # line {e.line}")

    lines += [
        "",
        "# Kept as text for cross-checking the hand-written implementation against",
        "# the source. Nothing evaluates these.",
        "[expressions]",
    ]
    for e in sorted(model.expressions, key=lambda x: x.name.lower()):
        lines.append(f'{e.name} = "{_toml_escape(e.value)}"  # line {e.line}')

    lines += ["", "[odes]"]
    for e in sorted(model.odes, key=lambda x: x.name.lower()):
        lines.append(f'{e.name} = "{_toml_escape(e.value)}"  # line {e.line}')

    lines += ["", "[xpp_settings]"]
    for i, s in enumerate(model.settings):
        lines.append(f'setting_{i} = "{_toml_escape(s)}"')

    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--source", required=True, help="source id in worm/data/sources.toml")
    parser.add_argument("--raw-dir", default="worm/data/raw", help="where the fetched file is")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--check",
        action="store_true",
        help="parse and report without writing, for CI",
    )
    args = parser.parse_args(argv)

    entry = dict(_source_entry(args.source))
    entry["source_id"] = args.source
    filename = str(entry.get("filename") or "")
    if not filename:
        raise SystemExit(f"{args.source} has no filename in sources.toml")

    path = Path(args.raw_dir) / filename
    if not path.exists():
        raise SystemExit(
            f"{path} not found. Fetch it first:\n"
            f"    python tools/fetch_datasets.py --dataset {args.source}"
        )

    data = path.read_bytes()
    digest = sha256(data).hexdigest()
    expected = str(entry.get("sha256") or "")
    if expected and digest != expected:
        raise SystemExit(
            f"SHA-256 mismatch for {path}\n  expected {expected}\n  got      {digest}\n"
            "The upstream file changed. Review it before regenerating anything."
        )

    model = parse_ode(data.decode("utf-8", errors="strict"), source_file=filename, digest=digest)
    print(model.summary())
    if not model.odes:
        raise SystemExit("no state variables found; this does not look like an XPPAUT model")

    text = render(model, entry)
    if args.check:
        current = args.out.read_text(encoding="utf-8") if args.out.exists() else ""
        # The import date changes every run and is not part of the model.
        strip = lambda s: "\n".join(  # noqa: E731
            ln for ln in s.splitlines() if not ln.startswith("imported = ")
        )
        if strip(current) != strip(text):
            raise SystemExit(f"{args.out} is stale; re-run without --check")
        print(f"{args.out} is up to date")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text, encoding="utf-8", newline="\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
