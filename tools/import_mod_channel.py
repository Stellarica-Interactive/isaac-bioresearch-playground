"""Parse NEURON ``.mod`` channel files into our own TOML, with line numbers.

    python tools/import_mod_channel.py kqt1 unc103 --fetch --out worm/neural/models/

A companion to :mod:`tools.import_neuron_model`, which does the same job for the
XPPAUT ``.ode`` sources of Nicoletti et al. 2019. This one reads the NMODL sources
of Nicoletti et al. 2024.

Why a parser rather than retyping
---------------------------------

Same reason as the ``.ode`` importer, and it has now been vindicated twice. Every
constant carries the file and line it came from, so any value in a generated TOML
can be checked against the published source in seconds, and nothing enters the
model because somebody believed they remembered it.

The 2024 sources make a second argument for it. Thirteen of the fifteen channels
share a gating *formula* with our 2019 import, but SHK-1 and SHL-1 carry different
fitted *constants* -- ``vashak`` is 20.4 mV in 2019 and 2.0 in 2024 (see
docs/model_assumptions.md 5W.5). A model assembled by matching channel names
against what we already hold would be silently wrong. Parsing every constant out
of the 2024 file makes that impossible rather than merely unlikely.

What is extracted and what is not
---------------------------------

Extracted: the ``PARAMETER`` block's named constants, the ``STATE`` variables, the
``BREAKPOINT`` current expression, and each ``FUNCTION``'s returned expression,
verbatim as text.

**Not** extracted: any translation of those expressions into executable Python.
Formulas are recorded as strings, exactly as ``import_neuron_model`` does, and
turning one into a gating function stays a deliberate, reviewed act in
``common/neural/conductance.py``. A parser that also generated the maths would be
the one place an error could enter without anyone reading it.

Units
-----

NMODL declares units per quantity and the 2024 files are inconsistent about them:
``gbar`` is commented ``(nS/cm2)`` in ``irk.mod`` while the drivers pass S/cm^2.
Units are captured as written and **not** converted here. Conversion belongs with
the per-cell conductances, where the surface area that makes it meaningful lives.
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

REPO = "martinanicoletti92/CelegansInterMotorNeuronsModels"
RAW = f"https://raw.githubusercontent.com/{REPO}/main/"

ASSIGN = re.compile(r"^\s*([A-Za-z_]\w*)\s*=\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*(\([^)]*\))?")
# The signature may contain nested parentheses -- ``FUNCTION minf(v (mV))`` -- so
# the argument list is matched up to the opening brace rather than to the first
# close paren. The stricter version found no functions at all, and produced a TOML
# carrying every constant and no kinetics.
FUNCTION = re.compile(r"FUNCTION\s+(\w+)\s*\([^{]*\{(.*?)\n\s*\}", re.S)
# The SLO ``iso`` variants put their kinetics in ``PROCEDURE rates(v, ca)`` rather
# than in one FUNCTION per quantity, so every assignment in the body is a formula.
PROCEDURE = re.compile(r"PROCEDURE\s+\w+\s*\([^{]*\{(.*?)\n\s*\}", re.S)
ASSIGNMENT = re.compile(r"^\s*([A-Za-z_]\w*)\s*=\s*(.+?)\s*$", re.M)
SUFFIX = re.compile(r"SUFFIX\s+(\w+)")
USEION = re.compile(r"USEION\s+(\w+)")


@dataclass
class Channel:
    """One ``.mod`` file, read but not interpreted."""

    name: str
    source: str
    ion: str | None = None
    parameters: dict[str, tuple[float, str | None, int]] = field(default_factory=dict)
    """``name -> (value, unit as written, 1-based line in the .mod file)``."""
    states: tuple[str, ...] = ()
    current: str | None = None
    formulas: dict[str, str] = field(default_factory=dict)
    procedure: list[tuple[str, str]] = field(default_factory=list)
    """Assignments from a ``PROCEDURE`` body, **in source order**.

    The SLO ``iso`` variants compute their kinetics this way, and later lines
    depend on names bound by earlier ones, so the order is part of the content.
    """

    @property
    def is_passive(self) -> bool:
        """No gating: the current depends on voltage alone.

        ``leak`` and ``nca`` are exactly ``gbar*(v - e)`` with no STATE block.
        Distinguishing these from a parse failure matters -- refusing to write
        them would be as wrong as writing a gated channel with no kinetics.
        """
        return not self.states

    @property
    def kinetics_missing(self) -> bool:
        """Has state variables but nothing that says how they evolve.

        Always a parser failure, never a real channel. The first version of this
        tool checked only for absent formulas, which refused ``leak`` and ``nca``
        for being passive and accepted nothing in their place.
        """
        return bool(self.states) and not (self.formulas or self.procedure)


def _block(text: str, name: str) -> str:
    """The body of a top-level NMODL block, or an empty string.

    Blocks may be written on one line. ``leak.mod`` and ``nca.mod`` are entirely
    ``BREAKPOINT { i = gbar*(v - e) }``, and an earlier pattern that required a
    newline before the closing brace silently returned nothing for them -- so the
    two passive channels imported with their constants and without the single
    equation that is their whole content.
    """
    match = re.search(rf"(?:^|\n)\s*{name}\s*\{{(.*?)\}}", text, re.S)
    return match.group(1) if match else ""


def parse(name: str, text: str) -> Channel:
    lines = text.splitlines()
    channel = Channel(name=name, source=f"{name}.mod")

    suffix = SUFFIX.search(text)
    if suffix:
        channel.name = suffix.group(1)
    ion = USEION.search(text)
    channel.ion = ion.group(1) if ion else None

    for raw in _block(text, "PARAMETER").splitlines():
        found = ASSIGN.match(raw)
        if not found:
            continue
        key, value, unit = found.groups()
        line = next((i + 1 for i, original in enumerate(lines) if original == raw), 0)
        channel.parameters[key] = (float(value), unit, line)

    channel.states = tuple(re.findall(r"[A-Za-z_]\w*", _block(text, "STATE")))

    current = re.search(r"^\s*i\w*\s*=\s*(.+)$", _block(text, "BREAKPOINT"), re.M)
    channel.current = current.group(1).strip() if current else None

    for function, body in FUNCTION.findall(text):
        expression = re.search(rf"^\s*{function}\s*=\s*(.+)$", body, re.M)
        if expression:
            channel.formulas[function] = expression.group(1).strip()

    # PROCEDURE bodies compute several quantities in sequence, and the order
    # matters: later lines use names bound by earlier ones. Recorded in source
    # order rather than sorted, because sorting them would lose that.
    for body in PROCEDURE.findall(text):
        for name_, expression in ASSIGNMENT.findall(body):
            if name_ in ("UNITSOFF", "UNITSON"):
                continue
            channel.procedure.append((name_, expression))

    return channel


def _number(value: float) -> str:
    text = repr(value)
    return text[:-2] if text.endswith(".0") else text


def to_toml(channel: Channel) -> str:
    quoted_states = ", ".join(f'"{s}"' for s in channel.states)
    out = [
        f"# Generated by tools/import_mod_channel.py from {channel.source}.",
        "# Nicoletti M, Chiodo L, Loppini A, Liu Q, Folli V, Ruocco G, Filippi S.",
        "# Biophysical modeling of the whole-cell dynamics of C. elegans motor and",
        "# interneurons families. PLoS ONE 19(3):e0298105 (2024).",
        "#",
        "# Line numbers refer to the published .mod file, so every constant can be",
        "# checked against the source. Regenerate rather than editing.",
        "",
        "[channel]",
        f'name = "{channel.name}"',
        f'source = "{channel.source}"',
        f'ion = "{channel.ion or ""}"',
        f"states = [{quoted_states}]",
    ]
    if channel.current:
        out.append(f'current = "{channel.current}"')
    out += ["", "[parameters]"]
    for key, (value, unit, line) in sorted(channel.parameters.items()):
        note = f"  # line {line}" + (f", {unit}" if unit else "")
        out.append(f"{key} = {_number(value)}{note}")
    out += ["", "[formulas]"]
    for key, expression in sorted(channel.formulas.items()):
        out.append(f'{key} = "{expression}"')
    if channel.procedure:
        out += [
            "",
            "# From a PROCEDURE body. `order` is the sequence the source evaluates",
            "# in, and it is load-bearing: later expressions use names bound by",
            "# earlier ones, so evaluating these in any other order is wrong.",
            "[procedure]",
            "order = [" + ", ".join(f'"{k}"' for k, _ in channel.procedure) + "]",
            "",
            "[procedure.expressions]",
        ]
        for key, expression in channel.procedure:
            out.append(f'{key} = "{expression}"')
    return "\n".join(out) + "\n"


def load(name: str, cache: Path, download: bool) -> str:
    target = cache / f"{name}.mod"
    if download and not target.exists():
        cache.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(RAW + f"{name}.mod", timeout=60) as response:
            target.write_bytes(response.read())
    return target.read_text(encoding="utf-8", errors="replace")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("channels", nargs="*", help="e.g. kqt1 unc103")
    parser.add_argument("--fetch", action="store_true", help="download if absent")
    parser.add_argument(
        "--cache",
        type=Path,
        default=Path("build") / "nicoletti2024",
        help="where .mod sources are kept. Not under version control: the "
        "repository carries no licence and is not redistributed.",
    )
    parser.add_argument("--out", type=Path, default=None, help="write TOML here")
    args = parser.parse_args()

    names = args.channels or [p.stem for p in sorted(args.cache.glob("*.mod"))]
    if not names:
        print("no channels named and none cached; use --fetch", file=sys.stderr)
        return 2

    failures = 0
    for name in names:
        try:
            text = load(name, args.cache, args.fetch)
        except (OSError, urllib.error.URLError) as error:
            print(f"{name}: {error}", file=sys.stderr)
            failures += 1
            continue
        channel = parse(name, text)
        if channel.kinetics_missing:
            # A channel with constants and no kinetics is worse than no channel:
            # it imports cleanly, carries authoritative-looking line numbers, and
            # models nothing. Refuse rather than write it.
            print(
                f"{name}: has states {channel.states} but no FUNCTION or PROCEDURE "
                f"body was parsed -- refusing to write kinetics-free channel",
                file=sys.stderr,
            )
            failures += 1
            continue
        rendered = to_toml(channel)
        if args.out:
            args.out.mkdir(parents=True, exist_ok=True)
            path = args.out / f"{channel.name}_nicoletti2024.toml"
            path.write_text(rendered, encoding="utf-8", newline="\n")
            print(
                f"{channel.name:12s} {len(channel.parameters):3d} parameters, "
                f"{len(channel.states)} states, {len(channel.formulas)} formulas"
                f"  -> {path}"
            )
        else:
            print(rendered)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
