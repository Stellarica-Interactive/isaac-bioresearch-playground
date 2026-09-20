"""Recover which channel each published conductance belongs to, from the code.

    python tools/assess_nicoletti2024.py            # from the cached copy
    python tools/assess_nicoletti2024.py --fetch    # download the drivers first

Each cell in Nicoletti et al. 2024 has its maximal conductances in a bare vector
``g0``. A comment above it names the entries, and for three of seven cells that
comment is wrong -- VB6 lists ten conductance names for thirteen values, VD5 lists
ten for eight, AVAL's is misspelled. `docs/model_assumptions.md` 5I concluded from
this that the mapping "cannot be recovered from the repository", and an email was
sent to the authors asking for it.

**The comment is not the mapping.** Each cell's voltage-clamp driver unpacks the
same vector into named NEURON mechanisms, by index:

    seg.slo1egl19.gbar = gVB6_scaled[0]
    seg.slo2egl19.gbar = gVB6_scaled[1]
    seg.slo1unc2.gbar  = gVB6_scaled[2]
    ...

That is unambiguous, it is executable rather than prose, and it is what the
published simulations actually ran. This module reads *that*, so the comments are
only used as a cross-check. All seven cells are recoverable, including VB6 and
VD5 -- the two wanted most, and the two 5I gave up on.

The lesson is narrow and worth keeping: a comment is a claim about code, and the
code is the fact. 5I compared a comment against a vector, found they disagreed,
and concluded the information was absent -- without looking for it anywhere else
in the repository. See 5W.

The repository carries **no licence** (verified via the GitHub API on 2026-09-20:
``license: null``), so nothing from it is redistributed here. Drivers are fetched
to a scratch directory outside version control, parsed for their numbers, and the
numbers are re-derived into our own implementation with the CC-BY paper cited --
the same treatment as the 2019 XPPAUT source.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path

REPO = "martinanicoletti92/CelegansInterMotorNeuronsModels"
RAW = f"https://raw.githubusercontent.com/{REPO}/main/"
TREE = f"https://api.github.com/repos/{REPO}/git/trees/main?recursive=1"

#: Tolerant on purpose. "conductances", "CONDUCTANCES", the misspelled
#: "coductances" (missing the n, not the o -- the first version of this pattern
#: allowed the wrong letter to be optional and still missed AVAL), and an optional
#: "in S/cm^2" all occur in the published files.
COMMENT = re.compile(r"#\s*co?n?ductances?[^:]*:\s*(.+)", re.IGNORECASE)
#: The vector is ``g0`` in five files and ``g`` in RIM's.
VECTOR = re.compile(r"^\s*g0?\s*=\s*\[([^\]]*)\]", re.MULTILINE)
#: ``gScm2(g0, surf, k)`` -- k is the index of the last conductance, so the
#: entries after it are reversal potentials and capacitances, not conductances.
INDEX = re.compile(r"gScm2\s*\(\s*g0?\s*,\s*[^,]+,\s*(\d+)\s*\)")
#: ``seg.<mechanism>.gbar = g<CELL>_scaled[<i>]`` in the clamp drivers. This is
#: the authoritative mapping: it is what ran.
GBAR = re.compile(r"seg\.(\w+)\.gbar\s*=\s*\w+\[(\d+)\]")
#: ``soma.cm = g<CELL>_scaled[<i>]`` -- pins the tail of the vector.
CM = re.compile(r"soma\.cm\s*=\s*\w+\[(\d+)\]")


@dataclass(frozen=True)
class Assessment:
    cell: str
    source: str
    names: tuple[str, ...]
    """Channel names as the comment claims them. A cross-check, never the source."""
    values: tuple[str, ...]
    last_conductance: int | None
    mapping: tuple[tuple[int, str], ...] = ()
    """``(index, mechanism)`` recovered from the clamp driver. The actual mapping."""

    @property
    def recovered(self) -> bool:
        """Every conductance slot is accounted for by the code."""
        if not self.mapping or self.last_conductance is None:
            return False
        return {i for i, _ in self.mapping} == set(range(self.last_conductance + 1))

    @property
    def comment_agrees(self) -> bool:
        """Whether the comment would have given the same answer.

        False for VB6, VD5 and AVAL. Recorded rather than ignored, because the
        comment being wrong is the fact that made 5I give up.
        """
        if not self.recovered:
            return False
        claimed = [n.lower() for n in self.names[: self.last_conductance + 1]]
        actual = [m.lower() for _, m in sorted(self.mapping)]
        return claimed == actual

    @property
    def verdict(self) -> str:
        if not self.mapping:
            return "no clamp driver found"
        if not self.recovered:
            got = len(self.mapping)
            want = (self.last_conductance or 0) + 1
            return f"INCOMPLETE: {got} assignments for {want} conductances"
        note = "" if self.comment_agrees else "  (comment disagrees; code wins)"
        return f"recovered: {len(self.mapping)} conductances{note}"


def assess(source: str, text: str, clamp: str) -> Assessment:
    """Read one cell. ``text`` is the driver holding ``g0``; ``clamp`` is the
    voltage-clamp module that unpacks it into named mechanisms."""
    comment = COMMENT.search(text)
    names = (
        tuple(n.strip() for n in comment.group(1).split(",") if n.strip())
        if comment
        else ()
    )
    vector = VECTOR.search(text)
    values = (
        tuple(v.strip() for v in vector.group(1).split(",") if v.strip())
        if vector
        else ()
    )
    index = INDEX.search(text)
    last = int(index.group(1)) if index else None
    if last is None and values:
        # RIM passes conductances already in S/cm^2, so there is no gScm2 call.
        # The vector still ends with eleak and cm.
        last = len(values) - 3
    mapping = tuple(
        sorted((int(i), mech) for mech, i in GBAR.findall(clamp))
    )
    return Assessment(
        cell=source.split("_")[0],
        source=source,
        names=names,
        values=values,
        last_conductance=last,
        mapping=mapping,
    )


def drivers(cache: Path, fetch: bool) -> dict[str, str]:
    cache.mkdir(parents=True, exist_ok=True)
    if fetch:
        with urllib.request.urlopen(TREE, timeout=60) as response:
            tree = json.load(response)
        names = [
            entry["path"]
            for entry in tree["tree"]
            if entry["type"] == "blob" and entry["path"].endswith(".py")
        ]
        for name in names:
            target = cache / name
            if not target.exists():
                with urllib.request.urlopen(RAW + name, timeout=60) as response:
                    target.write_bytes(response.read())
    return {
        path.name: path.read_text(encoding="utf-8", errors="replace")
        for path in sorted(cache.glob("*.py"))
        if not path.name.startswith("g_to_")
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("--fetch", action="store_true", help="download first")
    parser.add_argument(
        "--cache",
        type=Path,
        default=Path("build") / "nicoletti2024",
        help="where the drivers are kept. Not under version control: the "
        "repository has no licence and is not redistributed.",
    )
    args = parser.parse_args()

    found = drivers(args.cache, args.fetch)
    if not found:
        print(f"no drivers in {args.cache}; run with --fetch", file=sys.stderr)
        return 2

    rows = []
    for name, text in found.items():
        if not (COMMENT.search(text) or VECTOR.search(text)):
            continue
        cell = name.split("_")[0]
        clamp = next(
            (
                body
                for other, body in found.items()
                if other.startswith(cell) and "vclamp" in other
            ),
            "",
        )
        rows.append(assess(name, text, clamp))

    print(f"{'cell':6s} {'values':>6} {'cond':>5} {'mapped':>6}  verdict")
    for row in sorted(rows, key=lambda r: r.cell):
        print(
            f"{row.cell:6s} {len(row.values):6d} "
            f"{(row.last_conductance or 0) + 1:5d} {len(row.mapping):6d}"
            f"  {row.verdict}"
        )

    recovered = sorted(r.cell for r in rows if r.recovered)
    blocked = sorted(r.cell for r in rows if not r.recovered)
    disagree = sorted(r.cell for r in rows if r.recovered and not r.comment_agrees)
    print(f"{chr(10)}recoverable ({len(recovered)}): {', '.join(recovered)}")
    if blocked:
        print(f"blocked     ({len(blocked)}): {', '.join(blocked)}")
    print(f"comment wrong, code right: {', '.join(disagree) or 'none'}")

    for row in sorted(rows, key=lambda r: r.cell):
        if row.recovered:
            order = ", ".join(m for _, m in sorted(row.mapping))
            print(f"{chr(10)}  {row.cell}: {order}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
