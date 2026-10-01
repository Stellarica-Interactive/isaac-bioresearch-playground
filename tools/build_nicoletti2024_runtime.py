"""Assemble a runnable model TOML for a Nicoletti 2024 cell.

    python tools/build_nicoletti2024_runtime.py AVAL --out worm/neural/models/

Takes the channel kinetics from ``worm/neural/models/channels2024/`` and the
per-cell conductances from ``worm/neural/models/cells2024/`` -- both produced by
the importers, both carrying source line numbers -- and writes the single
``[parameters]`` table :func:`common.neural.conductance.load_parameters` expects.

Why the 2024 constants and not ours
-----------------------------------

Thirteen of the fifteen channels share a gating *formula* with the 2019 import,
and it is tempting to reuse the 2019 numbers for them. SHK-1 and SHL-1 were
refitted for the 2024 paper -- ``vashak`` moves from 20.4 mV to 2.0 -- so a cell
built that way would carry the right channels at the wrong densities, which is
§5W.2's failure. Constants therefore come from the 2024 sources, and a channel
whose kinetics this project has not implemented is refused rather than
substituted.

Units
-----

``g_to_Scm2.py`` computes ``g*1e-9/surf``, so the published vector is in
**nanosiemens** and our runtime uses nS directly: the conductances pass through
unconverted, and the surface area is needed only for the capacitance. The tail of
each vector is the leak reversal followed by specific capacitance in uF/cm^2,
which becomes pF as ``cm * surface * 1e6``.
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path

MODELS = Path("worm") / "neural" / "models"
#: Where the drivers cache, for the ion reversal potentials. NMODL declares
#: ``ek`` and ``eca`` through ``USEION`` without values, so they are not in the
#: channel files -- each cell's clamp driver sets them with ``seg.ek = -80``.
DRIVERS = Path("build") / "nicoletti2024"

#: Which of our implemented channel names each 2024 mechanism corresponds to.
#: Only mechanisms whose gating `common/neural/conductance.py` already computes
#: appear here; anything else must be implemented before a cell using it can be
#: built, and is reported rather than guessed at.
MECHANISM_TO_CHANNEL = {
    "egl19": "egl19",
    "irk": "kir",
    "nca": "nca",
    "leak": "leak",
    "cca1": "cca1",
    "unc2": "unc2",
    "egl2": "egl2",
    "shl1": "shal",
    "shk1": "shak",
}

#: Where a 2024 source uses one name for a quantity our implementation spells
#: several ways. Established by comparing the formulas character by character,
#: not by name similarity -- which is the §5W.2 failure.
#:
#: ``egl19``: the 2024 file carries a single ``shift`` that appears in all four
#: of its expressions::
#:
#:     tact:   (pdg1+(pdg2*exp(-(v-pdg3+shift)^2/...)))*ctm19
#:     tinact: pds1*(((pds2*pds3)/(1+exp((v-pds4+shift)/pds5)))+...)
#:
#: against our ``stau19`` in the first and ``shiftdps`` in the second, with
#: ``stm19`` and ``sth19`` in activation and inactivation. All four are 10 in the
#: 2019 import and ``shift`` is 10 here, so the expressions are identical.
ALIASES: dict[str, dict[str, tuple[str, ...]]] = {
    "egl19": {"shift": ("stm19", "sth19", "stau19", "shiftdps")},
    # NCA is a sodium leak; its reversal is the mechanism's own ``e``, which the
    # drivers do not override, so the .mod default of +30 mV stands. Our
    # implementation calls it ``ena``.
    "nca": {"e": ("ena",)},
}

#: Parameters that are per-mechanism and must never become global. Both
#: ``leak.mod`` and ``nca.mod`` name their reversal potential ``e``, so a shared
#: key silently resolves to whichever channel was read first. The leak reversal
#: comes from the cell's own vector instead -- ``seg.leak.e = g[4]`` in the
#: driver -- and NCA's is aliased above.
PER_MECHANISM_ONLY = frozenset({"e"})

#: Channels whose 2024 kinetics our implementation cannot express, with why.
#: These refuse a build rather than being approximated. Empty is the goal.
STRUCTURAL_MISMATCH: dict[str, str] = {}

#: Scale factors that appear as bare literals in a 2024 formula rather than as
#: named parameters, so the parser cannot extract them. Quoted from the source
#: alongside the value, because a divisor written into an equation is as much a
#: constant as one given a name -- and omitting it is not a rounding error.
#:
#: ``shl1.mod``::
#:
#:     mtau  = (ptmshal1/(exp(...)+exp(...))+ptmshal6)/2
#:     htauf = (pthfshal1/(1+exp(...))+pthfshal4)/3
#:     htaus = (pthsshal1/(1+exp(...))+pthsshal4)
#:
#: against the 2019 sources, which multiply activation by ``cashal`` and both
#: inactivation gates by a shared ``cshal``.
EMBEDDED_SCALES: dict[str, dict[str, float]] = {
    "shl1": {"cashal": 1.0 / 2.0, "cthfshal": 1.0 / 3.0, "cthsshal": 1.0},
}

#: Parameters present in a 2024 source whose effect our implementation does not
#: carry, together with the value at which omitting them is exact. A model is
#: refused if one of these is not at its neutral value, rather than quietly
#: dropping a term.
#:
#: ``ctm19`` multiplies EGL-19's activation time constant and is 1.
NEUTRAL_OMISSIONS: dict[str, dict[str, float]] = {
    "egl19": {"ctm19": 1.0},
}

#: Parameter name our implementation uses for each channel's conductance.
CONDUCTANCE_KEY = {
    "egl19": "gegl19",
    "kir": "gkir",
    "nca": "gnca",
    "leak": "gleak",
    "cca1": "gcca1",
    "unc2": "gunc2",
    "egl2": "gegl2",
    "shal": "gshal",
    "shak": "gshak",
}


def _load(path: Path) -> dict:
    return tomllib.loads(path.read_text(encoding="utf-8"))


def _reversals(cell: str) -> dict[str, float]:
    """``ek`` and ``eca`` as the cell's driver sets them.

    NMODL declares these through ``USEION ... READ ek`` with no value, so they
    are absent from the channel files and cannot be taken from there. Every
    driver sets them explicitly -- ``seg.ek = -80``, ``seg.eca = 60`` -- and that
    is the only place they exist.
    """
    found: dict[str, float] = {}
    for path in sorted(DRIVERS.glob(f"{cell.upper()}_*.py")):
        text = path.read_text(encoding="utf-8", errors="replace")
        for ion in ("ek", "eca"):
            match = re.search(rf"\.{ion}\s*=\s*(-?[\d.eE+]+)", text)
            if match:
                found.setdefault(ion, float(match.group(1)))
    missing = {"ek", "eca"} - set(found)
    if missing:
        raise KeyError(
            f"{cell}: no driver in {DRIVERS} sets {sorted(missing)}. These are "
            "ion reversal potentials and must come from the source, not a default."
        )
    return found


def build(cell: str) -> tuple[dict[str, float], tuple[str, ...], list[str]]:
    """Parameters, channel tuple, and any mechanisms we cannot yet model."""
    cell_data = _load(MODELS / "cells2024" / f"{cell.lower()}_nicoletti2024.toml")
    conductances = cell_data["conductances"]
    surface = float(cell_data["cell"]["surface_cm2"])
    tail = list(cell_data["tail"]["values"])

    unsupported = [m for m in conductances if m not in MECHANISM_TO_CHANNEL]
    channels = tuple(MECHANISM_TO_CHANNEL[m] for m in conductances if m in MECHANISM_TO_CHANNEL)

    parameters: dict[str, float] = {}
    for mechanism, value in conductances.items():
        channel = MECHANISM_TO_CHANNEL.get(mechanism)
        if channel is None:
            continue
        parameters[CONDUCTANCE_KEY[channel]] = float(value)
        source = MODELS / "channels2024" / f"{mechanism}_nicoletti2024.toml"
        if not source.exists():
            continue
        published = _load(source)["parameters"]
        aliases = ALIASES.get(mechanism, {})
        for key, expected in NEUTRAL_OMISSIONS.get(mechanism, {}).items():
            actual = published.get(key)
            if actual is not None and float(actual) != expected:
                unsupported.append(
                    f"{mechanism}.{key}={actual} (we omit this term, which is only "
                    f"exact at {expected})"
                )
        if mechanism in STRUCTURAL_MISMATCH:
            unsupported.append(f"{mechanism}: {STRUCTURAL_MISMATCH[mechanism]}")
        for key, value in EMBEDDED_SCALES.get(mechanism, {}).items():
            parameters.setdefault(key, value)
        for key, number in published.items():
            if key == "gbar":  # per-cell, taken from the cell vector above
                continue
            names = aliases.get(key)
            if names is None:
                if key in PER_MECHANISM_ONLY:
                    continue
                names = (key,)
            for name in names:
                parameters.setdefault(name, float(number))

    # Ion reversal potentials, from the cell's own driver rather than assumed.
    for name, value in _reversals(cell).items():
        parameters[name] = value

    # Leak reversal, then specific capacitance. The source does not name them, so
    # position is all there is -- which is why the importer records the tail
    # verbatim instead of guessing labels for it.
    parameters["eleak"] = float(tail[0])
    # `c`, not `c_m`: that is the name our implementation integrates with.
    parameters["c"] = float(tail[-1]) * surface * 1e6
    return parameters, channels, unsupported


def to_toml(cell: str, parameters: dict[str, float], channels: tuple[str, ...]) -> str:
    out = [
        f"# {cell}, Nicoletti et al. 2024. GENERATED by",
        "# tools/build_nicoletti2024_runtime.py from channels2024/ and cells2024/.",
        "# Regenerate rather than editing.",
        "#",
        "# Channel kinetics come from the 2024 .mod sources, NOT from our 2019",
        "# import: SHK-1 and SHL-1 were refitted for this paper, so reusing 2019",
        "# constants would give the right channels at the wrong densities.",
        "# See docs/model_assumptions.md 5W.5 and 5X.",
        "",
        "[source]",
        'source_id = "nicoletti_2024_motor_interneurons"',
        f'cell = "{cell}"',
        'url = "https://github.com/martinanicoletti92/CelegansInterMotorNeuronsModels"',
        'citation = "Nicoletti M, Chiodo L, Loppini A, Liu Q, Folli V, Ruocco G, '
        "Filippi S. Biophysical modeling of the whole-cell dynamics of C. elegans "
        'motor and interneurons families. PLoS ONE 19(3):e0298105 (2024)."',
        'doi = "10.1371/journal.pone.0298105"',
        'license = "Repository carries no licence (verified 2026-09-20). Not '
        'redistributed; constants re-derived and cited. Paper is CC-BY 4.0."',
        f"channels = [{', '.join(repr(c) for c in channels)}]".replace("'", '"'),
        "",
        "[parameters]",
    ]
    for key, value in sorted(parameters.items()):
        out.append(f"{key} = {value!r}")
    return "\n".join(out) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    parser.add_argument("cells", nargs="+", help="e.g. AVAL AVAR")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    failed = 0
    for cell in args.cells:
        parameters, channels, unsupported = build(cell)
        if unsupported:
            print(
                f"{cell}: needs channels this project has not implemented: "
                f"{', '.join(sorted(unsupported))} -- refusing to write a model "
                f"with them silently omitted",
                file=sys.stderr,
            )
            failed += 1
            continue
        print(
            f"{cell:5s} {len(channels)} channels ({', '.join(channels)}), "
            f"{len(parameters)} parameters, "
            f"c = {parameters['c']:.2f} pF, eleak = {parameters['eleak']:g} mV"
        )
        if args.out:
            args.out.mkdir(parents=True, exist_ok=True)
            path = args.out / f"{cell.lower()}_nicoletti2024.toml"
            path.write_text(to_toml(cell, parameters, channels), encoding="utf-8", newline="\n")
            print(f"      -> {path}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
