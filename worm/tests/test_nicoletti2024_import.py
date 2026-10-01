"""The imported Nicoletti 2024 channel and cell data.

These guard the *data*, not a simulation: the models are not runnable yet
(docs/model_assumptions.md 5X.3). That is exactly when this kind of test earns
its place, because an import with nobody consuming it has nothing else checking
whether it is coherent, and the failures in 5X.2 were all files that looked
correct and were missing their content.
"""

from __future__ import annotations

import math
import tomllib
from pathlib import Path

import pytest

MODELS = Path(__file__).resolve().parents[2] / "worm" / "neural" / "models"
CHANNELS = MODELS / "channels2024"
CELLS = MODELS / "cells2024"

#: The seven cells of the published model.
EXPECTED_CELLS = {"aiy", "aval", "avar", "rim", "va5", "vb6", "vd5"}

#: Channels with no gating: their current depends on voltage alone. Named so the
#: "has states but no kinetics" check cannot be satisfied by simply declaring
#: every channel passive.
PASSIVE = {"leak", "nca"}


def _load(path: Path) -> dict:
    return tomllib.loads(path.read_text(encoding="utf-8"))


def _channels() -> dict[str, dict]:
    return {_load(p)["channel"]["name"]: _load(p) for p in CHANNELS.glob("*.toml")}


def _cells() -> dict[str, dict]:
    return {p.stem.split("_")[0]: _load(p) for p in CELLS.glob("*.toml")}


def test_all_seven_cells_imported() -> None:
    """5I concluded four of seven were unrecoverable; 5W showed they are not."""
    assert set(_cells()) == EXPECTED_CELLS


def test_vb6_and_vd5_are_present() -> None:
    """Called out separately because they are the two cells 5I gave up on, and
    the two an unanswered email was sent to ask about."""
    cells = _cells()
    assert len(cells["vb6"]["conductances"]) == 13
    assert len(cells["vd5"]["conductances"]) == 8


def test_every_channel_a_cell_uses_was_imported() -> None:
    """A cell naming a channel we did not import is an unrunnable model, and the
    gap would otherwise only surface when somebody tried to build it."""
    available = set(_channels())
    for name, cell in _cells().items():
        missing = set(cell["conductances"]) - available
        assert not missing, f"{name} needs channels we have not imported: {missing}"


@pytest.mark.parametrize("channel", sorted(_channels()))
def test_a_gated_channel_has_kinetics(channel: str) -> None:
    """The failure mode of 5X.2: a file with every constant and no equations.

    Four separate parser bugs produced exactly this, each yielding a well-formed
    TOML with an authoritative provenance header.
    """
    data = _channels()[channel]
    kinetics = dict(data.get("formulas", {}))
    kinetics.update(data.get("procedure", {}).get("expressions", {}))
    # A calcium pool keeps its update in BREAKPOINT, including a conditional
    # that cannot be expressed as a named formula, so it is recorded verbatim
    # and implemented by hand. Counting it here keeps the check meaningful
    # without pretending the text has been parsed into expressions.
    if data.get("breakpoint", {}).get("source"):
        kinetics["breakpoint"] = data["breakpoint"]["source"]
    if channel in PASSIVE:
        assert not data["channel"]["states"]
        assert data["channel"]["current"], "a passive channel is still a current"
    else:
        assert data["channel"]["states"], f"{channel} declares no state"
        assert kinetics, f"{channel} has states but no kinetics"


@pytest.mark.parametrize("channel", sorted(_channels()))
def test_state_names_are_not_nmodl_keywords(channel: str) -> None:
    """``STATE { m FROM 0 TO 1 }`` is one state with bounds, not three.

    slo1iso and slo2iso were imported with states ('m', 'FROM', 'TO'), which
    would have allocated and integrated two phantom gating variables per cell.
    """
    states = _channels()[channel]["channel"]["states"]
    assert not {"FROM", "TO", "START"} & set(states)
    assert len(states) == len(set(states)), "a state declared twice"


@pytest.mark.parametrize("channel", sorted(_channels()))
def test_every_constant_is_finite(channel: str) -> None:
    for key, value in _channels()[channel]["parameters"].items():
        assert math.isfinite(value), f"{channel}.{key} is not finite"


def test_the_calcium_pool_update_is_captured_whole() -> None:
    """``cadiff.mod`` integrates calcium in BREAKPOINT behind a conditional.

    Read up to the first closing brace, the body truncated mid-conditional and
    lost both the 100 nM floor and the ``cai = ca`` that publishes the result --
    text that looked complete and was not. The block is brace-matched now.
    """
    channels = _channels()
    if "cadiff" not in channels:
        pytest.skip("cadiff not imported")
    source = channels["cadiff"]["breakpoint"]["source"]
    assert "ca = ca + (10000) * dt" in source, "the integration step"
    assert "if ( ca < 1e-4 )" in source, "the floor"
    assert "cai = ca" in source, "the line that publishes the result"
    assert source.count("{") == source.count("}"), "braces must balance"


def test_procedure_order_covers_its_expressions() -> None:
    """Order is content for the SLO iso variants: later expressions use names
    bound by earlier ones, so a missing entry is a broken evaluation sequence."""
    for name, data in _channels().items():
        procedure = data.get("procedure")
        if not procedure:
            continue
        assert set(procedure["order"]) == set(procedure["expressions"]), name
        assert len(procedure["order"]) == len(procedure["expressions"]), name


@pytest.mark.parametrize("cell", sorted(EXPECTED_CELLS))
def test_conductances_are_plausible_in_nanosiemens(cell: str) -> None:
    """The published vectors are in nS, despite ``irk.mod`` commenting ``gbar``
    as ``(nS/cm2)``: ``g_to_Scm2.py`` computes ``g*1e-9/surf``, so the value
    entering it is siemens-scaled-by-1e-9. Our runtime works in nS, so these are
    usable unconverted -- and a unit error would show up here as values orders of
    magnitude away from the 0.05-3.1 nS our 2019 import holds.
    """
    conductances = _cells()[cell]["conductances"]
    assert conductances, f"{cell} has no conductances"
    for name, value in conductances.items():
        assert 0.0 <= value <= 20.0, f"{cell}.{name} = {value} nS is off-scale"
    assert any(v > 0 for v in conductances.values()), "every channel at zero"


@pytest.mark.parametrize("cell", sorted(EXPECTED_CELLS))
def test_surface_area_recorded(cell: str) -> None:
    """Not needed to convert -- see above -- but recorded so that anyone who does
    need NEURON's per-area convention has it without re-reading the source."""
    area = _cells()[cell]["cell"]["surface_cm2"]
    assert 1e-8 < area < 1e-3, f"{cell} surface {area} cm2 is implausible"


def test_comment_disagreement_is_recorded_per_cell() -> None:
    """AIY, VB6 and VD5 have comments that contradict their own code.

    AIY is the dangerous one: its comment has the right COUNT and names ``irk``
    where the code assigns ``shl1``, so a count check passes and the model is
    wrong. 5I called it one of the safe cells on exactly that basis.
    """
    cells = _cells()
    disagree = {n for n, c in cells.items() if not c["cell"]["comment_agrees"]}
    assert disagree == {"aiy", "vb6", "vd5"}
    assert "shl1" in cells["aiy"]["conductances"]
    assert "irk" not in cells["aiy"]["conductances"]
