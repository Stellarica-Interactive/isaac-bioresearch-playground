"""Parsing NEURON .mod channel sources.

Every test here is a bug the parser actually had. Each one produced output that
looked fine: a TOML with the right channel name, authoritative-looking line
numbers, plausible constants -- and no kinetics, or the wrong ones. That is the
failure mode docs/model_assumptions.md 5W is about, so the tests assert on what
was silently absent rather than on what was present.
"""

from __future__ import annotations

import pytest

from tools.import_mod_channel import parse, to_toml

NESTED_SIGNATURE = """
TITLE test
NEURON {
    SUFFIX demo
    USEION k READ ek WRITE ik
    RANGE gbar
}
PARAMETER {
    gbar = 2.9 (S/cm2)
    va = -17.6053 (mV)
    ka = 9.5843 (mV)
}
STATE { m }
BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gbar*m*(v-ek)
}
FUNCTION minf(v (mV)) {
    UNITSOFF
    minf=1/(1+exp(-(v-va)/ka))
    UNITSON
}
"""

PASSIVE = """
NEURON { SUFFIX leak }
PARAMETER {
    gbar = 1 (S/cm2)
    e = -80 (mV)
}
BREAKPOINT { i = gbar*(v - e) }
"""

PROCEDURE_KINETICS = """
NEURON { SUFFIX slodemo }
PARAMETER { gbar = 0.11 (S/cm2) }
STATE { m }
BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gbar * m * (v - ek)
}
PROCEDURE rates(v (mV), ca (mM)) {
    UNITSOFF
    s0=1/(wyx-wxy)
    v0=s0*(log(wom/wop))
    minf =1/(1+exp(-(v-v0)/s0))
    UNITSON
}
"""

GATED_BUT_SILENT = """
NEURON { SUFFIX broken }
PARAMETER { gbar = 1 (S/cm2) }
STATE { m }
BREAKPOINT { ik = gbar*m*(v-ek) }
"""


def test_a_nested_parenthesis_in_the_signature_does_not_hide_the_function() -> None:
    """``FUNCTION minf(v (mV))`` -- the unit annotation is itself parenthesised.

    Matching the argument list as "up to the first close paren" ends inside
    ``(mV)``, so the pattern never reaches the brace and finds no functions at
    all. Every channel then imported with a complete PARAMETER block and an empty
    [formulas] section: a file that parses, validates, and models nothing.
    """
    channel = parse("demo", NESTED_SIGNATURE)
    assert channel.formulas == {"minf": "1/(1+exp(-(v-va)/ka))"}
    assert not channel.kinetics_missing


def test_a_passive_channel_is_not_mistaken_for_a_parse_failure() -> None:
    """leak and nca have no STATE and no FUNCTION, legitimately.

    The first version of the guard refused anything without formulas, which threw
    away the two channels that are simply ohmic.
    """
    channel = parse("leak", PASSIVE)
    assert channel.is_passive
    assert not channel.kinetics_missing
    assert channel.current == "gbar*(v - e)"
    assert channel.parameters["e"][0] == -80.0


def test_a_gated_channel_with_no_kinetics_is_refused() -> None:
    """States but nothing saying how they evolve is always a parser failure."""
    channel = parse("broken", GATED_BUT_SILENT)
    assert channel.states == ("m",)
    assert channel.kinetics_missing


def test_procedure_kinetics_are_captured_in_evaluation_order() -> None:
    """The SLO iso variants compute their gating in a PROCEDURE, not a FUNCTION.

    Order is content here: ``v0`` uses ``s0``, and ``minf`` uses both. Sorting
    these -- as the formulas table is sorted -- would silently produce
    expressions referring to names that are not yet bound.
    """
    channel = parse("slodemo", PROCEDURE_KINETICS)
    assert [name for name, _ in channel.procedure] == ["s0", "v0", "minf"]
    assert not channel.kinetics_missing

    rendered = to_toml(channel)
    assert "[procedure]" in rendered
    assert 'order = ["s0", "v0", "minf"]' in rendered


def test_the_procedure_survives_serialisation() -> None:
    """It once did not: the guard accepted the channel because the PROCEDURE had
    been parsed, and the writer omitted it, so the file on disk had states and an
    empty [formulas] block. The check and the output disagreed."""
    rendered = to_toml(parse("slodemo", PROCEDURE_KINETICS))
    assert "minf" in rendered.split("[procedure.expressions]")[1]


def test_constants_carry_the_line_they_came_from() -> None:
    """The point of parsing rather than retyping: any value can be checked."""
    channel = parse("demo", NESTED_SIGNATURE)
    value, unit, line = channel.parameters["va"]
    assert value == -17.6053
    assert unit == "(mV)"
    assert NESTED_SIGNATURE.splitlines()[line - 1].strip().startswith("va =")


def test_the_suffix_names_the_channel_not_the_filename() -> None:
    assert parse("whatever", NESTED_SIGNATURE).name == "demo"
    assert parse("whatever", NESTED_SIGNATURE).ion == "k"


@pytest.mark.parametrize("source", [NESTED_SIGNATURE, PASSIVE, PROCEDURE_KINETICS])
def test_rendered_toml_is_parseable(source: str) -> None:
    import tomllib

    tomllib.loads(to_toml(parse("x", source)))
