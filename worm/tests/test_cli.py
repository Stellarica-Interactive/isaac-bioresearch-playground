"""Smoke tests for the inspection CLI.

Every subcommand is exercised in ``--format json`` so that the machine-readable
output is actually kept working, and so that a number copied out of it always
carries its own provenance.
"""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest

from common.neural.synapses import UnknownSignPolicy
from tools.inspect_connectome import EXIT_FAIL, EXIT_OK, main

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
RUNNER = REPO_ROOT / "worm" / "isaac" / "run_connectome.py"
#: Both Isaac runners. run_body.py carried the same falsy-zero bug as
#: run_connectome.py -- --torque-scale 0, its own null control, silently ran the
#: default -- so neither is checked alone.
RUNNERS = (RUNNER, REPO_ROOT / "worm" / "isaac" / "run_body.py")

DATASET = "witvliet_2021_7"


def run(capsys: pytest.CaptureFixture[str], *argv: str) -> tuple[int, dict]:
    code = main(["--format", "json", *argv])
    out = capsys.readouterr().out
    return code, json.loads(out) if out.strip() else {}


class TestCommands:
    def test_list(self, capsys: pytest.CaptureFixture[str]) -> None:
        code, payload = run(capsys, "list")
        assert code == EXIT_OK
        assert any(d["id"] == DATASET for d in payload["datasets"])

    def test_summary(self, capsys: pytest.CaptureFixture[str]) -> None:
        code, payload = run(capsys, "--dataset", DATASET, "summary")
        assert code == EXIT_OK
        assert payload["totals"]["neurons"] == 181
        assert payload["scope"] == "head"

    def test_verify_passes_for_a_dataset_with_an_oracle(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        code, payload = run(capsys, "--dataset", DATASET, "verify")
        assert code == EXIT_OK
        assert payload["verified"] is True
        assert payload["fields"]["chemical_weight"] == {"published": 7467, "parsed": 7467}

    def test_validate(self, capsys: pytest.CaptureFixture[str]) -> None:
        code, payload = run(capsys, "--dataset", DATASET, "validate")
        assert code == EXIT_OK
        assert payload["violations"] == []

    def test_cell(self, capsys: pytest.CaptureFixture[str]) -> None:
        code, payload = run(capsys, "--dataset", DATASET, "cell", "ASHL")
        assert code == EXIT_OK
        assert payload["cell"]["roles"] == ["sensory"]
        assert payload["cell"]["neurotransmitters"] == ["glutamate"]
        assert payload["chemical_out"]

    def test_unknown_cell_fails_cleanly(self, capsys: pytest.CaptureFixture[str]) -> None:
        assert main(["--dataset", DATASET, "cell", "NOPE"]) == EXIT_FAIL

    def test_edges_filtering(self, capsys: pytest.CaptureFixture[str]) -> None:
        code, payload = run(
            capsys, "--dataset", DATASET, "edges", "--pre", "ASHL", "--type", "chemical"
        )
        assert code == EXIT_OK
        assert payload["edges"]
        assert all(e["pre"] == "ASHL" and e["type"] == "chemical" for e in payload["edges"])

    def test_role(self, capsys: pytest.CaptureFixture[str]) -> None:
        code, payload = run(capsys, "--dataset", DATASET, "role", "sensory")
        assert code == EXIT_OK
        assert "ASHL" in payload["cells"]

    def test_gaps_lists_unmeasured_fields(self, capsys: pytest.CaptureFixture[str]) -> None:
        code, payload = run(capsys, "--dataset", DATASET, "gaps")
        assert code == EXIT_OK
        assert payload["unmeasured_by_source"]

    def test_coverage(self, capsys: pytest.CaptureFixture[str]) -> None:
        code, payload = run(capsys, "--dataset", DATASET, "coverage")
        assert code == EXIT_OK
        assert all(c["matched"] == c["eligible"] for c in payload["coverage"])

    def test_compare(self, capsys: pytest.CaptureFixture[str]) -> None:
        code, payload = run(capsys, "compare", "witvliet_2021_7", "witvliet_2021_8")
        assert code == EXIT_OK
        assert payload["edges_shared"] > 0
        assert payload["cells_only_a"] or payload["cells_only_b"]

    @pytest.mark.parametrize(
        "metric", ["degree", "in-degree", "out-degree", "betweenness", "components", "reciprocity"]
    )
    def test_graph_metrics(self, capsys: pytest.CaptureFixture[str], metric: str) -> None:
        code, _ = run(capsys, "--dataset", DATASET, "graph", "--metric", metric, "--top", "5")
        assert code == EXIT_OK

    def test_annotations_none(self, capsys: pytest.CaptureFixture[str]) -> None:
        code, payload = run(capsys, "--dataset", DATASET, "--annotations", "none", "summary")
        assert code == EXIT_OK
        assert payload["roles"] == {}


class TestProvenanceStamp:
    @pytest.mark.parametrize("command", [["summary"], ["cell", "AVAL"], ["verify"], ["gaps"]])
    def test_json_output_always_carries_its_citation(
        self, capsys: pytest.CaptureFixture[str], command: list[str]
    ) -> None:
        _, payload = run(capsys, "--dataset", DATASET, *command)
        assert payload["dataset_id"] == DATASET
        assert "Witvliet" in payload["citation"]
        assert payload["schema_version"]


class TestRunConnectomeHelp:
    """``--help`` on the experiment runner, which was broken on main.

    argparse formats every help string with the ``%`` operator, so a literal
    percent sign must be doubled. ``--seed-wave`` said "97% one signal", and
    ``% o`` is a valid conversion -- the space flag applied to octal -- so
    argparse raised ``TypeError: %o format: an integer is required, not dict``
    and ``--help`` exited 1 for every flag in the file.

    It went unnoticed because the runner is always invoked with real arguments,
    and nothing imported it: it launches Isaac Sim at import time, so it can only
    be tested as a subprocess. It is worth the subprocess -- ``--help`` is how
    anyone finds the thirty-nine flags, several of which exist to keep results
    honest, and it returns in about 0.2 s because argparse exits before Isaac
    starts.
    """

    @staticmethod
    def _help() -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(RUNNER), "--help"],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=REPO_ROOT,
        )

    def test_help_succeeds(self) -> None:
        result = self._help()
        assert result.returncode == 0, result.stderr[-2000:]

    def test_every_help_string_is_rendered(self) -> None:
        """A bare ``%`` would fail before this, but a stray ``%%`` would survive
        into the output, so check the rendered text rather than only the code."""
        out = self._help().stdout
        assert "%%" not in out
        for flag in ("--unknown-sign", "--seed-wave", "--lesion", "--drag-ratio"):
            assert flag in out, f"{flag} missing from --help"

    def test_unknown_sign_offers_every_policy(self) -> None:
        """The flag exists so the choice is explicit; a missing option would
        quietly make one of them unreachable, as --drag-ratio 0 once was."""
        out = self._help()
        for policy in UnknownSignPolicy:
            assert policy.value in out.stdout


@pytest.mark.parametrize("runner", RUNNERS, ids=lambda p: p.name)
def test_both_isaac_runners_render_their_help(runner: pathlib.Path) -> None:
    """Neither launches Isaac to answer --help, so this costs about 0.2 s each."""
    result = subprocess.run(
        [sys.executable, str(runner), "--help"],
        capture_output=True,
        text=True,
        timeout=120,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, result.stderr[-2000:]
    assert "%%" not in result.stdout
    assert "--torque-scale" in result.stdout
