"""Smoke tests for the inspection CLI.

Every subcommand is exercised in ``--format json`` so that the machine-readable
output is actually kept working, and so that a number copied out of it always
carries its own provenance.
"""

from __future__ import annotations

import json

import pytest

from tools.inspect_connectome import EXIT_FAIL, EXIT_OK, main

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
