"""Citations are data, and are checked like data.

The network-dependent verification lives in CI (``tools/verify_citations.py``,
which resolves every DOI against CrossRef). These tests are the offline half:
they check that citations are *present, well-formed and checkable*, so that a new
one cannot be added in a shape the verifier would silently skip.

The motivating incident: a citation for Fenyves et al. 2020 was once written from
memory with five co-authors who are not on the paper. It read entirely plausibly.
That is exactly the failure this project's rules exist to prevent, and free text
is where it hides.
"""

from __future__ import annotations

import re

import pytest

from tools.generate_neuron_reference import CURATED
from tools.verify_citations import KNOWN_CROSSREF_ERRATA, collect
from worm.importers.sources import source_manifest

DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
YEAR_RE = re.compile(r"\(\d{4}\)")


class TestCuratedNotes:
    def test_every_note_has_a_resolvable_doi(self) -> None:
        """Free text alone cannot be machine-checked; a DOI can."""
        missing = sorted(k for k, note in CURATED.items() if not note.doi)
        assert not missing, f"curated notes with no DOI: {missing}"

    def test_dois_are_well_formed(self) -> None:
        bad = sorted(k for k, note in CURATED.items() if not DOI_RE.match(note.doi))
        assert not bad, f"malformed DOIs: {bad}"

    def test_every_note_carries_a_year(self) -> None:
        bad = sorted(k for k, note in CURATED.items() if not YEAR_RE.search(note.source))
        assert not bad, f"sources with no year, so not really citations: {bad}"


class TestSourceManifest:
    def test_every_published_source_has_a_doi(self) -> None:
        """A source we cannot resolve is a source nobody can check."""
        missing = [
            key
            for key, entry in source_manifest().items()
            if isinstance(entry, dict)
            and entry.get("citation")
            and not entry.get("doi")
            # WormAtlas is a living web resource with no DOI; it cites Cook 2019,
            # which does have one, and that is the checkable part.
            and key != "openworm_cect_all_cell_info"
        ]
        assert not missing, f"sources with a citation but no DOI: {missing}"

    def test_manifest_dois_are_well_formed(self) -> None:
        for key, entry in source_manifest().items():
            if isinstance(entry, dict) and entry.get("doi"):
                assert DOI_RE.match(str(entry["doi"])), key


class TestCollection:
    def test_the_verifier_finds_a_substantial_number_of_citations(self) -> None:
        """Guards against a parser change silently checking almost nothing."""
        citations = collect()
        assert len(citations) > 60

    def test_most_collected_citations_carry_a_doi(self) -> None:
        citations = collect()
        with_doi = [c for c in citations if c.doi]
        assert len(with_doi) / len(citations) > 0.7

    @pytest.mark.parametrize("doi", sorted(KNOWN_CROSSREF_ERRATA))
    def test_every_erratum_explains_itself(self, doi: str) -> None:
        """An excused mismatch must say why, or it is just a silenced failure."""
        assert DOI_RE.match(doi)
        for field, reason in KNOWN_CROSSREF_ERRATA[doi].items():
            assert field.strip()
            assert len(reason) > 80, f"{doi}/{field}: reason is too thin to audit"
