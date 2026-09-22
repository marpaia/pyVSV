"""Europe PMC ingest writes a review queue and does not promote designs."""

import json
from collections.abc import Iterator
from pathlib import Path

from vsv.ingest import ingest


def _keys(value: object) -> Iterator[str]:
    if isinstance(value, dict):
        yield from value
        for child in value.values():
            yield from _keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from _keys(child)


def test_ingest_queues_partial_proposals_and_omits_nucleotide_runs(tmp_path: Path) -> None:
    calls = {"n": 0}

    def fetch(url: str) -> dict[str, object]:
        calls["n"] += 1
        assert "cursorMark" in url
        if calls["n"] > 1:
            return {"resultList": {"result": []}, "nextCursorMark": "end"}
        abstract = "Recombinant VSV ebola glycoprotein vaccine " + ("ACGT" * 15)
        return {
            "resultList": {
                "result": [
                    {
                        "pmid": "999001",
                        "doi": "10.1000/example",
                        "title": "A recombinant VSV ebola vaccine vector",
                        "abstractText": abstract,
                    },
                    {
                        "pmid": "999002",
                        "doi": "",
                        "title": "Vesicular stomatitis virus ecology",
                        "abstractText": "Field isolates, with no recombinant vector.",
                    },
                    {
                        "pmid": "999003",
                        "doi": "10.1000/niv",
                        "title": "Recombinant vesicular stomatitis Nipah fusion vaccine vector",
                        "abstractText": "A recombinant VSV vector expressing Nipah fusion protein.",
                    },
                ]
            },
            "nextCursorMark": "page-2",
        }

    report = ingest(tmp_path, fetch=fetch)
    assert (report.papers, report.proposals, report.cursor) == (3, 2, "page-2")
    paper = json.loads((tmp_path / "papers" / "999001.json").read_text())
    assert "[sequence omitted]" in paper["abstract"]
    assert "ACGTACGT" not in paper["abstract"]
    proposal = json.loads((tmp_path / "queue" / "999001.json").read_text())
    design = proposal["design"]
    assert proposal["status"] == "pending"
    assert design["id"] == "pmid_999001"
    assert design["replication"] == "unknown"
    assert design["parent_status"] == "unknown"
    assert design["modules"][0]["locus_status"] == "unknown"
    assert design["modules"][0]["antigen_class"] == "filovirus-glycoprotein"
    assert "sequence" not in set(_keys(proposal))
    assert "ACGTACGT" not in json.dumps(proposal)
    fusion = json.loads((tmp_path / "queue" / "999003.json").read_text())
    assert fusion["design"]["modules"][0]["antigen_class"] == "henipavirus-fusion"
    assert not (tmp_path / "queue" / "999002.json").exists()
    assert not (tmp_path / "designs").exists()

    again = ingest(tmp_path, fetch=fetch)
    assert (again.papers, again.proposals) == (0, 0)


def test_since_ignores_the_saved_cursor(tmp_path: Path) -> None:
    (tmp_path / "cursor.json").write_text(json.dumps({"cursor": "saved-page"}))
    seen: list[str] = []

    def fetch(url: str) -> dict[str, object]:
        seen.append(url)
        return {"resultList": {"result": []}, "nextCursorMark": "next"}

    ingest(tmp_path, fetch=fetch, since="2024-01-01")
    assert "saved-page" not in seen[0]
    assert "FIRST_PDATE" in seen[0]
    assert "2024-01-01" in seen[0]
    assert "cursorMark=%2A" in seen[0]
