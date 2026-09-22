"""Fetch newly published VSV vector papers into a review queue.

Accepted designs stay in ``data/designs``. This module only writes paper records and pending
proposals under ``data/ingest``. Nucleotide runs in titles and abstracts are omitted before
anything is saved. Proposals are not promoted automatically.
"""

from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from vsv.sequence import InputError

INGEST_ROOT = Path(__file__).resolve().parents[2] / "data" / "ingest"
DEFAULT_QUERY = (
    'TITLE_ABS:"vesicular stomatitis" AND TITLE_ABS:(vaccine OR vector) AND TITLE_ABS:recombinant'
)
_EUROPE_PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
_NUCLEOTIDE_RUN = re.compile(r"[ACGTUacgtu]{40,}")
_ANTIGEN_CUES = (
    ("nipah", "fusion", "henipavirus-fusion"),
    ("nipah", "glycoprotein", "henipavirus-glycoprotein"),
    ("lassa", "", "arenavirus-gpc"),
    ("marburg", "", "marburgvirus-glycoprotein"),
    ("ebola", "", "filovirus-glycoprotein"),
    ("zaire", "", "filovirus-glycoprotein"),
)


@dataclass(frozen=True, slots=True)
class IngestReport:
    """How many new papers and pending proposals one run wrote."""

    papers: int
    proposals: int
    cursor: str


def redact(text: str) -> str:
    """Remove long nucleotide runs from text that will be stored."""
    return _NUCLEOTIDE_RUN.sub("[sequence omitted]", text)


def europepmc_search(url: str) -> dict[str, object]:
    """GET a Europe PMC search URL and return the JSON object."""
    request = urllib.request.Request(url, headers={"User-Agent": "vsv-ingest"})
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.loads(response.read().decode())
    if not isinstance(payload, dict):
        raise InputError("Europe PMC returned a response that is not an object.")
    return payload


def _search_url(query: str, cursor: str, since: str | None) -> str:
    terms = query if since is None else f"{query} AND FIRST_PDATE:[{since} TO 3000-01-01]"
    return (
        _EUROPE_PMC
        + "?"
        + urllib.parse.urlencode(
            {
                "query": terms,
                "format": "json",
                "pageSize": "25",
                "cursorMark": cursor,
                "resultType": "core",
            }
        )
    )


def _antigen_class(text: str) -> str | None:
    lowered = text.lower()
    for pathogen, kind, antigen_class in _ANTIGEN_CUES:
        if pathogen in lowered and (not kind or kind in lowered):
            return antigen_class
    return None


def _proposal(paper: dict[str, object]) -> dict[str, object] | None:
    title = redact(str(paper.get("title") or ""))
    abstract = redact(str(paper.get("abstractText") or ""))
    sentence = title or abstract[:240]
    antigen_class = _antigen_class(f"{title} {abstract}")
    if antigen_class is None:
        return None
    pmid = str(paper.get("pmid") or paper.get("id") or "")
    if not pmid:
        return None
    doi = str(paper.get("doi") or "")
    citations = [f"https://doi.org/{doi}"] if doi else []
    return {
        "status": "pending",
        "pmid": pmid,
        "source_sentence": sentence,
        "design": {
            "id": f"pmid_{pmid}",
            "name": sentence[:80],
            "backbone": "VSV",
            "replication": "unknown",
            "parent_status": "unknown",
            "modules": [
                {
                    "id": "antigen_1",
                    "kind": "antigen",
                    "name": antigen_class,
                    "locus_status": "unknown",
                    "antigen_class": antigen_class,
                }
            ],
            "citations": citations,
            "unmodeled": [
                {"key": "abstract", "value": abstract[:500], "source": sentence},
            ],
        },
    }


def _paper_record(paper: dict[str, object], retrieved: str) -> dict[str, object] | None:
    pmid = str(paper.get("pmid") or "")
    if not pmid:
        return None
    doi = str(paper.get("doi") or "")
    return {
        "pmid": pmid,
        "doi": doi,
        "title": redact(str(paper.get("title") or "")),
        "abstract": redact(str(paper.get("abstractText") or "")),
        "retrieved": retrieved,
        "source": f"https://europepmc.org/article/MED/{pmid}",
    }


def ingest(
    dest: Path | None = None,
    *,
    fetch: Callable[[str], dict[str, object]] = europepmc_search,
    query: str = DEFAULT_QUERY,
    since: str | None = None,
    pages: int = 1,
) -> IngestReport:
    """Fetch one or more Europe PMC pages and write papers plus pending proposals.

    ``pages`` limits a single run. The stored cursor is the next page, so a later run continues.
    Existing paper files are left in place.
    """
    root = INGEST_ROOT if dest is None else dest
    papers_dir = root / "papers"
    queue_dir = root / "queue"
    papers_dir.mkdir(parents=True, exist_ok=True)
    queue_dir.mkdir(parents=True, exist_ok=True)
    cursor_path = root / "cursor.json"
    cursor = "*"
    if cursor_path.is_file() and since is None:
        saved = json.loads(cursor_path.read_text())
        if isinstance(saved, dict) and isinstance(saved.get("cursor"), str):
            cursor = saved["cursor"]
    retrieved = date.today().isoformat()
    written_papers = 0
    written_proposals = 0
    for _page in range(pages):
        payload = fetch(_search_url(query, cursor, since))
        result_list = payload.get("resultList")
        results: list[object] = []
        if isinstance(result_list, dict) and isinstance(result_list.get("result"), list):
            results = result_list["result"]
        for item in results:
            if not isinstance(item, dict):
                continue
            record = _paper_record(item, retrieved)
            if record is None:
                continue
            paper_path = papers_dir / f"{record['pmid']}.json"
            if not paper_path.exists():
                paper_path.write_text(json.dumps(record, indent=2) + "\n")
                written_papers += 1
            proposal = _proposal(item)
            if proposal is None:
                continue
            queue_path = queue_dir / f"{record['pmid']}.json"
            if not queue_path.exists():
                queue_path.write_text(json.dumps(proposal, indent=2) + "\n")
                written_proposals += 1
        next_cursor = payload.get("nextCursorMark")
        if not isinstance(next_cursor, str) or next_cursor == cursor or not results:
            cursor = next_cursor if isinstance(next_cursor, str) else cursor
            break
        cursor = next_cursor
    cursor_path.write_text(json.dumps({"cursor": cursor, "retrieved": retrieved}, indent=2) + "\n")
    return IngestReport(written_papers, written_proposals, cursor)
