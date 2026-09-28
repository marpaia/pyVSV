"""Synchronize the PubMed collection and enrich it from Europe PMC."""

from __future__ import annotations

import urllib.parse
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.documents.text import clean_text
from vsv.literature.http import Fetch, live_fetch
from vsv.literature.records import Paper
from vsv.literature.storage.corpus import (
    PubMedSnapshot,
    digest_pmids,
    load_papers,
    write_ledger,
    write_papers,
)
from vsv.literature.storage.paths import corpus_dir

Log = Callable[[str], None]

PUBMED_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EUROPEPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# Title/abstract phrase plus the Indiana and New Jersey MeSH headings.
# A later ingest diffs this term against PubMed and records the count in the ledger.
PUBMED_TERM = (
    '"vesicular stomatitis virus"[tiab] OR '
    '"Vesicular stomatitis Indiana virus"[mh] OR '
    '"Vesicular stomatitis New Jersey virus"[mh]'
)

# Papers that say only "VSV" and are not already in PUBMED_TERM.
# The count is recorded. The hits are not part of the corpus.
RESIDUAL_TERM = f'"VSV"[tiab] NOT ({PUBMED_TERM})'

_PAGE = 1000


@dataclass(frozen=True, slots=True)
class SyncReport:
    """What one complete ingest changed."""

    pubmed_count: int
    stored: int
    new: int
    removed: int
    europepmc_misses: int
    residual_count: int


def sync_pubmed(
    root: Path,
    *,
    fetch: Fetch = live_fetch,
    page_size: int = 9000,
    batch_size: int = 80,
    log: Log | None = None,
) -> SyncReport:
    """Replace the pubmed snapshot with the current PubMed result.

    Files are written under ``root/corpus``. An existing Europe PMC record is kept
    when its PMID is still in the result.
    """
    return ingest_corpus(
        corpus_dir(root),
        fetch=fetch,
        page_size=page_size,
        batch_size=batch_size,
        log=log,
    )


def ingest_corpus(
    root: Path,
    *,
    fetch: Fetch = live_fetch,
    page_size: int = 9000,
    batch_size: int = 80,
    log: Log | None = None,
) -> SyncReport:
    """Replace the corpus with the current PubMed result for the frozen query.

    Existing Europe PMC records are kept for PMIDs that are still in the result. New PMIDs are
    fetched. PMIDs that left the result are dropped.
    """
    say = log or _silent
    count, pmids = pubmed_ids(PUBMED_TERM, fetch, page_size=page_size)
    say(f"pubmed {count}")
    residual = pubmed_count(RESIDUAL_TERM, fetch)
    previous: dict[str, Paper] = {}
    try:
        previous = {paper.pmid: paper for paper in load_papers(root)}
    except InputError:
        previous = {}
    wanted = set(pmids)
    removed = sorted(set(previous) - wanted, key=int)
    new_ids = [pmid for pmid in pmids if pmid not in previous]
    retry = [pmid for pmid in pmids if pmid in previous and not previous[pmid].in_europepmc]
    fresh = _papers_from_europepmc([*new_ids, *retry], fetch, batch_size=batch_size, log=say)
    found = {paper.pmid: paper for paper in fresh}
    retrieved = date.today().isoformat()
    stored = []
    for pmid in pmids:
        if pmid in found:
            stored.append(found[pmid])
        elif pmid in previous:
            stored.append(previous[pmid])
        else:
            stored.append(_stub(pmid, retrieved))
    write_papers(stored, root)
    all_misses = tuple(
        sorted(
            {
                *(paper.pmid for paper in stored if not paper.in_europepmc),
            },
            key=int,
        )
    )
    write_ledger(
        PubMedSnapshot(
            PUBMED_TERM,
            RESIDUAL_TERM,
            retrieved,
            count,
            residual,
            tuple(pmids),
            all_misses,
            digest_pmids(pmids),
        ),
        root,
    )
    return SyncReport(count, len(stored), len(new_ids), len(removed), len(all_misses), residual)


def pubmed_url(term: str, retstart: int, retmax: int) -> str:
    """An E-utilities esearch URL for ``term``."""
    return (
        PUBMED_SEARCH
        + "?"
        + urllib.parse.urlencode(
            {
                "db": "pubmed",
                "retmode": "json",
                "term": term,
                "retstart": str(retstart),
                "retmax": str(retmax),
                "tool": "vsv",
            }
        )
    )


def pubmed_ids(
    term: str,
    fetch: Fetch,
    *,
    page_size: int = 9000,
) -> tuple[int, list[str]]:
    """Every PMID for ``term``, and the PubMed count."""
    start = 0
    total = 0
    found: list[str] = []
    while True:
        payload = fetch(pubmed_url(term, start, page_size), None)
        count, page = _esearch(payload)
        total = count
        found.extend(page)
        if not page or len(found) >= total:
            break
        start += len(page)
    unique = sorted(set(found), key=int)
    if len(unique) != total:
        raise InputError(f"PubMed returned {len(unique)} identifiers for a count of {total}.")
    return total, unique


def pubmed_count(term: str, fetch: Fetch) -> int:
    """The PubMed count for ``term``, without storing identifiers."""
    count, _ids = _esearch(fetch(pubmed_url(term, 0, 0), None))
    return count


def europepmc_url(pmids: Sequence[str], cursor: str) -> str:
    """Search URL for one page of Europe PMC records for ``pmids``.

    The search service accepts GET. A batch stays small enough for the query string.
    """
    for pmid in pmids:
        if not pmid.isdigit():
            raise InputError(f"Invalid PMID {pmid!r}.")
    joined = " OR ".join(f"EXT_ID:{pmid}" for pmid in pmids)
    return (
        EUROPEPMC_SEARCH
        + "?"
        + urllib.parse.urlencode(
            {
                "query": f"({joined}) AND SRC:MED",
                "format": "json",
                "resultType": "core",
                "pageSize": str(_PAGE),
                "cursorMark": cursor,
            }
        )
    )


def europepmc_records(
    pmids: Sequence[str],
    fetch: Fetch,
    *,
    batch_size: int = 80,
    log: Callable[[str], None] | None = None,
) -> list[dict[str, object]]:
    """Core Europe PMC records for ``pmids``. Missing PMIDs are absent from the list.

    A short page is requested again in smaller batches. A PMID that is still absent from a
    one-identifier request is a real miss.
    """
    records: list[dict[str, object]] = []
    done = 0
    for batch in _batches(pmids, batch_size):
        cursor = "*"
        while True:
            payload = fetch(europepmc_url(batch, cursor), None)
            page, cursor, stop = _europepmc_page(payload, cursor)
            records.extend(page)
            if stop:
                break
        done += len(batch)
        if log is not None:
            log(f"europepmc {done}/{len(pmids)}")
    found = {str(record.get("pmid") or "") for record in records}
    missing = [pmid for pmid in pmids if pmid not in found]
    if missing and batch_size > 1:
        if log is not None:
            log(f"retry {len(missing)}")
        records.extend(
            europepmc_records(missing, fetch, batch_size=max(1, batch_size // 4), log=log)
        )
    return records


def _papers_from_europepmc(
    pmids: list[str],
    fetch: Fetch,
    *,
    batch_size: int,
    log: Log,
) -> list[Paper]:
    if not pmids:
        return []
    retrieved = date.today().isoformat()
    records = europepmc_records(pmids, fetch, batch_size=batch_size, log=log)
    papers = []
    for record in records:
        paper = _paper(record, retrieved)
        if paper is not None:
            papers.append(paper)
    return papers


def _paper(record: dict[str, object], retrieved: str) -> Paper | None:
    pmid = str(record.get("pmid") or "")
    if not pmid.isdigit():
        return None
    doi = record.get("doi")
    pmcid = record.get("pmcid")
    year = record.get("pubYear")
    return Paper(
        pmid,
        doi if isinstance(doi, str) else "",
        pmcid if isinstance(pmcid, str) else "",
        str(year) if isinstance(year, int | str) and str(year).isdigit() else "",
        tuple(_named_list(record.get("pubTypeList"), "pubType")),
        tuple(_mesh(record.get("meshHeadingList"))),
        clean_text(str(record.get("title") or "")),
        clean_text(str(record.get("abstractText") or "")),
        str(record.get("isOpenAccess") or "").casefold() == "y",
        f"https://europepmc.org/article/MED/{pmid}",
        retrieved,
        True,
    )


def _stub(pmid: str, retrieved: str) -> Paper:
    return Paper(
        pmid,
        "",
        "",
        "",
        (),
        (),
        "",
        "",
        False,
        f"https://pubmed.ncbi.nlm.nih.gov/{pmid}",
        retrieved,
        False,
    )


def _named_list(value: object, key: str) -> list[str]:
    if not isinstance(value, dict):
        return []
    items = value.get(key)
    if isinstance(items, str):
        return [items]
    if isinstance(items, list):
        return [item for item in items if isinstance(item, str)]
    return []


def _mesh(value: object) -> list[str]:
    if not isinstance(value, dict):
        return []
    headings = value.get("meshHeading")
    if isinstance(headings, dict):
        headings = [headings]
    if not isinstance(headings, list):
        return []
    names: list[str] = []
    for heading in headings:
        if not isinstance(heading, dict):
            continue
        name = heading.get("descriptorName")
        if isinstance(name, str) and name not in names:
            names.append(name)
    return names


def _esearch(payload: Mapping[str, object]) -> tuple[int, list[str]]:
    result = payload.get("esearchresult")
    if not isinstance(result, dict):
        raise InputError("PubMed returned a response without an esearch result.")
    count = result.get("count")
    raw_ids = result.get("idlist", [])
    if not isinstance(count, str) or not count.isdigit():
        raise InputError("PubMed returned a response without a count.")
    if not isinstance(raw_ids, list) or not all(isinstance(item, str) for item in raw_ids):
        raise InputError("PubMed returned a response without a PMID list.")
    return int(count), list(raw_ids)


def _europepmc_page(
    payload: Mapping[str, object], cursor: str
) -> tuple[list[dict[str, object]], str, bool]:
    result_list = payload.get("resultList")
    results: list[dict[str, object]] = []
    if isinstance(result_list, dict) and isinstance(result_list.get("result"), list):
        for item in result_list["result"]:
            if isinstance(item, dict):
                results.append(item)
    nxt = payload.get("nextCursorMark")
    if not isinstance(nxt, str) or nxt == cursor or not results:
        return results, cursor, True
    return results, nxt, False


def _batches(pmids: Sequence[str], size: int) -> list[tuple[str, ...]]:
    return [tuple(pmids[start : start + size]) for start in range(0, len(pmids), size)]


def _silent(_message: str) -> None:
    return None
