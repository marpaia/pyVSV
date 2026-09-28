"""Discovery of collections and the snapshots that record their members."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.catalog import CensusIndex
from vsv.literature.papers import load_papers
from vsv.literature.pipeline import FetchReport, ingest_corpus
from vsv.literature.records import Paper
from vsv.literature.run import make_index
from vsv.literature.sources import Fetch, FetchJson, FetchText, live_fetch
from vsv.literature.store import census_dir, corpus_dir

Log = Callable[[str], None]


def sync_pubmed(
    root: Path,
    *,
    fetch: Fetch = live_fetch,
    page_size: int = 9000,
    batch_size: int = 80,
    log: Log | None = None,
) -> FetchReport:
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


def sync_review_references(
    root: Path,
    fetch_json: FetchJson,
    fetch_text: FetchText,
    email: str,
    retrieved: str,
    papers: dict[str, Paper] | None = None,
) -> CensusIndex:
    """Refresh review-reference membership from the review bibliographies.

    A PMID that is absent from ``papers`` still receives a bibliographic record.
    Whether it also belongs to the pubmed collection is recorded separately.
    """
    known = papers
    if known is None:
        try:
            known = {paper.pmid: paper for paper in load_papers(corpus_dir(root))}
        except InputError:
            known = {}
    return make_index(census_dir(root), known, fetch_json, fetch_text, email, retrieved)
