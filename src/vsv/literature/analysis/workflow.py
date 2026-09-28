"""Run classification and extraction across a stored corpus."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.analysis.cluster import clusters_from, unclustered_pmids
from vsv.literature.analysis.extract import describe
from vsv.literature.analysis.records import VECTOR_KINDS, PaperExtraction
from vsv.literature.documents.acquire import fulltext_url
from vsv.literature.documents.jats import methods_text
from vsv.literature.http import FetchText, live_fetch_text
from vsv.literature.storage.corpus import load_papers, write_extractions
from vsv.literature.storage.paths import FULLTEXT_NAME

Log = Callable[[str], None]


@dataclass(frozen=True, slots=True)
class ExtractionReport:
    """What one local extraction wrote."""

    papers: int
    vectors: int
    unclustered: int
    fulltext: int


def extract_corpus(
    root: Path,
    *,
    fulltext: bool = False,
    fetch_text: FetchText = live_fetch_text,
    log: Log | None = None,
) -> ExtractionReport:
    """Classify every local paper and extract fields for vector papers.

    A saved article next to the corpus, at ``census/acquired/{pmid}.txt``, is used when
    its text is long enough to be a body. With ``fulltext``, open-access methods are
    read from the local cache or fetched once for vector papers that have no saved article.
    """
    say = log or _silent
    papers = load_papers(root)
    extractions: list[PaperExtraction] = []
    fetched = 0
    for index, paper in enumerate(papers, start=1):
        methods = ""
        article = _saved_article(root, paper.pmid)
        if article:
            preview = describe(paper)
            if preview.kind in VECTOR_KINDS:
                methods = article
                fetched += 1
        elif fulltext and paper.open_access and paper.pmcid:
            preview = describe(paper)
            if preview.kind in VECTOR_KINDS:
                methods = _methods(root, paper.pmcid, fetch_text)
                if methods:
                    fetched += 1
        extractions.append(describe(paper, methods))
        if index % 500 == 0:
            say(f"extracted {index}/{len(papers)}")
    write_extractions(extractions, root)
    clusters = clusters_from(extractions)
    vectors = sum(1 for item in extractions if item.kind in VECTOR_KINDS)
    return ExtractionReport(
        len(extractions),
        vectors,
        len(unclustered_pmids(extractions, clusters)),
        fetched,
    )


def _saved_article(root: Path, pmid: str) -> str:
    """Plain text of a saved article for ``pmid``. Empty when the file is missing or short."""
    path = root.parent / "census" / "acquired" / f"{pmid}.txt"
    if not path.is_file():
        return ""
    text = path.read_text(encoding="utf-8")
    if len(" ".join(text.split())) < 2000:
        return ""
    return text


def _methods(root: Path, pmcid: str, fetch_text: FetchText) -> str:
    cache = root / FULLTEXT_NAME
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{pmcid}.xml"
    if path.is_file():
        xml = path.read_text(encoding="utf-8")
    else:
        try:
            xml = fetch_text(fulltext_url(pmcid))
        except InputError:
            return ""
        path.write_text(xml, encoding="utf-8")
    return methods_text(xml)


def _silent(_message: str) -> None:
    return None
