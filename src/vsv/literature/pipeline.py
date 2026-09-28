"""Fetch the PubMed corpus and extract design fields from the local files."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.audit import audit_messages
from vsv.literature.cluster import clusters_from, unclustered_pmids
from vsv.literature.extract import describe
from vsv.literature.fetch import (
    Fetch,
    FetchText,
    clean_text,
    europepmc_records,
    fulltext_url,
    live_fetch,
    live_fetch_text,
    methods_text,
    pubmed_count,
    pubmed_ids,
)
from vsv.literature.papers import (
    digest_pmids,
    load_audit,
    load_extractions,
    load_ledger,
    load_papers,
    write_extractions,
    write_ledger,
    write_papers,
)
from vsv.literature.query import FULLTEXT_NAME, PUBMED_TERM, RESIDUAL_TERM, VECTOR_KINDS
from vsv.literature.records import Extraction, Ledger, Paper

Log = Callable[[str], None]


@dataclass(frozen=True, slots=True)
class FetchReport:
    """What one complete ingest changed."""

    pubmed_count: int
    stored: int
    new: int
    removed: int
    europepmc_misses: int
    residual_count: int


@dataclass(frozen=True, slots=True)
class ExtractReport:
    """What one local extraction wrote."""

    papers: int
    vectors: int
    unclustered: int
    fulltext: int


def ingest_corpus(
    root: Path,
    *,
    fetch: Fetch = live_fetch,
    page_size: int = 9000,
    batch_size: int = 80,
    log: Log | None = None,
) -> FetchReport:
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
        Ledger(
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
    return FetchReport(count, len(stored), len(new_ids), len(removed), len(all_misses), residual)


def extract_corpus(
    root: Path,
    *,
    fulltext: bool = False,
    fetch_text: FetchText = live_fetch_text,
    log: Log | None = None,
) -> ExtractReport:
    """Classify every local paper and extract fields for vector papers.

    With ``fulltext``, open-access methods are read from the local cache or fetched once.
    """
    say = log or _silent
    papers = load_papers(root)
    extractions: list[Extraction] = []
    fetched = 0
    for index, paper in enumerate(papers, start=1):
        methods = ""
        if fulltext and paper.open_access and paper.pmcid:
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
    return ExtractReport(
        len(extractions),
        vectors,
        len(unclustered_pmids(extractions, clusters)),
        fetched,
    )


def coverage_report(root: Path) -> tuple[str, bool]:
    """Text of the coverage ledger and whether every PMID is stored."""
    ledger = load_ledger(root)
    papers = load_papers(root)
    stored = {paper.pmid for paper in papers}
    expected = set(ledger.pmids)
    missing = sorted(expected - stored, key=int)
    extra = sorted(stored - expected, key=int)
    digest_ok = digest_pmids(ledger.pmids) == ledger.digest
    lines = [
        f"query  {ledger.query}",
        f"retrieved  {ledger.retrieved}",
        f"pubmed  {ledger.pubmed_count}",
        f"stored  {len(stored)}",
        f"europepmc misses  {len(ledger.europepmc_misses)}",
        f"open access  {sum(1 for paper in papers if paper.open_access)}",
        f"residual  {ledger.residual_count}",
        f"digest  {'ok' if digest_ok else 'mismatch'}",
    ]
    ok = (
        digest_ok
        and not missing
        and not extra
        and len(stored) == ledger.pubmed_count
        and len(ledger.pmids) == ledger.pubmed_count
    )
    if missing:
        lines.append("missing  " + ", ".join(missing[:20]))
        ok = False
    if extra:
        lines.append("extra  " + ", ".join(extra[:20]))
        ok = False
    try:
        extractions = load_extractions(root)
    except InputError:
        lines.append("extractions  absent")
        return "\n".join(lines), ok
    counts: dict[str, int] = {}
    for item in extractions:
        counts[item.kind] = counts.get(item.kind, 0) + 1
    lines.append("extractions  " + str(len(extractions)))
    for kind in sorted(counts):
        lines.append(f"  {kind}  {counts[kind]}")
    clusters = clusters_from(extractions)
    lines.append(f"clusters  {len(clusters)}")
    lines.append(f"unclustered vectors  {len(unclustered_pmids(extractions, clusters))}")
    extracted = {item.pmid for item in extractions}
    uncovered = sorted(stored - extracted, key=int)
    if uncovered:
        lines.append(f"papers without an extraction  {len(uncovered)}")
        ok = False
    labels = load_audit(root)
    if labels is None:
        lines.append("audit  unlabeled")
    else:
        messages = audit_messages(extractions, labels)
        lines.append(f"audit  {len(labels)} labeled, {len(messages)} misses")
        lines.extend(f"  {message}" for message in messages[:40])
    return "\n".join(lines), ok


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
