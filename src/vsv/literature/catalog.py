"""The catalog assembled from stored papers, snapshots, documents, and text."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from vsv.literature.discovery.reviews import REVIEWS
from vsv.literature.records import (
    PUBMED,
    REVIEW_REFERENCES,
    Citation,
    Collection,
    Document,
    DocumentText,
    Membership,
    Paper,
)
from vsv.literature.storage.corpus import load_ledger, load_papers, papers_path
from vsv.literature.storage.documents import (
    ManifestRow,
    load_documents,
    manifest_path,
    read_manifest,
)
from vsv.literature.storage.paths import census_dir, corpus_dir, default_root
from vsv.literature.storage.references import CitedPaper, index_path, read_index
from vsv.literature.storage.text import load_texts

_REVIEW_RULE = "A PubMed identifier cited by one of the construction-paper reviews."
_REVIEW_PROVENANCE = "; ".join(f"{review.label} (PMID {review.pmid})" for review in REVIEWS)


@dataclass(frozen=True, slots=True)
class Catalog:
    """Papers, collections, citations, documents, and readings in one workspace."""

    papers: tuple[Paper, ...]
    collections: tuple[Collection, ...]
    citations: tuple[Citation, ...]
    documents: tuple[Document, ...]
    readings: tuple[DocumentText, ...]
    corpus_pmids: frozenset[str]

    def named(self, name: str) -> Collection | None:
        """The collection ``name``, or ``None`` when the workspace has no such snapshot."""
        for collection in self.collections:
            if collection.name == name:
                return collection
        return None

    def paper(self, pmid: str) -> Paper | None:
        """The bibliographic record for ``pmid``, or ``None``."""
        for paper in self.papers:
            if paper.pmid == pmid:
                return paper
        return None

    def collections_of(self, pmid: str) -> tuple[str, ...]:
        """Collection names that include ``pmid``, in catalog order."""
        return tuple(
            collection.name
            for collection in self.collections
            if any(member.pmid == pmid for member in collection.members)
        )


def load_catalog(root: Path | None = None) -> Catalog:
    """Read the workspace into one catalog.

    A paper stored outside the PubMed snapshot keeps its bibliographic record.
    Collection membership is read from each snapshot and is not a field of the paper.
    """
    workspace = default_root() if root is None else root
    by_pmid, corpus_pmids = _papers(workspace)
    collections: list[Collection] = []
    citations: list[Citation] = []
    seen: set[tuple[str, str]] = set()

    ledger_file = corpus_dir(workspace) / "ledger.json"
    if ledger_file.is_file():
        ledger = load_ledger(corpus_dir(workspace))
        collections.append(
            Collection(
                PUBMED,
                ledger.query,
                "PubMed esearch",
                tuple(Membership(pmid, ledger.retrieved) for pmid in ledger.pmids),
            )
        )
    elif corpus_pmids:
        retrieved = ""
        if by_pmid:
            retrieved = next(iter(by_pmid.values())).retrieved
        collections.append(
            Collection(
                PUBMED,
                "",
                "papers.jsonl",
                tuple(Membership(pmid, retrieved) for pmid in sorted(corpus_pmids, key=int)),
            )
        )

    index_file = index_path(census_dir(workspace))
    if index_file.is_file():
        index = read_index(census_dir(workspace))
        members: list[Membership] = []
        for cited in index.papers:
            members.append(Membership(cited.pmid, index.retrieved))
            if cited.pmid not in by_pmid:
                by_pmid[cited.pmid] = _paper_from_cited(cited, index.retrieved)
            for source in cited.cited_by:
                edge = (source, cited.pmid)
                if edge not in seen:
                    seen.add(edge)
                    citations.append(Citation(source, cited.pmid))
        collections.append(
            Collection(REVIEW_REFERENCES, _REVIEW_RULE, _REVIEW_PROVENANCE, tuple(members))
        )

    for row in _manifest_rows(workspace):
        if row.pmid and row.pmid not in by_pmid:
            by_pmid[row.pmid] = _paper_from_row(row)
        for source in row.cited_by:
            edge = (source, row.pmid)
            if edge not in seen and source.isdigit() and row.pmid.isdigit():
                seen.add(edge)
                citations.append(Citation(source, row.pmid))

    documents = load_documents(workspace, by_pmid)
    for document in documents:
        if document.pmid and document.pmid not in by_pmid:
            by_pmid[document.pmid] = _paper_from_document(document)
    readings = load_texts(workspace)
    papers = tuple(by_pmid[pmid] for pmid in sorted(by_pmid, key=int))
    citations.sort(key=lambda item: (int(item.source), int(item.target)))
    return Catalog(
        papers,
        tuple(collections),
        tuple(citations),
        documents,
        readings,
        frozenset(corpus_pmids),
    )


def _papers(root: Path) -> tuple[dict[str, Paper], set[str]]:
    path = papers_path(corpus_dir(root))
    if not path.is_file():
        return {}, set()
    papers = load_papers(corpus_dir(root))
    return {paper.pmid: paper for paper in papers}, {paper.pmid for paper in papers}


def _paper_from_cited(cited: CitedPaper, retrieved: str) -> Paper:
    return Paper(
        cited.pmid,
        cited.doi,
        cited.pmcid,
        cited.year,
        (),
        (),
        cited.title,
        "",
        cited.open_access,
        f"https://pubmed.ncbi.nlm.nih.gov/{cited.pmid}/",
        retrieved,
        False,
    )


def _manifest_rows(root: Path) -> tuple[ManifestRow, ...]:
    manifest = manifest_path(census_dir(root))
    if not manifest.is_file():
        return ()
    return read_manifest(census_dir(root))


def _paper_from_row(row: ManifestRow) -> Paper:
    return Paper(
        row.pmid,
        row.doi,
        row.pmcid,
        row.year,
        (),
        (),
        row.title,
        "",
        False,
        row.source,
        row.retrieved,
        False,
    )


def _paper_from_document(document: Document) -> Paper:
    return Paper(
        document.pmid,
        "",
        "",
        "",
        (),
        (),
        "",
        "",
        False,
        document.source,
        document.retrieved,
        False,
    )
