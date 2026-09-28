"""Workspace paths and the catalog assembled from the stored files."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.catalog import CitedPaper, ManifestRow
from vsv.literature.index import index_path, manifest_path, read_index, read_manifest
from vsv.literature.papers import load_ledger, load_papers, papers_path
from vsv.literature.query import CORPUS_ROOT
from vsv.literature.records import (
    PUBMED,
    REVIEW_REFERENCES,
    Acquisition,
    Block,
    Citation,
    Collection,
    Document,
    Membership,
    Paper,
    Reading,
)
from vsv.literature.reviews import REVIEWS

_REVIEW_RULE = "A PubMed identifier cited by one of the construction-paper reviews."
_REVIEW_PROVENANCE = "; ".join(f"{review.label} (PMID {review.pmid})" for review in REVIEWS)


def default_root() -> Path:
    """The data directory that holds the corpus and the census."""
    return CORPUS_ROOT.parent


def corpus_dir(root: Path) -> Path:
    """PubMed collection files under ``root``."""
    return root / "corpus"


def census_dir(root: Path) -> Path:
    """Review-reference files under ``root``."""
    return root / "census"


def fulltext_dir(root: Path) -> Path:
    """Open-access XML files for the workspace."""
    return corpus_dir(root) / "fulltext"


def acquired_dir(root: Path) -> Path:
    """Publisher files saved for review references."""
    return census_dir(root) / "acquired"


def sections_dir(root: Path) -> Path:
    """Parsed readings for the workspace."""
    return census_dir(root) / "sections"


@dataclass(frozen=True, slots=True)
class Catalog:
    """Papers, collections, citations, documents, and readings in one workspace."""

    papers: tuple[Paper, ...]
    collections: tuple[Collection, ...]
    citations: tuple[Citation, ...]
    documents: tuple[Document, ...]
    readings: tuple[Reading, ...]
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

    documents = _documents(workspace, by_pmid)
    for document in documents:
        if document.pmid and document.pmid not in by_pmid:
            by_pmid[document.pmid] = _paper_from_document(document)
    readings = _readings(workspace)
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


def document_from_row(row: ManifestRow) -> Document:
    """Split a manifest status string into format, method, and availability."""
    kind, availability, method = _decode(row.status, row.sha256, row.source)
    location = _location(row, kind, availability)
    attempt = Acquisition(row.retrieved, row.source, method, availability, row.detail)
    return Document(
        row.pmid,
        kind,
        location,
        row.sha256,
        row.nbytes,
        row.source,
        row.retrieved,
        availability,
        row.license,
        row.body_characters,
        method,
        row.detail,
        (attempt,),
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


def _documents(root: Path, papers: dict[str, Paper]) -> tuple[Document, ...]:
    found: list[Document] = []
    covered: set[str] = set()
    manifest = manifest_path(census_dir(root))
    if manifest.is_file():
        for row in read_manifest(census_dir(root)):
            document = document_from_row(row)
            found.append(document)
            if document.location:
                covered.add(document.location)
    acquired = acquired_dir(root)
    if acquired.is_dir():
        for kind, suffix in (("html", ".html"), ("pdf", ".pdf")):
            for path in sorted(acquired.glob(f"*{suffix}")):
                location = f"census/acquired/{path.name}"
                if location in covered or not path.stem.isdigit():
                    continue
                found.append(_file_document(path, path.stem, kind, location, "browser"))
                covered.add(location)
    fulltext = fulltext_dir(root)
    if fulltext.is_dir():
        by_pmcid = {paper.pmcid: paper.pmid for paper in papers.values() if paper.pmcid}
        for path in sorted(fulltext.glob("*.xml")):
            location = f"corpus/fulltext/{path.name}"
            if location in covered:
                continue
            pmid = by_pmcid.get(path.stem, "")
            found.append(_file_document(path, pmid, "xml", location, "download"))
            covered.add(location)
    found.sort(key=lambda item: (int(item.pmid or 0), item.location))
    return tuple(found)


def _file_document(path: Path, pmid: str, kind: str, location: str, method: str) -> Document:
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    attempt = Acquisition("", "", method, "available", "")
    return Document(
        pmid,
        kind,
        location,
        digest,
        path.stat().st_size,
        "",
        "",
        "available",
        "",
        0,
        method,
        "",
        (attempt,),
    )


def _readings(root: Path) -> tuple[Reading, ...]:
    directory = sections_dir(root)
    if not directory.is_dir():
        return ()
    readings: list[Reading] = []
    for path in sorted(directory.glob("*.json")):
        readings.append(_reading(path))
    readings.sort(key=lambda item: int(item.pmid) if item.pmid.isdigit() else 0)
    return tuple(readings)


def _reading(path: Path) -> Reading:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise InputError(f"{path}: {error.msg}.") from error
    if not isinstance(raw, dict):
        raise InputError(f"{path}: reading must be an object.")
    pmid = raw.get("pmid")
    blocks = raw.get("blocks")
    if not isinstance(pmid, str) or not isinstance(blocks, list):
        raise InputError(f"{path}: reading needs pmid and blocks.")
    version = raw.get("sha256")
    parsed: list[Block] = []
    for item in blocks:
        if not isinstance(item, dict):
            raise InputError(f"{path}: each block must be an object.")
        parsed.append(
            Block(
                _string(item, "id", path),
                _string(item, "section", path),
                _string(item, "kind", path),
                _string(item, "text", path),
            )
        )
    digest = version if isinstance(version, str) else ""
    return Reading(pmid, digest, tuple(parsed))


def _string(raw: dict[str, object], key: str, path: Path) -> str:
    value = raw.get(key)
    if isinstance(value, str):
        return value
    raise InputError(f"{path}: {key} must be a string.")


def _decode(status: str, sha256: str, source: str) -> tuple[str, str, str]:
    if status == "xml-body":
        method = "cache" if source == "cache" else "download"
        return "xml", "available", method
    if status == "browser-html":
        return "html", "available", "browser"
    if status == "browser-pdf":
        return "pdf", "available", "browser"
    if status == "front-matter-only":
        return "xml", "partial", "download"
    if status == "needs-login":
        return "", "needs-login", "browser"
    if status == "unavailable":
        kind = "html" if sha256 else ""
        return kind, "unavailable", "browser"
    if status == "error":
        return "", "error", ""
    if status == "no-pmc":
        return "", "unavailable", ""
    return "", "unavailable", ""


def _location(row: ManifestRow, kind: str, availability: str) -> str:
    if kind == "xml" and availability == "available" and row.pmcid:
        return f"corpus/fulltext/{row.pmcid}.xml"
    if kind == "html" and row.pmid:
        return f"census/acquired/{row.pmid}.html"
    if kind == "pdf" and row.pmid:
        return f"census/acquired/{row.pmid}.pdf"
    return ""
