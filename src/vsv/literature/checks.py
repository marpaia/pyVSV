"""Workspace integrity: membership, snapshot completeness, files, and citations."""

from __future__ import annotations

import hashlib
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from vsv.literature.analysis.evaluate import extraction_coverage
from vsv.literature.catalog import Catalog, load_catalog
from vsv.literature.records import PUBMED
from vsv.literature.storage.corpus import digest_pmids, ledger_path, load_ledger, load_papers
from vsv.literature.storage.paths import corpus_dir


@dataclass(frozen=True, slots=True)
class CheckReport:
    """Problems found in one workspace. An empty list is a consistent catalog."""

    messages: tuple[str, ...]

    @property
    def ok(self) -> bool:
        """Whether the workspace has no consistency problems."""
        return not self.messages


def check_workspace(root: Path) -> CheckReport:
    """Compare snapshots, files, and citations in ``root``.

    A paper that belongs to one collection and not another is consistent.
    A member without a bibliographic record, a broken checksum, and a missing
    file for an available document are not.
    """
    catalog = load_catalog(root)
    messages: list[str] = []
    messages.extend(_pubmed(root, catalog))
    known = {paper.pmid for paper in catalog.papers}
    for collection in catalog.collections:
        if collection.name == PUBMED:
            continue
        missing = [member.pmid for member in collection.members if member.pmid not in known]
        if missing:
            messages.append(_pmids(f"{collection.name}: missing", missing))
    for citation in catalog.citations:
        if citation.target not in known:
            messages.append(
                f"citation {citation.source} -> {citation.target}: no paper for the cited pmid"
            )
    messages.extend(_documents(root, catalog))
    for reading in catalog.readings:
        if reading.pmid not in known:
            messages.append(f"reading {reading.pmid}: no paper")
    return CheckReport(tuple(messages))


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
    analysis, analysis_ok = extraction_coverage(root, stored)
    if analysis is None:
        lines.append("extractions  absent")
        return "\n".join(lines), ok
    lines.extend(analysis)
    return "\n".join(lines), ok and analysis_ok


def _pubmed(root: Path, catalog: Catalog) -> list[str]:
    path = ledger_path(corpus_dir(root))
    if not path.is_file():
        return []
    ledger = load_ledger(corpus_dir(root))
    messages: list[str] = []
    if digest_pmids(ledger.pmids) != ledger.digest:
        messages.append("pubmed: digest mismatch")
    if len(ledger.pmids) != ledger.pubmed_count:
        messages.append("pubmed: count mismatch")
    expected = set(ledger.pmids)
    missing = sorted(expected - catalog.corpus_pmids, key=int)
    extra = sorted(catalog.corpus_pmids - expected, key=int)
    if missing:
        messages.append(_pmids("pubmed: missing", missing))
    if extra:
        messages.append(_pmids("pubmed: extra", extra))
    return messages


def _documents(root: Path, catalog: Catalog) -> list[str]:
    messages: list[str] = []
    base = root.resolve()
    for document in catalog.documents:
        if not document.pmid:
            messages.append(f"document {document.location or '?'}: no paper")
            continue
        if document.availability != "available":
            continue
        if not document.location:
            messages.append(f"document {document.pmid}: available without a file")
            continue
        path = (root / document.location).resolve()
        if not path.is_relative_to(base):
            messages.append(f"document {document.pmid}: location leaves the workspace")
            continue
        if not path.is_file():
            messages.append(f"document {document.pmid}: missing {document.location}")
            continue
        if document.sha256 and hashlib.sha256(path.read_bytes()).hexdigest() != document.sha256:
            messages.append(f"document {document.pmid}: checksum mismatch")
    return messages


def _pmids(prefix: str, pmids: Sequence[str]) -> str:
    shown = ", ".join(pmids[:20])
    if len(pmids) > 20:
        return f"{prefix} {shown} and {len(pmids) - 20} more"
    return f"{prefix} {shown}"
