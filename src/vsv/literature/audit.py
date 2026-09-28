"""Hand-labeled sample used to measure classification misses, and workspace checks."""

from __future__ import annotations

import hashlib
import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.papers import digest_pmids, ledger_path, load_ledger
from vsv.literature.query import AUDIT_SEED, AUDIT_SIZE, VECTOR_KINDS
from vsv.literature.records import PUBMED, Extraction, Paper
from vsv.literature.store import Catalog, corpus_dir, load_catalog


def stratified_pmids(
    papers: Sequence[Paper],
    *,
    size: int = AUDIT_SIZE,
    seed: int = AUDIT_SEED,
) -> list[str]:
    """``size`` PMIDs spread across decades of publication.

    The draw is deterministic for a seed. Papers with no year share one bucket.
    """
    buckets: dict[str, list[str]] = {}
    for paper in papers:
        decade = (
            f"{paper.year[:3]}0" if len(paper.year) >= 3 and paper.year[:3].isdigit() else "unknown"
        )
        buckets.setdefault(decade, []).append(paper.pmid)
    rng = random.Random(seed)
    for pmids in buckets.values():
        rng.shuffle(pmids)
    order = sorted(buckets)
    chosen: list[str] = []
    while len(chosen) < size and any(buckets.values()):
        for decade in order:
            pending = buckets[decade]
            if not pending or len(chosen) >= size:
                continue
            chosen.append(pending.pop())
    return chosen


def audit_messages(
    extractions: Sequence[Extraction],
    labels: Sequence[Mapping[str, object]],
) -> list[str]:
    """Kind mismatches and vector papers the classifier missed."""
    by_pmid = {item.pmid: item for item in extractions}
    messages: list[str] = []
    for label in labels:
        pmid = label.get("pmid")
        kind = label.get("kind")
        vector = label.get("vector")
        if not isinstance(pmid, str) or not isinstance(kind, str) or not isinstance(vector, bool):
            raise InputError("An audit label needs pmid, kind, and vector.")
        found = by_pmid.get(pmid)
        if found is None:
            messages.append(f"{pmid}: labeled paper is not in the extractions.")
            continue
        if found.kind != kind:
            messages.append(f"{pmid}: labeled {kind}, classified {found.kind}.")
        if vector and found.kind not in VECTOR_KINDS:
            messages.append(f"{pmid}: labeled vector, classified {found.kind}.")
    return messages


@dataclass(frozen=True, slots=True)
class Consistency:
    """Problems found in one workspace. An empty list is a consistent catalog."""

    messages: tuple[str, ...]

    @property
    def ok(self) -> bool:
        """Whether the workspace has no consistency problems."""
        return not self.messages


def check_workspace(root: Path) -> Consistency:
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
    return Consistency(tuple(messages))


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
