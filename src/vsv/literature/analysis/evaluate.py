"""Audit sampling, comparison with manual labels, and analysis coverage."""

from __future__ import annotations

import random
from collections.abc import Mapping, Sequence
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.analysis.cluster import clusters_from, unclustered_pmids
from vsv.literature.analysis.records import VECTOR_KINDS, PaperExtraction
from vsv.literature.records import Paper
from vsv.literature.storage.corpus import load_audit, load_extractions

AUDIT_SEED = 20260921
AUDIT_SIZE = 100


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
    extractions: Sequence[PaperExtraction],
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


def extraction_coverage(root: Path, stored: set[str]) -> tuple[list[str] | None, bool]:
    """Lines describing extraction, cluster, and audit coverage.

    ``None`` means the extraction file is absent. The boolean is false when a stored
    paper has no extraction.
    """
    try:
        extractions = load_extractions(root)
    except InputError:
        return None, True
    lines: list[str] = []
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
    ok = True
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
    return lines, ok
