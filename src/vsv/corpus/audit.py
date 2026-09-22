"""Hand-labeled sample used to measure classification misses."""

from __future__ import annotations

import random
from collections.abc import Mapping, Sequence

from vsv.corpus.models import Extraction, Paper
from vsv.corpus.query import AUDIT_SEED, AUDIT_SIZE, VECTOR_KINDS
from vsv.errors import InputError


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
