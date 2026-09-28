"""Group extracted papers into architectures.

The cluster key is the stated edit, the antigen, and the stated locus. Unknown locus stays in
the key, so a paper that states a locus does not join one that does not. Unknown backbone and
replication are gaps. Two stated backbones stay apart. Papers with no antigen, or with a
conflicting edit, stay unclustered.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence

from vsv.design.evidence import EvidenceLevel, Host, Outcome
from vsv.literature.analysis.records import VECTOR_KINDS, PaperCluster, PaperExtraction

_HOST_LEVEL = {
    Host.CELL_LINE: EvidenceLevel.IN_VITRO,
    Host.MOUSE: EvidenceLevel.RODENT,
    Host.HAMSTER: EvidenceLevel.RODENT,
    Host.NONHUMAN_PRIMATE: EvidenceLevel.NONHUMAN_PRIMATE,
    Host.HUMAN: EvidenceLevel.HUMAN,
}


def clusters_from(extractions: Sequence[PaperExtraction]) -> list[PaperCluster]:
    """Clusters for vector papers that state an antigen and an unambiguous edit."""
    grouped: dict[tuple[str, ...], list[PaperExtraction]] = {}
    for item in extractions:
        key = _key(item)
        if key is None:
            continue
        grouped.setdefault(key, []).append(item)
    clusters = [
        _cluster(key, items)
        for key, members in grouped.items()
        for items in _backbone_groups(members)
    ]
    return _unique_ids(clusters)


def unclustered_pmids(
    extractions: Sequence[PaperExtraction], clusters: Sequence[PaperCluster]
) -> list[str]:
    """Vector papers that did not join a cluster."""
    placed = {pmid for cluster in clusters for pmid in cluster.pmids}
    return [
        item.pmid for item in extractions if item.kind in VECTOR_KINDS and item.pmid not in placed
    ]


def _key(item: PaperExtraction) -> tuple[str, ...] | None:
    if item.kind not in VECTOR_KINDS or "edit" in item.conflicts or "antigen" in item.conflicts:
        return None
    if not item.edit_activity.known:
        return None
    antigen = (
        item.antigen_slug.value if item.antigen_slug.known else _token(item.antigen_name.value)
    )
    if antigen == "unknown" or not item.antigen_name.known:
        return None
    return (
        item.edit_activity.value,
        item.edit_target.value,
        antigen,
        item.locus.value,
    )


def _backbone_groups(items: list[PaperExtraction]) -> list[list[PaperExtraction]]:
    """Keep two stated backbones apart. An unstated backbone joins a single stated one."""
    stated = {item.backbone.value for item in items if item.backbone.known}
    if len(stated) <= 1:
        return [items]
    groups: dict[str, list[PaperExtraction]] = {value: [] for value in stated}
    groups["unknown"] = []
    for item in items:
        key = item.backbone.value if item.backbone.known else "unknown"
        groups[key].append(item)
    return [group for group in groups.values() if group]


def _stated(items: list[PaperExtraction], value: Callable[[PaperExtraction], str]) -> str:
    stated = {value(item) for item in items}
    stated.discard("unknown")
    if len(stated) == 1:
        return next(iter(stated))
    return "unknown"


def _cluster(key: tuple[str, ...], items: list[PaperExtraction]) -> PaperCluster:
    activity, target, antigen, locus = key
    backbone = _stated(items, lambda item: item.backbone.value)
    replication = _stated(items, lambda item: item.replication.value)
    named = next((item.antigen_name.value for item in items if item.antigen_name.known), antigen)
    slug = next((item.antigen_slug.value for item in items if item.antigen_slug.known), "")
    hosts: list[str] = []
    phenotypes: list[str] = []
    level = EvidenceLevel.NONE
    for item in items:
        for claim in item.claims:
            if claim.outcome != Outcome.REPORTED.value:
                continue
            if claim.host not in hosts:
                hosts.append(claim.host)
            if claim.phenotype not in phenotypes:
                phenotypes.append(claim.phenotype)
            try:
                candidate = _HOST_LEVEL.get(Host(claim.host), EvidenceLevel.NONE)
            except ValueError:
                candidate = EvidenceLevel.NONE
            if candidate.rank > level.rank:
                level = candidate
    gaps = ["parent-unknown"]
    if any(not item.backbone.known for item in items):
        gaps.append("backbone-unknown")
    if replication == "unknown" or any(not item.replication.known for item in items):
        gaps.append("replication-unknown")
    if locus == "unknown":
        gaps.append("locus-unknown")
    pmids = tuple(sorted((item.pmid for item in items), key=int))
    return PaperCluster(
        id=_ident("c", backbone, replication, activity, target, antigen, "locus", locus),
        name=f"rVSV {named}",
        backbone=backbone,
        replication=replication,
        edit_activity=activity,
        edit_target=target,
        antigen_name=named,
        antigen_slug=slug or (antigen if "-" in antigen else ""),
        locus=locus,
        pmids=pmids,
        hosts=tuple(hosts),
        phenotypes=tuple(phenotypes),
        evidence=level.value,
        gaps=tuple(gaps),
        match=None,
    )


def _unique_ids(clusters: list[PaperCluster]) -> list[PaperCluster]:
    seen: dict[str, int] = {}
    unique: list[PaperCluster] = []
    for cluster in sorted(clusters, key=lambda item: item.id):
        count = seen.get(cluster.id, 0) + 1
        seen[cluster.id] = count
        if count == 1:
            unique.append(cluster)
            continue
        unique.append(
            PaperCluster(
                f"{cluster.id}_{count}",
                cluster.name,
                cluster.backbone,
                cluster.replication,
                cluster.edit_activity,
                cluster.edit_target,
                cluster.antigen_name,
                cluster.antigen_slug,
                cluster.locus,
                cluster.pmids,
                cluster.hosts,
                cluster.phenotypes,
                cluster.evidence,
                cluster.gaps,
                cluster.match,
            )
        )
    return unique


def _ident(*parts: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", "_".join(parts))
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    if not cleaned or cleaned[:1].isdigit():
        cleaned = f"c_{cleaned}"
    return cleaned[:120]


def _token(value: str) -> str:
    if value == "unknown":
        return "unknown"
    return _ident(value).casefold()
