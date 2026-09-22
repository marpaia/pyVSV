"""Group extracted papers into architectures.

The cluster key is backbone, replication, edit, antigen, and stated locus. Unknown values stay
in the key, so a paper that states a locus does not join one that does not. Papers with no
antigen, or with a conflicting edit, stay unclustered.
"""

from __future__ import annotations

import json
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from vsv.corpus.models import Extraction, Paper
from vsv.corpus.query import VECTOR_KINDS
from vsv.design.architecture import Architecture, Replication
from vsv.design.edits import Delete, Insert, Replace
from vsv.design.evidence import (
    Assay,
    Design,
    EvidenceLevel,
    Host,
    Implementation,
    Observation,
    Outcome,
    Phenotype,
    Stage,
    is_citation,
)
from vsv.design.library import load_design, load_library
from vsv.design.modules import Antigen, Leader, Module, Reporter, Trailer, ViralGene
from vsv.errors import InputError
from vsv.genes import Gene
from vsv.ingest import redact

_HOST_LEVEL = {
    Host.CELL_LINE: EvidenceLevel.IN_VITRO,
    Host.MOUSE: EvidenceLevel.RODENT,
    Host.HAMSTER: EvidenceLevel.RODENT,
    Host.NONHUMAN_PRIMATE: EvidenceLevel.NONHUMAN_PRIMATE,
    Host.HUMAN: EvidenceLevel.HUMAN,
}
_BACKBONES = {
    "indiana": "VSV-Indiana",
    "new-jersey": "VSV-NewJersey",
    "unknown": "VSV",
}


@dataclass(frozen=True, slots=True)
class Cluster:
    """Papers that share one architecture key."""

    id: str
    name: str
    backbone: str
    replication: str
    edit_activity: str
    edit_target: str
    antigen_name: str
    antigen_slug: str
    locus: str
    pmids: tuple[str, ...]
    hosts: tuple[str, ...]
    phenotypes: tuple[str, ...]
    evidence: str
    gaps: tuple[str, ...]
    match: str | None


def clusters_from(
    extractions: Sequence[Extraction],
    curated: Sequence[Design] | None = None,
) -> list[Cluster]:
    """Clusters for vector papers that state an antigen and an unambiguous edit.

    ``curated`` defaults to the design library.
    """
    records = load_library() if curated is None else curated
    grouped: dict[tuple[str, ...], list[Extraction]] = {}
    for item in extractions:
        key = _key(item)
        if key is None:
            continue
        grouped.setdefault(key, []).append(item)
    clusters = [_cluster(key, items) for key, items in grouped.items()]
    linked = [_with_match(cluster, records) for cluster in clusters]
    return _unique_ids(linked)


def unclustered_pmids(extractions: Sequence[Extraction], clusters: Sequence[Cluster]) -> list[str]:
    """Vector papers that did not join a cluster."""
    placed = {pmid for cluster in clusters for pmid in cluster.pmids}
    return [
        item.pmid for item in extractions if item.kind in VECTOR_KINDS and item.pmid not in placed
    ]


def write_design(
    cluster: Cluster,
    extractions: Sequence[Extraction],
    papers: Sequence[Paper],
    dest: Path,
) -> None:
    """Write one accepted cluster as a design directory.

    The directory is not created when it already exists. Nucleotide runs in quotes are removed.
    """
    if dest.exists():
        raise InputError(f"{dest}: design directory already exists.")
    if cluster.edit_activity == "unknown" or not cluster.antigen_slug:
        raise InputError(f"{cluster.id}: promoting a cluster needs a known edit and antigen class.")
    by_pmid = {item.pmid: item for item in extractions}
    paper_by_pmid = {paper.pmid: paper for paper in papers}
    members = [by_pmid[pmid] for pmid in cluster.pmids if pmid in by_pmid]
    architecture = _architecture(cluster, members)
    implementations, observations = _evidence(cluster, members)
    design = {
        "id": architecture.id,
        "name": architecture.name,
        "backbone": architecture.backbone,
        "replication": architecture.replication.value,
        "modules": [_module_json(module) for module in architecture.modules],
        "parent_id": architecture.parent_id,
        "parent_status": "known" if architecture.parent_id else "unknown",
        "citations": list(architecture.citations),
        "derivation": _derivation_json(architecture),
        "implementations": [
            {
                "id": item.id,
                "stage": item.stage.value,
                "citation": item.citation,
                "note": item.note,
            }
            for item in implementations
        ],
    }
    observation_rows = [
        {
            "id": item.id,
            "implementation_id": item.implementation_id,
            "phenotype": item.phenotype.value,
            "host": item.host.value,
            "assay": item.assay.value,
            "outcome": item.outcome.value,
            "citation": item.citation,
            "note": item.note,
        }
        for item in observations
    ]
    dest.mkdir(parents=True)
    (dest / "design.json").write_text(json.dumps(design, indent=2) + "\n", encoding="utf-8")
    (dest / "observations.json").write_text(
        json.dumps(observation_rows, indent=2) + "\n",
        encoding="utf-8",
    )
    (dest / "PROVENANCE.md").write_text(
        _provenance(cluster, members, paper_by_pmid),
        encoding="utf-8",
    )
    load_design(dest)


def _key(item: Extraction) -> tuple[str, ...] | None:
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
        item.backbone.value,
        item.replication.value,
        item.edit_activity.value,
        item.edit_target.value,
        antigen,
        item.locus.value,
    )


def _cluster(key: tuple[str, ...], items: list[Extraction]) -> Cluster:
    backbone, replication, activity, target, antigen, locus = key
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
    if replication == "unknown":
        gaps.append("replication-unknown")
    if locus == "unknown":
        gaps.append("locus-unknown")
    pmids = tuple(sorted((item.pmid for item in items), key=int))
    return Cluster(
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


def _with_match(cluster: Cluster, curated: Sequence[Design]) -> Cluster:
    for record in curated:
        architecture = record.architecture
        if _matches(cluster, architecture):
            return Cluster(
                cluster.id,
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
                architecture.id,
            )
    return cluster


def _matches(cluster: Cluster, architecture: Architecture) -> bool:
    slug = cluster.antigen_slug
    if not slug or slug not in architecture.antigen_classes:
        return False
    edit = architecture.edit
    if isinstance(edit, Replace) and cluster.edit_activity == "replace-module":
        same = edit.target == cluster.edit_target
    elif isinstance(edit, Insert) and cluster.edit_activity == "insert-module":
        same = edit.after == cluster.edit_target
    elif isinstance(edit, Delete) and cluster.edit_activity == "delete-module":
        same = edit.target == cluster.edit_target
    else:
        return False
    if not same:
        return False
    return (
        cluster.replication == "unknown"
        or architecture.replication is Replication.UNKNOWN
        or cluster.replication == architecture.replication.value
    )


def _unique_ids(clusters: list[Cluster]) -> list[Cluster]:
    seen: dict[str, int] = {}
    unique: list[Cluster] = []
    for cluster in sorted(clusters, key=lambda item: item.id):
        count = seen.get(cluster.id, 0) + 1
        seen[cluster.id] = count
        if count == 1:
            unique.append(cluster)
            continue
        unique.append(
            Cluster(
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


def _architecture(cluster: Cluster, members: Sequence[Extraction]) -> Architecture:
    payload: Module
    slug = cluster.antigen_slug
    if slug.startswith("reporter-"):
        payload = Reporter("reporter", cluster.antigen_name)
    elif slug:
        payload = Antigen("antigen", cluster.antigen_name, slug)
    else:
        raise InputError(f"{cluster.id}: promoting a cluster needs an antigen class slug.")
    modules: list[Module] = [Leader()]
    if cluster.edit_activity == "replace-module":
        for index, gene in enumerate(Gene):
            if gene.value == cluster.edit_target:
                modules.append(payload)
            else:
                modules.append(ViralGene(gene.value, gene.value, locus=index))
        edit: Replace | Insert | Delete = Replace(cluster.edit_target, payload.id)
    elif cluster.edit_activity == "delete-module":
        kept = [gene for gene in Gene if gene.value != cluster.edit_target]
        for index, gene in enumerate(kept):
            modules.append(ViralGene(gene.value, gene.value, locus=index))
        edit = Delete(cluster.edit_target)
    elif cluster.edit_activity == "insert-module":
        if cluster.edit_target == "leader":
            modules.append(payload)
        for index, gene in enumerate(Gene):
            modules.append(ViralGene(gene.value, gene.value, locus=index))
            if gene.value == cluster.edit_target:
                modules.append(payload)
        edit = Insert(cluster.edit_target, payload.id)
    else:
        raise InputError(f"{cluster.id}: unknown edit.")
    placed = any(module.id == payload.id for module in modules)
    if cluster.edit_activity != "delete-module" and not placed:
        raise InputError(f"{cluster.id}: the edit does not place {payload.id}.")
    modules.append(Trailer())
    replication = {
        "replication-competent": Replication.COMPETENT,
        "single-cycle": Replication.SINGLE_CYCLE,
    }.get(cluster.replication, Replication.UNKNOWN)
    citations = tuple(_citations(members))
    parent = "vsv_newjersey" if cluster.backbone == "new-jersey" else "vsv_indiana"
    return Architecture(
        cluster.id,
        cluster.name,
        _BACKBONES.get(cluster.backbone, "VSV"),
        replication,
        tuple(modules),
        parent=parent,
        citations=citations,
        edit=edit,
    )


def _citations(members: Sequence[Extraction]) -> list[str]:
    citations: list[str] = []
    for item in members:
        citation = f"https://pubmed.ncbi.nlm.nih.gov/{item.pmid}"
        if citation not in citations and is_citation(citation):
            citations.append(citation)
        if len(citations) == 25:
            break
    if not citations:
        raise InputError("A promoted design needs a citation.")
    return citations


def _evidence(
    cluster: Cluster, members: Sequence[Extraction]
) -> tuple[list[Implementation], list[Observation]]:
    implementations: list[Implementation] = []
    observations: list[Observation] = []
    for item in members:
        if not item.stage.known:
            continue
        implementation_id = f"pmid_{item.pmid}"
        note = redact(item.edit_activity.quote or item.antigen_name.quote)
        implementations.append(
            Implementation(
                implementation_id,
                Stage(item.stage.value),
                f"https://pubmed.ncbi.nlm.nih.gov/{item.pmid}",
                note,
            )
        )
        for index, claim in enumerate(item.claims, start=1):
            observations.append(
                Observation(
                    f"{implementation_id}_{index}",
                    implementation_id,
                    Phenotype(claim.phenotype),
                    Host(claim.host),
                    Assay(claim.assay),
                    Outcome(claim.outcome),
                    f"https://pubmed.ncbi.nlm.nih.gov/{item.pmid}",
                    redact(claim.quote),
                )
            )
    return implementations, observations


def _module_json(module: Module) -> dict[str, object]:
    payload: dict[str, object] = {
        "id": module.id,
        "kind": module.kind.value,
        "name": module.name,
    }
    if module.locus is not None:
        payload["locus"] = module.locus
    elif module.locus_status.value == "unknown":
        payload["locus_status"] = "unknown"
    antigen_class = module.antigen_class
    if antigen_class is not None:
        payload["antigen_class"] = antigen_class
    return payload


def _derivation_json(architecture: Architecture) -> dict[str, str] | None:
    edit = architecture.edit
    if isinstance(edit, Replace):
        return {"activity": "replace-module", "target": edit.target, "module": edit.module}
    if isinstance(edit, Insert):
        return {"activity": "insert-module", "target": edit.after, "module": edit.module}
    if isinstance(edit, Delete):
        return {"activity": "delete-module", "target": edit.target}
    return None


def _provenance(
    cluster: Cluster,
    members: Sequence[Extraction],
    papers: dict[str, Paper],
) -> str:
    lines = [
        f"# {cluster.name}",
        "",
        (
            f"{len(cluster.pmids)} PubMed records share this backbone, replication class, "
            "edit, antigen, and locus. The edit is recorded on the canonical chassis because "
            "an edit names its parent. A field without a sentence below is unknown. "
            "No nucleotide or amino-acid sequence is stored."
        ),
        "",
        f"Edit: {cluster.edit_activity} {cluster.edit_target}.",
        f"Antigen: {cluster.antigen_name}.",
        f"Locus: {cluster.locus}.",
        f"Replication: {cluster.replication}.",
        f"Backbone: {cluster.backbone}.",
    ]
    if cluster.match is not None:
        lines.append(f"Curated design: {cluster.match}.")
    lines.extend(["", "## Sentences", ""])
    shown = 0
    for item in members:
        if shown >= 12:
            break
        quote = item.edit_activity.quote or item.antigen_name.quote
        if not quote:
            continue
        lines.append(f"PMID {item.pmid}. {redact(quote)}")
        lines.append("")
        shown += 1
    lines.extend(["## Records", ""])
    for pmid in cluster.pmids:
        paper = papers.get(pmid)
        source = paper.source if paper is not None else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}"
        lines.append(f"- PMID {pmid}. {source}")
    lines.append("")
    return "\n".join(lines)


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
