"""Match literature clusters to curated designs and export an accepted cluster."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import sbol3

from vsv.design.architecture import Architecture, Replication
from vsv.design.document import design_to_document
from vsv.design.edits import Delete, Insert, Replace
from vsv.design.evidence import (
    Assay,
    Design,
    Host,
    Implementation,
    Observation,
    Outcome,
    Phenotype,
    Stage,
    is_citation,
)
from vsv.design.library import load_design, load_library
from vsv.design.modules import Antigen, Leader, Module, Reporter, Trailer, ViralGene, at_locus
from vsv.errors import InputError
from vsv.genes import Gene
from vsv.ingest import redact
from vsv.literature.analysis.cluster import clusters_from as grouped_clusters
from vsv.literature.analysis.records import PaperCluster, PaperExtraction
from vsv.literature.records import Paper

_BACKBONES = {
    "indiana": "VSV-Indiana",
    "new-jersey": "VSV-NewJersey",
    "unknown": "VSV",
}


def clusters_from(
    extractions: Sequence[PaperExtraction],
    curated: Sequence[Design] | None = None,
) -> list[PaperCluster]:
    """Clusters for vector papers that state an antigen and an unambiguous edit.

    ``curated`` defaults to the design library. ``match`` names the curated architecture
    when the edit, antigen, locus, backbone, and replication agree.
    """
    records = load_library() if curated is None else curated
    return [_with_match(cluster, records) for cluster in grouped_clusters(extractions)]


def write_design(
    cluster: PaperCluster,
    extractions: Sequence[PaperExtraction],
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
    document = design_to_document(Design(architecture, tuple(observations), tuple(implementations)))
    dest.mkdir(parents=True)
    (dest / "design.ttl").write_text(document.write_string(sbol3.TURTLE), encoding="utf-8")
    (dest / "PROVENANCE.md").write_text(
        _provenance(cluster, members, paper_by_pmid),
        encoding="utf-8",
    )
    load_design(dest)


def _with_match(cluster: PaperCluster, curated: Sequence[Design]) -> PaperCluster:
    for record in curated:
        architecture = record.architecture
        if _matches(cluster, architecture):
            return PaperCluster(
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


def _matches(cluster: PaperCluster, architecture: Architecture) -> bool:
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
    if cluster.locus != "unknown":
        module_id = edit.target if isinstance(edit, Delete) else edit.module
        module = next((item for item in architecture.modules if item.id == module_id), None)
        if module is None or module.locus is None or str(module.locus) != cluster.locus:
            return False
    expected = {"VSV-Indiana": "indiana", "VSV-NewJersey": "new-jersey"}.get(architecture.backbone)
    if cluster.backbone not in {"unknown", expected}:
        return False
    return (
        cluster.replication == "unknown"
        or architecture.replication is Replication.UNKNOWN
        or cluster.replication == architecture.replication.value
    )


class _Plan:
    """Modules of one promoted plan, with the next transcription index."""

    def __init__(self) -> None:
        self.modules: list[Module] = [Leader()]
        self._locus = 0

    def place(self, module: Module) -> None:
        """Append ``module`` at the next locus."""
        self.modules.append(at_locus(module, self._locus))
        self._locus += 1

    def finish(self) -> tuple[Module, ...]:
        """The module list, closed by the trailer."""
        self.modules.append(Trailer())
        return tuple(self.modules)


def _architecture(cluster: PaperCluster, members: Sequence[PaperExtraction]) -> Architecture:
    payload: Module
    slug = cluster.antigen_slug
    if slug.startswith("reporter-"):
        payload = Reporter("reporter", cluster.antigen_name)
    elif slug:
        payload = Antigen("antigen", cluster.antigen_name, slug)
    else:
        raise InputError(f"{cluster.id}: promoting a cluster needs an antigen class slug.")
    plan = _Plan()
    if cluster.edit_activity == "replace-module":
        for gene in Gene:
            if gene.value == cluster.edit_target:
                plan.place(payload)
            else:
                plan.place(ViralGene(gene.value, gene.value))
        edit: Replace | Insert | Delete = Replace(cluster.edit_target, payload.id)
    elif cluster.edit_activity == "delete-module":
        for gene in Gene:
            if gene.value != cluster.edit_target:
                plan.place(ViralGene(gene.value, gene.value))
        edit = Delete(cluster.edit_target)
    elif cluster.edit_activity == "insert-module":
        if cluster.edit_target == "leader":
            plan.place(payload)
        for gene in Gene:
            plan.place(ViralGene(gene.value, gene.value))
            if gene.value == cluster.edit_target:
                plan.place(payload)
        edit = Insert(cluster.edit_target, payload.id)
    else:
        raise InputError(f"{cluster.id}: unknown edit.")
    placed = {module.id for module in plan.modules}
    if cluster.edit_activity != "delete-module" and payload.id not in placed:
        raise InputError(f"{cluster.id}: the edit does not place {payload.id}.")
    modules = plan.finish()
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


def _citations(members: Sequence[PaperExtraction]) -> list[str]:
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
    cluster: PaperCluster, members: Sequence[PaperExtraction]
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


def _provenance(
    cluster: PaperCluster,
    members: Sequence[PaperExtraction],
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
