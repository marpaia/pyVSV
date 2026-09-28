"""Commands for the curated VSV design library."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import click

from vsv.cli.corpus import cluster_edit, resolve_cluster
from vsv.cli.render import cli_name, indent, row, table, widths, wrap
from vsv.design import (
    KNOWLEDGE_GAPS,
    LIBRARY_ROOT,
    Architecture,
    Design,
    EvidenceLevel,
    Host,
    LocusStatus,
    Module,
    Observation,
    Phenotype,
    Plan,
    load_design,
    load_library,
    plans,
    read,
    read_plasmid,
    select,
)
from vsv.design import check as check_design
from vsv.design.edits import Delete, Edit, Insert, Replace
from vsv.errors import InputError
from vsv.genome.payload import AntigenClass
from vsv.ingest import DEFAULT_QUERY, INGEST_ROOT
from vsv.ingest import ingest as fetch_papers
from vsv.literature.cluster import Cluster, clusters_from
from vsv.literature.papers import load_extractions
from vsv.literature.pipeline import ingest_corpus
from vsv.literature.query import CORPUS_ROOT

_EVIDENCE = {
    "none": EvidenceLevel.NONE,
    "in-vitro": EvidenceLevel.IN_VITRO,
    "rodent": EvidenceLevel.RODENT,
    "preclinical": EvidenceLevel.RODENT,
    "nonhuman-primate": EvidenceLevel.NONHUMAN_PRIMATE,
    "human": EvidenceLevel.HUMAN,
}


class EvidenceName(click.ParamType[EvidenceLevel]):
    """An evidence threshold, including the ``preclinical`` alias for rodent."""

    name = "evidence"

    def convert(
        self, value: object, param: click.Parameter | None, ctx: click.Context | None
    ) -> EvidenceLevel:
        if isinstance(value, EvidenceLevel):
            return value
        text = str(value)
        try:
            return _EVIDENCE[text]
        except KeyError:
            self.fail(f"unknown evidence {text!r}.", param, ctx)


def _phenotype(value: str | None) -> Phenotype | None:
    return None if value is None else Phenotype(value)


def _host(value: str | None) -> Host | None:
    return None if value is None else Host(value)


def _design_rows(records: Sequence[Design], clusters: Sequence[Cluster]) -> str:
    rows: list[list[str]] = []
    for record in records:
        card = record.card
        edit = record.architecture.edit
        rows.append(
            [
                record.architecture.id,
                "curated",
                _derivation_text(edit) if edit is not None else "-",
                ",".join(card.antigen_classes) or "-",
                card.replication.value,
                card.evidence_level.value,
                "-",
                _gap_text(record),
            ]
        )
    for cluster in clusters:
        status = f"confirmed:{cluster.match}" if cluster.match else "extracted"
        rows.append(
            [
                cluster.id,
                status,
                cluster_edit(cluster),
                cluster.antigen_slug or cluster.antigen_name,
                cluster.replication,
                cluster.evidence,
                str(len(cluster.pmids)),
                ",".join(cluster.gaps) or "-",
            ]
        )
    if not rows:
        return "no designs match"
    return table(
        ["id", "status", "edit", "antigens", "replication", "evidence", "papers", "gaps"],
        rows,
    )


def _clusters(
    curated: Sequence[Design],
    antigen_class: str | None,
    phenotype: str | None,
    host: str | None,
    evidence: EvidenceLevel | None,
) -> list[Cluster]:
    try:
        extractions = load_extractions()
    except InputError:
        return []
    found = []
    for cluster in clusters_from(extractions, curated):
        if antigen_class not in {None, cluster.antigen_slug, cluster.antigen_name}:
            continue
        if phenotype is not None and phenotype not in cluster.phenotypes:
            continue
        if host is not None and host not in cluster.hosts:
            continue
        if evidence is not None and EvidenceLevel(cluster.evidence).rank < evidence.rank:
            continue
        found.append(cluster)
    return found


def _gap_text(record: Design) -> str:
    names = [
        item.rule
        for item in check_design(record.architecture, record.observations)
        if item.rule in KNOWLEDGE_GAPS
    ]
    return ",".join(names) or "-"


def _parent_text(architecture: Architecture) -> str:
    if architecture.parent_id is not None:
        return architecture.parent_id
    return str(architecture.parent)


def _derivation_text(edit: Edit) -> str:
    match edit:
        case Delete(target=target, replication=replication):
            text = f"delete {cli_name(target)}"
        case Insert(after=after, module=module, replication=replication):
            text = f"insert after {cli_name(after)} with {cli_name(module)}"
        case Replace(target=target, module=module, replication=replication):
            text = f"replace {cli_name(target)}"
            shown = cli_name(module)
            if shown != cli_name(target):
                text += f" with {shown}"
    if replication is not None:
        text += f", {replication.value}"
    return text


def _fact_lines(pairs: Sequence[tuple[str, str]]) -> list[str]:
    width = max(len(key) for key, _value in pairs)
    return [f"{key:<{width}}  {value}" for key, value in pairs]


def _locus(module: Module) -> str:
    if module.locus_status is LocusStatus.NOT_TRANSCRIBED:
        return "-"
    if module.locus_status is LocusStatus.UNKNOWN or module.locus is None:
        return "unknown"
    return str(module.locus)


def _module_table(modules: Sequence[Module]) -> str:
    rows = [[_locus(module), cli_name(module.id), module.antigen_class or ""] for module in modules]
    return table(["locus", "module", "class"], rows)


def _observation_lines(observations: Sequence[Observation], citation: str) -> list[str]:
    if not observations:
        return ["  no observations"]
    header = ["phenotype", "host", "assay", "outcome"]
    rows = [
        [item.phenotype.value, item.host.value, item.assay.value, item.outcome.value]
        for item in observations
    ]
    column_widths = widths(header, rows)
    lines = ["  " + row(header, column_widths)]
    for item, cells in zip(observations, rows, strict=True):
        lines.append("  " + row(cells, column_widths))
        if item.citation != citation:
            lines.append(f"    {item.citation}")
        if item.note:
            lines.extend(wrap(item.note, "    "))
        lines.append("")
    return lines[:-1]


def _show(record: Design) -> str:
    architecture = record.architecture
    facts = [
        ("backbone", architecture.backbone),
        ("replication", architecture.replication.value),
        ("parent", _parent_text(architecture)),
    ]
    if architecture.edit is not None:
        facts.append(("derivation", _derivation_text(architecture.edit)))
    if record.realization:
        facts.append(("realization", record.realization))
    if record.model_source:
        facts.append(("model", record.model_source))
    lines = [
        architecture.id,
        architecture.name,
        "",
        *_fact_lines(facts),
        "",
        "modules",
        indent(_module_table(architecture.modules)),
    ]
    known = {item.id for item in record.implementations}
    grouped: dict[str, list[Observation]] = {item_id: [] for item_id in known}
    for observation in record.observations:
        if observation.implementation_id in known:
            grouped[observation.implementation_id].append(observation)
    if not record.implementations:
        lines.extend(["", "implementations", "  none"])
    for implementation in record.implementations:
        lines.extend(
            [
                "",
                f"{implementation.stage.value}  {implementation.id}",
                f"  {implementation.citation}",
            ]
        )
        if implementation.note:
            lines.extend(wrap(implementation.note, "  "))
        lines.append("")
        lines.extend(_observation_lines(grouped[implementation.id], implementation.citation))
    if record.unmodeled:
        key_width = max(len(field.key) for field in record.unmodeled)
        lines.extend(["", "unmodeled"])
        for field in record.unmodeled:
            lines.append(f"  {field.key:<{key_width}}  {field.value}")
            if field.source:
                lines.extend(wrap(field.source, " " * (key_width + 4)))
    loose = [
        observation.id
        for observation in record.observations
        if observation.implementation_id not in known
    ]
    if loose:
        lines.extend(["", "observations without an implementation", "  " + ", ".join(loose)])
    gaps = [
        item.message
        for item in check_design(architecture, record.observations)
        if item.rule in KNOWLEDGE_GAPS
    ]
    if gaps:
        lines.extend(["", "gaps", *[f"  {gap}" for gap in gaps]])
    return "\n".join(lines)


def _resolve_show(identifier: str) -> Design:
    path = Path(identifier)
    if path.is_dir():
        return load_design(path)
    matches = [record for record in load_library() if record.architecture.id == identifier]
    if len(matches) != 1:
        raise InputError(f"{identifier}: design not found.")
    return matches[0]


@click.command()
@click.option("--antigen-class", default=None, help="Antigen class slug.")
@click.option(
    "--phenotype",
    type=click.Choice([item.value for item in Phenotype]),
    default=None,
    help="Observation phenotype.",
)
@click.option(
    "--host",
    type=click.Choice([item.value for item in Host]),
    default=None,
    help="Observation host.",
)
@click.option(
    "--evidence",
    type=EvidenceName(),
    default=None,
    help="Minimum evidence: in-vitro, rodent, preclinical, nonhuman-primate, or human.",
)
@click.option(
    "--library",
    type=click.Path(path_type=Path),
    default=LIBRARY_ROOT,
    show_default=True,
    help="Design library.",
)
def designs(
    antigen_class: str | None,
    phenotype: str | None,
    host: str | None,
    evidence: EvidenceLevel | None,
    library: Path,
) -> None:
    """List curated and extracted VSV architectures."""
    curated = load_library(library)
    records = select(
        curated,
        antigen_class=antigen_class,
        phenotype=_phenotype(phenotype),
        host=_host(host),
        min_evidence=evidence,
    )
    click.echo(_design_rows(records, _clusters(curated, antigen_class, phenotype, host, evidence)))


@click.command()
@click.argument("design", type=click.Path(path_type=Path))
@click.pass_context
def check(ctx: click.Context, design: Path) -> None:
    """Report design-rule violations for one curated design."""
    record = load_design(design)
    violations = check_design(record.architecture, record.observations)
    structural = [item for item in violations if item.rule not in KNOWLEDGE_GAPS]
    gaps = [item for item in violations if item.rule in KNOWLEDGE_GAPS]
    if not structural and not gaps:
        click.echo(f"{record.architecture.id}: ok")
        return
    blocks = []
    if structural:
        blocks.append(
            table(["rule", "message"], [(item.rule, item.message) for item in structural])
        )
    if gaps:
        blocks.append(table(["gap", "message"], [(item.rule, item.message) for item in gaps]))
    click.echo("\n".join(blocks))
    if structural:
        ctx.exit(1)


@click.command()
@click.argument("plan", type=click.Path(path_type=Path))
@click.option("--antigen-class", default=None, help="Antigen class slug, including subclasses.")
@click.option(
    "--phenotype",
    type=click.Choice([item.value for item in Phenotype]),
    default=None,
    help="Observation phenotype.",
)
@click.option(
    "--host",
    type=click.Choice([item.value for item in Host]),
    default=None,
    help="Observation host.",
)
def precedents(
    plan: Path, antigen_class: str | None, phenotype: str | None, host: str | None
) -> None:
    """List published plans with this plan's arrangement."""
    candidate = _read_plan(plan)
    antigen = None if antigen_class is None else AntigenClass(antigen_class)
    rows: list[list[str]] = []
    for published in candidate.precedents(
        plans(),
        antigen_class=antigen,
        phenotype=_phenotype(phenotype),
        host=_host(host),
    ):
        for result in published.observations:
            if phenotype is not None and result.phenotype is not Phenotype(phenotype):
                continue
            if host is not None and result.host is not Host(host):
                continue
            rows.append(
                [
                    published.id,
                    result.phenotype.value,
                    result.host.value,
                    result.assay.value,
                    result.citation,
                ]
            )
    if not rows:
        click.echo("no precedents")
        return
    click.echo(table(["id", "phenotype", "host", "assay", "citation"], rows))


@click.command()
@click.argument("plasmid", type=click.Path(path_type=Path))
def plasmid(plasmid: Path) -> None:
    """Show the backbone, slot, and insert in a plasmid file."""
    record = read_plasmid(plasmid)
    click.echo(
        "\n".join(
            _fact_lines(
                [
                    ("plan", record.plan_id),
                    ("edit", record.edit),
                    ("slot", record.slot),
                    ("insert", record.insert_id),
                    ("backbone", record.backbone.iri),
                    ("product", record.product.iri),
                ]
            )
        )
    )


def _read_plan(path: Path) -> Plan:
    turtle = path / "design.ttl" if path.is_dir() else path
    return read(turtle)


@click.command()
@click.argument("design")
def show(design: str) -> None:
    """Show one curated design or one extracted cluster."""
    try:
        record = _resolve_show(design)
    except InputError:
        text = resolve_cluster(design)
        if text is None:
            raise
        click.echo(text)
        return
    click.echo(_show(record))


@click.command()
@click.option(
    "--dest",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory. Defaults to data/ingest, or data/corpus with --all.",
)
@click.option("--since", default=None, help="Only papers first published on or after this date.")
@click.option("--query", default=DEFAULT_QUERY, help="Europe PMC query.")
@click.option("--all", "complete", is_flag=True, help="Fetch the frozen PubMed VSV corpus.")
def ingest(dest: Path | None, since: str | None, query: str, complete: bool) -> None:
    """Fetch papers into the review queue, or the whole PubMed corpus with --all."""
    if complete:
        if since is not None or query != DEFAULT_QUERY:
            raise InputError("vsv ingest --all uses the frozen PubMed corpus query.")
        root = CORPUS_ROOT if dest is None else dest
        report = ingest_corpus(root, log=click.echo)
        click.echo(
            "\n".join(
                [
                    f"pubmed: {report.pubmed_count}",
                    f"stored: {report.stored}",
                    f"new: {report.new}",
                    f"removed: {report.removed}",
                    f"europepmc misses: {report.europepmc_misses}",
                    f"residual: {report.residual_count}",
                ]
            )
        )
        return
    queue = fetch_papers(INGEST_ROOT if dest is None else dest, query=query, since=since)
    click.echo(f"papers: {queue.papers}\nproposals: {queue.proposals}\ncursor: {queue.cursor}")
