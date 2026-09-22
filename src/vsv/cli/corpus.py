"""Commands for the PubMed VSV corpus."""

from __future__ import annotations

from pathlib import Path

import click

from vsv.cli.render import cli_name, wrap
from vsv.corpus.cluster import Cluster, clusters_from, write_design
from vsv.corpus.models import Extraction
from vsv.corpus.pipeline import coverage_report, extract_corpus
from vsv.corpus.query import CORPUS_ROOT, VECTOR_KINDS
from vsv.corpus.store import load_extractions, load_papers
from vsv.errors import InputError


def cluster_edit(cluster: Cluster) -> str:
    """The edit as shown on the design list."""
    target = cli_name(cluster.edit_target)
    if cluster.edit_activity == "replace-module":
        return f"replace {target}"
    if cluster.edit_activity == "insert-module":
        return f"insert after {target}"
    if cluster.edit_activity == "delete-module":
        return f"delete {target}"
    return "-"


def format_cluster(cluster: Cluster, extractions: list[Extraction]) -> str:
    """One extracted architecture, its quotes, and its PMIDs."""
    by_pmid = {item.pmid: item for item in extractions}
    status = f"confirmed {cluster.match}" if cluster.match else "extracted"
    lines = [
        cluster.id,
        cluster.name,
        "",
        f"status        {status}",
        f"backbone      {cluster.backbone}",
        f"replication   {cluster.replication}",
        f"edit          {cluster_edit(cluster)}",
        f"antigen       {cluster.antigen_name}",
        f"class         {cluster.antigen_slug or '-'}",
        f"locus         {cluster.locus}",
        f"evidence      {cluster.evidence}",
        f"papers        {len(cluster.pmids)}",
    ]
    if cluster.gaps:
        lines.extend(["", "gaps", *[f"  {gap}" for gap in cluster.gaps]])
    lines.extend(["", "sentences"])
    shown = 0
    for pmid in cluster.pmids:
        item = by_pmid.get(pmid)
        if item is None:
            continue
        quote = item.edit_activity.quote or item.antigen_slug.quote or item.antigen_name.quote
        if not quote:
            continue
        lines.append("")
        lines.append(f"  PMID {pmid}")
        lines.extend(wrap(quote, "  "))
        shown += 1
        if shown == 8:
            break
    lines.extend(["", "pmids", "  " + ", ".join(cluster.pmids[:30])])
    if len(cluster.pmids) > 30:
        lines.append(f"  and {len(cluster.pmids) - 30} more")
    return "\n".join(lines)


def resolve_cluster(identifier: str) -> str | None:
    """Formatted cluster text, or ``None`` when the corpus has no such id."""
    try:
        extractions = load_extractions()
    except InputError:
        return None
    for cluster in clusters_from(extractions):
        if cluster.id == identifier:
            return format_cluster(cluster, extractions)
    return None


@click.group(invoke_without_command=True)
@click.option(
    "--dest",
    type=click.Path(path_type=Path),
    default=None,
    help="Corpus directory. Defaults to data/corpus.",
)
@click.pass_context
def corpus(ctx: click.Context, dest: Path | None) -> None:
    """Show whether the local corpus holds the PubMed VSV result."""
    ctx.obj = CORPUS_ROOT if dest is None else dest
    if ctx.invoked_subcommand is None:
        text, ok = coverage_report(ctx.obj)
        click.echo(text)
        if not ok:
            ctx.exit(1)


@corpus.command("extract")
@click.option("--fulltext", is_flag=True, help="Read open-access methods for vector papers.")
@click.pass_context
def extract_literature(ctx: click.Context, fulltext: bool) -> None:
    """Classify every stored paper and extract quote-backed design fields."""
    root = ctx.obj if isinstance(ctx.obj, Path) else CORPUS_ROOT
    report = extract_corpus(root, fulltext=fulltext, log=click.echo)
    click.echo(
        "\n".join(
            [
                f"papers: {report.papers}",
                f"vectors: {report.vectors}",
                f"unclustered: {report.unclustered}",
                f"fulltext: {report.fulltext}",
            ]
        )
    )


@click.command()
@click.option("--class", "kind", default=None, help="Paper class.")
@click.option("--year", default=None, help="Publication year.")
@click.option("--mesh", default=None, help="MeSH term substring.")
@click.option("--text", default=None, help="Title or abstract substring.")
@click.option("--antigen", default=None, help="Antigen name or slug substring.")
@click.option("--unresolved", is_flag=True, help="Papers the classifier could not place.")
@click.option("--unclustered", is_flag=True, help="Vector papers that did not join a design.")
@click.option("--limit", default=40, show_default=True, help="Maximum rows. 0 lists every match.")
@click.option(
    "--dest",
    type=click.Path(path_type=Path),
    default=None,
    help="Corpus directory. Defaults to data/corpus.",
)
def papers(
    kind: str | None,
    year: str | None,
    mesh: str | None,
    text: str | None,
    antigen: str | None,
    unresolved: bool,
    unclustered: bool,
    limit: int,
    dest: Path | None,
) -> None:
    """Search the local VSV corpus."""
    root = CORPUS_ROOT if dest is None else dest
    records = load_papers(root)
    extractions: dict[str, Extraction] = {}
    if kind or unresolved or unclustered or antigen:
        extractions = {item.pmid: item for item in load_extractions(root)}
    placed: set[str] = set()
    if unclustered:
        placed = {
            pmid for cluster in clusters_from(list(extractions.values())) for pmid in cluster.pmids
        }
    rows = []
    for paper in records:
        item = extractions.get(paper.pmid)
        if kind is not None and (item is None or item.kind != kind):
            continue
        if unresolved and (item is None or item.status != "unresolved"):
            continue
        if unclustered and (item is None or item.kind not in VECTOR_KINDS or paper.pmid in placed):
            continue
        if year is not None and paper.year != year:
            continue
        if mesh is not None and not any(mesh.casefold() in term.casefold() for term in paper.mesh):
            continue
        if (
            text is not None
            and text.casefold() not in f"{paper.title}\n{paper.abstract}".casefold()
        ):
            continue
        if antigen is not None:
            if item is None:
                continue
            blob = f"{item.antigen_name.value}\n{item.antigen_slug.value}"
            if antigen.casefold() not in blob.casefold():
                continue
        label = item.kind if item is not None else "-"
        rows.append((paper.year or "-", paper.pmid, label, paper.title))
    rows.sort(key=lambda row: (row[0], int(row[1])), reverse=True)
    total = len(rows)
    if limit > 0:
        rows = rows[:limit]
    if not rows:
        click.echo("no papers match")
        return
    width = max(len(row[2]) for row in rows)
    for row_year, pmid, label, title in rows:
        click.echo(f"{row_year}  {pmid}  {label.ljust(width)}  {title}")
    if limit > 0 and total > limit:
        click.echo(f"showing {len(rows)} of {total}")


@click.command()
@click.argument("cluster_id")
@click.option(
    "--dest",
    type=click.Path(path_type=Path),
    default=None,
    help="Design directory to write. Defaults to data/designs/<id>.",
)
@click.option(
    "--corpus",
    "corpus_root",
    type=click.Path(path_type=Path),
    default=None,
    help="Corpus directory. Defaults to data/corpus.",
)
def promote(cluster_id: str, dest: Path | None, corpus_root: Path | None) -> None:
    """Copy one extracted cluster into a curated design directory."""
    root = CORPUS_ROOT if corpus_root is None else corpus_root
    extractions = load_extractions(root)
    papers_ = load_papers(root)
    found = next((item for item in clusters_from(extractions) if item.id == cluster_id), None)
    if found is None:
        raise InputError(f"{cluster_id}: cluster not found.")
    target = dest if dest is not None else Path("data/designs") / found.id
    write_design(found, extractions, papers_, target)
    click.echo(str(target))
