"""Commands for the construction-paper bibliography and its article files."""

from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import click

from vsv.errors import InputError
from vsv.literature.fetch import live_fetch, live_fetch_text
from vsv.literature.index import CENSUS_ROOT, read_manifest
from vsv.literature.papers import load_papers
from vsv.literature.query import CORPUS_ROOT, FULLTEXT_NAME
from vsv.literature.run import download_xml, make_index
from vsv.literature.section import section_xml, write_sections


def _live_json(url: str) -> dict[str, object]:
    return live_fetch(url, None)


def _email() -> str:
    email = os.environ.get("VSV_NCBI_EMAIL", "").strip()
    if not email:
        raise InputError("Set VSV_NCBI_EMAIL to a contact address for NCBI requests.")
    return email


def _dest(ctx: click.Context) -> Path:
    dest = ctx.obj
    return dest if isinstance(dest, Path) else CENSUS_ROOT


@click.group(invoke_without_command=True)
@click.option(
    "--dest",
    type=click.Path(path_type=Path),
    default=None,
    help="Census directory. Defaults to data/census.",
)
@click.pass_context
def census(ctx: click.Context, dest: Path | None) -> None:
    """Build the construction-paper bibliography and download article XML."""
    ctx.obj = CENSUS_ROOT if dest is None else dest
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@census.command("index")
@click.pass_context
def census_index(ctx: click.Context) -> None:
    """Write the review bibliography to index.json."""
    email = _email()
    papers = {paper.pmid: paper for paper in load_papers()}
    index = make_index(
        _dest(ctx),
        papers,
        _live_json,
        live_fetch_text,
        email,
        date.today().isoformat(),
    )
    in_corpus = sum(1 for paper in index.papers if paper.in_corpus)
    with_pmcid = sum(1 for paper in index.papers if paper.pmcid)
    click.echo(f"papers: {len(index.papers)}")
    click.echo(f"in corpus: {in_corpus}")
    click.echo(f"with pmcid: {with_pmcid}")


@census.command("fetch")
@click.option(
    "--limit",
    default=0,
    show_default=True,
    help="Stop after this many papers. 0 fetches every paper.",
)
@click.pass_context
def census_fetch(ctx: click.Context, limit: int) -> None:
    """Download article XML for the bibliography into data/corpus/fulltext."""
    counts = download_xml(
        _dest(ctx),
        CORPUS_ROOT / FULLTEXT_NAME,
        live_fetch_text,
        _email(),
        date.today().isoformat(),
        log=click.echo,
        limit=limit,
    )
    for status in ("xml-body", "front-matter-only", "no-pmc", "error"):
        click.echo(f"{status}: {counts.get(status, 0)}")


@census.command("section")
@click.pass_context
def census_section(ctx: click.Context) -> None:
    """Write quotable blocks for every article XML file in the manifest."""
    dest = _dest(ctx)
    written = 0
    for row in read_manifest(dest):
        if row.status != "xml-body" or not row.pmcid:
            continue
        xml_path = CORPUS_ROOT / FULLTEXT_NAME / f"{row.pmcid}.xml"
        blocks = section_xml(xml_path.read_text(encoding="utf-8"))
        write_sections(dest / "sections" / f"{row.pmid}.json", row.pmid, blocks)
        written += 1
    click.echo(f"sections: {written}")
