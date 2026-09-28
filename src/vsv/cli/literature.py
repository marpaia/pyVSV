"""Commands for discovering, acquiring, and tracking papers."""

from __future__ import annotations

import json
import os
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import click

from vsv.cli.corpus import echo_papers
from vsv.design.library import load_library
from vsv.errors import InputError
from vsv.literature.catalog import Catalog, load_catalog
from vsv.literature.checks import check_workspace
from vsv.literature.discovery.pubmed import sync_pubmed
from vsv.literature.discovery.references import sync_review_references
from vsv.literature.documents.acquire import fetch_collection
from vsv.literature.documents.parse import parse_collection
from vsv.literature.http import live_fetch, live_fetch_text
from vsv.literature.records import PUBMED, REVIEW_REFERENCES, Document, Paper
from vsv.literature.search import matching_papers
from vsv.literature.storage.corpus import load_extractions
from vsv.literature.storage.paths import corpus_dir, default_root

_STATUS_ORDER = (
    "xml-body",
    "front-matter-only",
    "no-pmc",
    "browser-html",
    "browser-pdf",
    "unavailable",
    "needs-login",
    "error",
)


@dataclass(frozen=True, slots=True)
class _Call:
    root: Path
    as_json: bool


def _workspace(function: Callable[..., None]) -> Callable[..., None]:
    """Add the shared ``--root`` and ``--json`` options."""
    function = click.option("--json", "as_json", is_flag=True, help="Print JSON.")(function)
    return click.option(
        "--root",
        type=click.Path(path_type=Path),
        default=None,
        help="Workspace directory. Defaults to data/.",
    )(function)


def _call(ctx: click.Context, root: Path | None, as_json: bool) -> _Call:
    parent = ctx.obj if isinstance(ctx.obj, _Call) else _Call(default_root(), False)
    chosen = parent.root if root is None else root
    return _Call(chosen, as_json or parent.as_json)


def _emit(payload: object, as_json: bool, text: str) -> None:
    if as_json:
        click.echo(json.dumps(payload, indent=2, ensure_ascii=False))
        return
    if text:
        click.echo(text)


def _email() -> str:
    email = os.environ.get("VSV_NCBI_EMAIL", "").strip()
    if not email:
        raise InputError("Set VSV_NCBI_EMAIL to a contact address for NCBI requests.")
    return email


def _live_json(url: str) -> dict[str, object]:
    return live_fetch(url, None)


def _aliases(plan: str | None) -> tuple[str, ...]:
    if plan is None:
        return ()
    matched = [record for record in load_library() if record.architecture.id == plan]
    if len(matched) != 1:
        raise InputError(f"{plan}: design not found.")
    return matched[0].aliases


@click.group(invoke_without_command=True)
@_workspace
@click.pass_context
def literature(ctx: click.Context, root: Path | None, as_json: bool) -> None:
    """Discover, find, inspect, acquire, and track papers."""
    ctx.obj = _call(ctx, root, as_json)
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@literature.command("status")
@_workspace
@click.pass_context
def literature_status(ctx: click.Context, root: Path | None, as_json: bool) -> None:
    """Show stored records and work that is still open."""
    call = _call(ctx, root, as_json)
    payload = _status(load_catalog(call.root), call.root)
    _emit(payload, call.as_json, _status_text(payload))


@literature.command("collections")
@_workspace
@click.pass_context
def literature_collections(ctx: click.Context, root: Path | None, as_json: bool) -> None:
    """List collection snapshots in the workspace."""
    call = _call(ctx, root, as_json)
    payload = _collection_rows(load_catalog(call.root))
    if call.as_json:
        _emit(payload, True, "")
        return
    if not payload:
        click.echo("no collections")
        return
    click.echo(_collections_text(payload))


@literature.group("sync", invoke_without_command=True)
@click.pass_context
def sync(ctx: click.Context) -> None:
    """Update bibliographic records and collection membership."""
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@sync.command("pubmed")
@_workspace
@click.pass_context
def sync_pubmed_command(ctx: click.Context, root: Path | None, as_json: bool) -> None:
    """Fetch the frozen PubMed VSV result into the pubmed collection."""
    call = _call(ctx, root, as_json)
    report = sync_pubmed(call.root, log=None if call.as_json else click.echo)
    payload = {
        "pubmed": report.pubmed_count,
        "stored": report.stored,
        "new": report.new,
        "removed": report.removed,
        "europepmc_misses": report.europepmc_misses,
        "residual": report.residual_count,
    }
    text = "\n".join(
        [
            f"pubmed: {report.pubmed_count}",
            f"stored: {report.stored}",
            f"new: {report.new}",
            f"removed: {report.removed}",
            f"europepmc misses: {report.europepmc_misses}",
            f"residual: {report.residual_count}",
        ]
    )
    _emit(payload, call.as_json, text)


@sync.command("review-references")
@_workspace
@click.pass_context
def sync_reviews_command(ctx: click.Context, root: Path | None, as_json: bool) -> None:
    """Fetch references cited by the construction-paper reviews."""
    call = _call(ctx, root, as_json)
    email = _email()
    sync_review_references(
        call.root,
        _live_json,
        live_fetch_text,
        email,
        date.today().isoformat(),
    )
    catalog = load_catalog(call.root)
    review = catalog.named(REVIEW_REFERENCES)
    pubmed = catalog.named(PUBMED)
    review_ids = set(review.pmids) if review is not None else set()
    pubmed_ids = set(pubmed.pmids) if pubmed is not None else set()
    papers = len(review_ids)
    in_pubmed = len(review_ids & pubmed_ids)
    with_pmcid = sum(1 for pmid in review_ids if (paper := catalog.paper(pmid)) and paper.pmcid)
    payload = {"papers": papers, "in_pubmed": in_pubmed, "with_pmcid": with_pmcid}
    text = "\n".join(
        [
            f"papers: {papers}",
            f"in pubmed: {in_pubmed}",
            f"with pmcid: {with_pmcid}",
        ]
    )
    _emit(payload, call.as_json, text)


@literature.command("papers")
@_workspace
@click.option("--collection", default=None, help="Collection name. Defaults to every paper.")
@click.option("--class", "kind", default=None, help="Paper class.")
@click.option("--year", default=None, help="Publication year.")
@click.option("--mesh", default=None, help="MeSH term substring.")
@click.option("--text", default=None, help="Title or abstract substring.")
@click.option("--antigen", default=None, help="Antigen name or slug substring.")
@click.option("--unresolved", is_flag=True, help="Papers the classifier could not place.")
@click.option("--unclustered", is_flag=True, help="Vector papers that did not join a design.")
@click.option("--plan", default=None, help="Curated design id. Matches the plan's aliases.")
@click.option("--limit", default=40, show_default=True, help="Maximum rows. 0 lists every match.")
@click.pass_context
def literature_papers(
    ctx: click.Context,
    root: Path | None,
    as_json: bool,
    collection: str | None,
    kind: str | None,
    year: str | None,
    mesh: str | None,
    text: str | None,
    antigen: str | None,
    unresolved: bool,
    unclustered: bool,
    plan: str | None,
    limit: int,
) -> None:
    """Search the catalog. A collection limits the search to that snapshot."""
    call = _call(ctx, root, as_json)
    catalog = load_catalog(call.root)
    members = _members(catalog, collection)
    try:
        extractions = load_extractions(corpus_dir(call.root))
    except InputError:
        extractions = None
    matched = matching_papers(
        catalog.papers,
        extractions,
        members=members,
        year=year,
        mesh=mesh,
        text=text,
        kind=kind,
        antigen=antigen,
        unresolved=unresolved,
        unclustered=unclustered,
        aliases=_aliases(plan),
        use_aliases=plan is not None,
        show_kind=extractions is not None,
    )
    if call.as_json:
        shown = matched[:limit] if limit > 0 else matched
        payload = {
            "total": len(matched),
            "papers": [
                {
                    "year": paper.year,
                    "pmid": paper.pmid,
                    "kind": label,
                    "title": paper.title,
                    "collections": list(catalog.collections_of(paper.pmid)),
                }
                for paper, label in shown
            ],
        }
        _emit(payload, True, "")
        return
    echo_papers(matched, limit)


@literature.command("show")
@_workspace
@click.argument("pmid")
@click.pass_context
def literature_show(ctx: click.Context, root: Path | None, as_json: bool, pmid: str) -> None:
    """Show one paper, the collections it belongs to, and its documents."""
    call = _call(ctx, root, as_json)
    catalog = load_catalog(call.root)
    paper = catalog.paper(pmid)
    if paper is None:
        raise InputError(f"{pmid}: paper not found.")
    payload = _show_payload(catalog, paper)
    _emit(payload, call.as_json, _show_text(payload))


@literature.command("fetch")
@_workspace
@click.option("--collection", required=True, help="Collection whose documents to acquire.")
@click.option(
    "--limit",
    default=0,
    show_default=True,
    help="Stop after this many papers. 0 fetches every paper in the collection.",
)
@click.pass_context
def literature_fetch(
    ctx: click.Context,
    root: Path | None,
    as_json: bool,
    collection: str,
    limit: int,
) -> None:
    """Acquire documents for a collection. Existing files outside the run stay."""
    call = _call(ctx, root, as_json)
    email = _email() if collection == REVIEW_REFERENCES else ""
    counts = fetch_collection(
        call.root,
        collection,
        fetch_text=live_fetch_text,
        email=email,
        retrieved=date.today().isoformat(),
        log=None if call.as_json else click.echo,
        limit=limit,
    )
    payload = {"collection": collection, "counts": _count_map(counts)}
    _emit(payload, call.as_json, _counts_text(counts))


@literature.command("parse")
@_workspace
@click.option("--collection", required=True, help="Collection whose local documents to parse.")
@click.pass_context
def literature_parse(
    ctx: click.Context,
    root: Path | None,
    as_json: bool,
    collection: str,
) -> None:
    """Parse local documents in a collection into text blocks."""
    call = _call(ctx, root, as_json)
    written = parse_collection(call.root, collection)
    _emit({"readings": written}, call.as_json, f"readings: {written}")


@literature.command("check")
@_workspace
@click.pass_context
def literature_check(ctx: click.Context, root: Path | None, as_json: bool) -> None:
    """Validate the workspace and return a failing status when it does not agree."""
    call = _call(ctx, root, as_json)
    result = check_workspace(call.root)
    payload = {"ok": result.ok, "messages": list(result.messages)}
    text = "ok" if result.ok else "\n".join(result.messages)
    _emit(payload, call.as_json, text)
    if not result.ok:
        ctx.exit(1)


def _members(catalog: Catalog, collection: str | None) -> set[str] | None:
    if collection is None:
        return None
    selected = catalog.named(collection)
    if selected is None:
        raise InputError(f"{collection}: unknown collection.")
    return set(selected.pmids)


def _status(catalog: Catalog, root: Path) -> dict[str, object]:
    pubmed = catalog.named(PUBMED)
    review = catalog.named(REVIEW_REFERENCES)
    pubmed_ids = set(pubmed.pmids) if pubmed is not None else set()
    review_ids = set(review.pmids) if review is not None else set()
    documented = {document.pmid for document in catalog.documents}
    read = {reading.pmid for reading in catalog.readings}
    available = [
        document
        for document in catalog.documents
        if document.availability == "available" and document.pmid
    ]
    return {
        "root": str(root),
        "papers": len(catalog.papers),
        "collections": _collection_rows(catalog),
        "in_both": len(pubmed_ids & review_ids),
        "citations": len(catalog.citations),
        "documents": _document_tally(catalog.documents),
        "readings": len(catalog.readings),
        "pending": {
            "no_document": sum(1 for pmid in review_ids if pmid not in documented),
            "retry": sum(
                1
                for document in catalog.documents
                if document.availability in {"error", "needs-login"}
            ),
            "no_reading": sum(1 for document in available if document.pmid not in read),
        },
    }


def _collection_rows(catalog: Catalog) -> list[dict[str, object]]:
    return [
        {
            "name": collection.name,
            "members": len(collection.members),
            "retrieved": collection.retrieved,
            "rule": collection.rule,
            "provenance": collection.provenance,
        }
        for collection in catalog.collections
    ]


def _document_tally(documents: Sequence[Document]) -> dict[str, object]:
    available = {"xml": 0, "html": 0, "pdf": 0}
    counts = {"partial": 0, "unavailable": 0, "needs-login": 0, "error": 0}
    for document in documents:
        if document.availability == "available":
            if document.format in available:
                available[document.format] += 1
            continue
        if document.availability in counts:
            counts[document.availability] += 1
            continue
        counts["unavailable"] += 1
    return {"available": available, **counts}


def _status_text(payload: Mapping[str, object]) -> str:
    collections = payload["collections"]
    documents = payload["documents"]
    pending = payload["pending"]
    lines = [
        f"root  {payload['root']}",
        f"papers  {payload['papers']}",
        f"in both  {payload['in_both']}",
        f"citations  {payload['citations']}",
    ]
    if isinstance(collections, list) and collections:
        lines.append("collections")
        for item in collections:
            if isinstance(item, dict):
                lines.append(f"  {item['name']}  {item['members']}  {item['retrieved']}")
    if isinstance(documents, dict):
        available = documents.get("available")
        lines.append("documents")
        if isinstance(available, dict):
            for name in ("xml", "html", "pdf"):
                lines.append(f"  {name}  {available.get(name, 0)}")
        for name in ("partial", "unavailable", "needs-login", "error"):
            lines.append(f"  {name}  {documents.get(name, 0)}")
    lines.append(f"readings  {payload['readings']}")
    if isinstance(pending, dict):
        lines.append("pending")
        lines.append(f"  no document  {pending['no_document']}")
        lines.append(f"  retry  {pending['retry']}")
        lines.append(f"  no reading  {pending['no_reading']}")
    return "\n".join(lines)


def _collections_text(rows: Sequence[Mapping[str, object]]) -> str:
    blocks: list[str] = []
    for row in rows:
        blocks.append(
            "\n".join(
                [
                    str(row["name"]),
                    f"  members     {row['members']}",
                    f"  retrieved   {row['retrieved']}",
                    f"  rule        {row['rule']}",
                    f"  provenance  {row['provenance']}",
                ]
            )
        )
    return "\n\n".join(blocks)


def _show_payload(catalog: Catalog, paper: Paper) -> dict[str, object]:
    cited_by = sorted(
        (item.source for item in catalog.citations if item.target == paper.pmid),
        key=int,
    )
    cites = sorted(
        (item.target for item in catalog.citations if item.source == paper.pmid),
        key=int,
    )
    documents = [
        {
            "format": document.format,
            "availability": document.availability,
            "method": document.method,
            "location": document.location,
            "sha256": document.sha256,
            "source": document.source,
            "retrieved": document.retrieved,
            "detail": document.detail,
        }
        for document in catalog.documents
        if document.pmid == paper.pmid
    ]
    readings = [
        {"document": reading.document, "blocks": len(reading.blocks)}
        for reading in catalog.readings
        if reading.pmid == paper.pmid
    ]
    return {
        "pmid": paper.pmid,
        "doi": paper.doi,
        "pmcid": paper.pmcid,
        "year": paper.year,
        "title": paper.title,
        "open_access": paper.open_access,
        "collections": list(catalog.collections_of(paper.pmid)),
        "cited_by": cited_by,
        "cites": cites,
        "documents": documents,
        "readings": readings,
    }


def _show_text(payload: Mapping[str, object]) -> str:
    collections = payload["collections"]
    cited_by = payload["cited_by"]
    lines = [
        f"pmid          {payload['pmid']}",
        f"title         {payload['title']}",
        f"year          {payload['year']}",
        f"doi           {payload['doi']}",
        f"pmcid         {payload['pmcid']}",
        f"collections   {', '.join(collections) if isinstance(collections, list) else ''}",
        f"cited by      {', '.join(cited_by) if isinstance(cited_by, list) else ''}",
    ]
    documents = payload["documents"]
    lines.append("documents")
    if isinstance(documents, list) and documents:
        for item in documents:
            if not isinstance(item, dict):
                continue
            lines.append(
                f"  {item['format'] or '-'}  {item['availability']}  "
                f"{item['method'] or '-'}  {item['location']}"
            )
    else:
        lines.append("  none")
    readings = payload["readings"]
    if isinstance(readings, list):
        blocks = sum(item["blocks"] for item in readings if isinstance(item, dict))
        lines.append(f"readings      {len(readings)}  {blocks} blocks")
    return "\n".join(lines)


def _count_map(counts: Mapping[str, int]) -> dict[str, int]:
    names = [name for name in _STATUS_ORDER if counts.get(name, 0)]
    names.extend(name for name in sorted(counts) if name not in names and counts[name])
    return {name: counts[name] for name in names}


def _counts_text(counts: Mapping[str, int]) -> str:
    shown = _count_map(counts)
    if not shown:
        return "xml-body: 0"
    return "\n".join(f"{name}: {shown[name]}" for name in shown)
