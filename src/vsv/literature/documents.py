"""Acquisition records and the local files they point at."""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.acquire import acquire
from vsv.literature.catalog import ManifestRow
from vsv.literature.fetch import fulltext_url
from vsv.literature.index import manifest_path, read_index, read_manifest, write_manifest
from vsv.literature.papers import load_papers
from vsv.literature.records import PUBMED, REVIEW_REFERENCES
from vsv.literature.store import census_dir, corpus_dir, fulltext_dir
from vsv.literature.xmlbody import BODY_MIN, body_characters

Log = Callable[[str], None]
FetchText = Callable[[str], str]
_KEPT = frozenset({"xml-body", "browser-html", "browser-pdf"})


def download_xml(
    dest: Path,
    fulltext: Path,
    fetch_text: FetchText,
    email: str,
    retrieved: str,
    log: Log | None = None,
    limit: int = 0,
) -> Counter[str]:
    """Download XML for the index at ``dest`` and update ``manifest.jsonl``.

    Rows this run does not process stay in the manifest. An attempt that stores
    no file leaves an existing acquired row in place.
    """
    if not email:
        raise InputError("Set VSV_NCBI_EMAIL to a contact address for NCBI requests.")
    index = read_index(dest)
    selected = index.papers if limit <= 0 else index.papers[:limit]
    existing = _existing(dest)
    say = log or _silent
    counts: Counter[str] = Counter()
    for paper in selected:
        try:
            attempt = acquire(paper, fulltext, fetch_text, email, retrieved)
        except (InputError, UnicodeError, OSError) as error:
            attempt = ManifestRow(
                paper.pmid,
                paper.pmcid,
                paper.doi,
                paper.title,
                paper.year,
                "",
                retrieved,
                "",
                0,
                0,
                "",
                "error",
                str(error),
                paper.cited_by,
            )
        stored = merge_row(existing.get(paper.pmid), attempt)
        existing[paper.pmid] = stored
        counts[stored.status] += 1
        say(f"{paper.pmid} {stored.status}")
    write_manifest(dest, list(existing.values()))
    return counts


def merge_row(existing: ManifestRow | None, attempt: ManifestRow) -> ManifestRow:
    """Keep a stored file when the new attempt did not write one."""
    if existing is None:
        return attempt
    if attempt.sha256:
        return attempt
    if existing.sha256 or existing.status in _KEPT:
        return existing
    return attempt


def fetch_collection(
    root: Path,
    collection: str,
    *,
    fetch_text: FetchText,
    email: str,
    retrieved: str,
    log: Log | None = None,
    limit: int = 0,
) -> Counter[str]:
    """Acquire documents for ``collection``. Files outside the run are left in place."""
    if collection == REVIEW_REFERENCES:
        return download_xml(
            census_dir(root),
            fulltext_dir(root),
            fetch_text,
            email,
            retrieved,
            log=log,
            limit=limit,
        )
    if collection == PUBMED:
        return _fetch_open_access(root, fetch_text, log=log, limit=limit)
    raise InputError(f"{collection}: unknown collection.")


def _existing(dest: Path) -> dict[str, ManifestRow]:
    if not manifest_path(dest).is_file():
        return {}
    return {row.pmid: row for row in read_manifest(dest)}


def _fetch_open_access(
    root: Path,
    fetch_text: FetchText,
    *,
    log: Log | None,
    limit: int,
) -> Counter[str]:
    say = log or _silent
    papers = [paper for paper in load_papers(corpus_dir(root)) if paper.open_access and paper.pmcid]
    if limit > 0:
        papers = papers[:limit]
    counts: Counter[str] = Counter()
    cache = fulltext_dir(root)
    for paper in papers:
        path = cache / f"{paper.pmcid}.xml"
        if path.is_file() and _usable(path.read_text(encoding="utf-8")):
            counts["xml-body"] += 1
            say(f"{paper.pmid} xml-body")
            continue
        url = fulltext_url(paper.pmcid)
        try:
            xml = fetch_text(url)
        except InputError:
            counts["error"] += 1
            say(f"{paper.pmid} error")
            continue
        if _usable(xml):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(xml, encoding="utf-8")
            counts["xml-body"] += 1
            say(f"{paper.pmid} xml-body")
        else:
            counts["front-matter-only"] += 1
            say(f"{paper.pmid} front-matter-only")
    return counts


def _usable(xml: str) -> bool:
    count = body_characters(xml)
    return count is not None and count >= BODY_MIN


def _silent(_message: str) -> None:
    return None
