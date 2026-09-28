"""Download article files, choose among sources, and keep files already acquired."""

from __future__ import annotations

import hashlib
from collections import Counter
from collections.abc import Callable
from pathlib import Path
from urllib.parse import quote

from vsv.errors import InputError
from vsv.literature.documents.jats import BODY_MIN, body_characters, license_sentence
from vsv.literature.http import FetchText
from vsv.literature.records import PUBMED, REVIEW_REFERENCES
from vsv.literature.storage.corpus import load_papers
from vsv.literature.storage.documents import (
    ManifestRow,
    manifest_path,
    read_manifest,
    write_manifest,
)
from vsv.literature.storage.paths import census_dir, corpus_dir, fulltext_dir
from vsv.literature.storage.references import CitedPaper, read_index

Log = Callable[[str], None]

EUROPEPMC_FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
OAI = "https://pmc.ncbi.nlm.nih.gov/api/oai/v1/mh/"
_KEPT = frozenset({"xml-body", "browser-html", "browser-pdf"})


def fulltext_url(pmcid: str) -> str:
    """Europe PMC full-text XML URL for a PMC identifier."""
    identifier = pmcid if pmcid.startswith("PMC") else f"PMC{pmcid}"
    return f"{EUROPEPMC_FULLTEXT}/{identifier}/fullTextXML"


def europepmc_fulltext_url(pmcid: str) -> str:
    """Europe PMC open-access XML URL."""
    return fulltext_url(pmcid)


def efetch_pmc_url(pmcid: str, email: str) -> str:
    """E-utilities efetch URL for one PMC article."""
    numeric = pmcid.removeprefix("PMC")
    return f"{EFETCH}?db=pmc&id={numeric}&retmode=xml&tool=vsv&email={quote(email)}"


def oai_url(pmcid: str) -> str:
    """PMC OAI-PMH URL for one article. Used when efetch does not return XML."""
    numeric = pmcid.removeprefix("PMC")
    identifier = f"oai:pubmedcentral.nih.gov:{numeric}"
    return f"{OAI}?verb=GetRecord&identifier={quote(identifier)}&metadataPrefix=pmc"


def acquire(
    paper: CitedPaper,
    fulltext: Path,
    fetch_text: FetchText,
    email: str,
    retrieved: str,
) -> ManifestRow:
    """Save a usable article body, or a manifest row that says why none was saved.

    A file already at ``fulltext/{PMCID}.xml`` whose body passes the length test is kept
    and not downloaded again. Europe PMC is tried first, then NCBI efetch. OAI-PMH is
    tried only when efetch does not return XML. Front matter is not written over a file.
    """
    if not paper.pmcid:
        return _row(paper, "", retrieved, "", 0, 0, "", "no-pmc", "")
    path = fulltext / f"{paper.pmcid}.xml"
    if path.is_file():
        cached = path.read_text(encoding="utf-8")
        cached_body = body_characters(cached)
        if cached_body is not None and cached_body >= BODY_MIN:
            return _stored(paper, path, "cache", retrieved, cached, cached_body)
    europe = _pull(fetch_text, europepmc_fulltext_url(paper.pmcid))
    if _usable(europe):
        return _write(paper, path, europepmc_fulltext_url(paper.pmcid), retrieved, europe or "")
    efetch_target = efetch_pmc_url(paper.pmcid, email)
    efetch = _pull(fetch_text, efetch_target)
    if _usable(efetch):
        return _write(paper, path, efetch_target, retrieved, efetch or "")
    if efetch is not None and body_characters(efetch) is not None:
        return _row(paper, efetch_target, retrieved, "", 0, 0, "", "front-matter-only", "")
    oai_target = oai_url(paper.pmcid)
    oai = _pull(fetch_text, oai_target)
    if _usable(oai):
        return _write(paper, path, oai_target, retrieved, oai or "")
    if europe is None and efetch is None and oai is None:
        return _row(paper, efetch_target, retrieved, "", 0, 0, "", "error", "download failed")
    return _row(paper, efetch_target, retrieved, "", 0, 0, "", "front-matter-only", "")


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


def _usable(xml: str | None) -> bool:
    if xml is None:
        return False
    count = body_characters(xml)
    return count is not None and count >= BODY_MIN


def _pull(fetch_text: FetchText, url: str) -> str | None:
    """GET ``url`` twice at most. ``None`` when both attempts fail or are not XML."""
    text = ""
    for _attempt in range(2):
        try:
            text = fetch_text(url)
        except InputError:
            text = ""
            continue
        if body_characters(text) is not None:
            return text
    return None


def _write(paper: CitedPaper, path: Path, source: str, retrieved: str, xml: str) -> ManifestRow:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(xml, encoding="utf-8")
    count = body_characters(xml) or 0
    return _stored(paper, path, source, retrieved, xml, count)


def _stored(
    paper: CitedPaper,
    path: Path,
    source: str,
    retrieved: str,
    xml: str,
    count: int,
) -> ManifestRow:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return _row(
        paper,
        source,
        retrieved,
        digest,
        path.stat().st_size,
        count,
        license_sentence(xml),
        "xml-body",
        "",
    )


def _row(
    paper: CitedPaper,
    source: str,
    retrieved: str,
    digest: str,
    nbytes: int,
    count: int,
    license_text: str,
    status: str,
    detail: str,
) -> ManifestRow:
    return ManifestRow(
        paper.pmid,
        paper.pmcid,
        paper.doi,
        paper.title,
        paper.year,
        source,
        retrieved,
        digest,
        nbytes,
        count,
        license_text,
        status,
        detail,
        paper.cited_by,
    )


def _silent(_message: str) -> None:
    return None
