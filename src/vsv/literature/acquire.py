"""Download article XML for a cited PMID, or record why it was not kept."""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.bibliography import efetch_pmc_url, europepmc_fulltext_url, oai_url
from vsv.literature.catalog import CitedPaper, ManifestRow
from vsv.literature.xmlbody import BODY_MIN, body_characters, license_sentence

FetchText = Callable[[str], str]


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
