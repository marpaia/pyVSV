"""PubMed identifiers and Europe PMC metadata for the VSV corpus.

Live calls live in :func:`live_fetch` and :func:`live_fetch_text`. Tests pass their own
callables. Nothing in this module runs at import.
"""

from __future__ import annotations

import html
import json
import re
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections.abc import Callable, Mapping, Sequence

try:
    import certifi
except ImportError:  # a system trust store is enough when certifi is not installed
    certifi = None  # type: ignore[assignment]

from vsv.corpus.query import EUROPEPMC_FULLTEXT, EUROPEPMC_SEARCH, PUBMED_SEARCH
from vsv.errors import InputError
from vsv.ingest import redact

Fetch = Callable[[str, bytes | None], dict[str, object]]
FetchText = Callable[[str], str]

_USER_AGENT = "vsv-corpus"
_PAGE = 1000


def pubmed_url(term: str, retstart: int, retmax: int) -> str:
    """An E-utilities esearch URL for ``term``."""
    return (
        PUBMED_SEARCH
        + "?"
        + urllib.parse.urlencode(
            {
                "db": "pubmed",
                "retmode": "json",
                "term": term,
                "retstart": str(retstart),
                "retmax": str(retmax),
                "tool": "vsv",
            }
        )
    )


def pubmed_ids(
    term: str,
    fetch: Fetch,
    *,
    page_size: int = 9000,
) -> tuple[int, list[str]]:
    """Every PMID for ``term``, and the PubMed count."""
    start = 0
    total = 0
    found: list[str] = []
    while True:
        payload = fetch(pubmed_url(term, start, page_size), None)
        count, page = _esearch(payload)
        total = count
        found.extend(page)
        if not page or len(found) >= total:
            break
        start += len(page)
    unique = sorted(set(found), key=int)
    if len(unique) != total:
        raise InputError(f"PubMed returned {len(unique)} identifiers for a count of {total}.")
    return total, unique


def pubmed_count(term: str, fetch: Fetch) -> int:
    """The PubMed count for ``term``, without storing identifiers."""
    count, _ids = _esearch(fetch(pubmed_url(term, 0, 0), None))
    return count


def europepmc_url(pmids: Sequence[str], cursor: str) -> str:
    """Search URL for one page of Europe PMC records for ``pmids``.

    The search service accepts GET. A batch stays small enough for the query string.
    """
    for pmid in pmids:
        if not pmid.isdigit():
            raise InputError(f"Invalid PMID {pmid!r}.")
    joined = " OR ".join(f"EXT_ID:{pmid}" for pmid in pmids)
    return (
        EUROPEPMC_SEARCH
        + "?"
        + urllib.parse.urlencode(
            {
                "query": f"({joined}) AND SRC:MED",
                "format": "json",
                "resultType": "core",
                "pageSize": str(_PAGE),
                "cursorMark": cursor,
            }
        )
    )


def europepmc_records(
    pmids: Sequence[str],
    fetch: Fetch,
    *,
    batch_size: int = 80,
    log: Callable[[str], None] | None = None,
) -> list[dict[str, object]]:
    """Core Europe PMC records for ``pmids``. Missing PMIDs are absent from the list.

    A short page is requested again in smaller batches. A PMID that is still absent from a
    one-identifier request is a real miss.
    """
    records: list[dict[str, object]] = []
    done = 0
    for batch in _batches(pmids, batch_size):
        cursor = "*"
        while True:
            payload = fetch(europepmc_url(batch, cursor), None)
            page, cursor, stop = _europepmc_page(payload, cursor)
            records.extend(page)
            if stop:
                break
        done += len(batch)
        if log is not None:
            log(f"europepmc {done}/{len(pmids)}")
    found = {str(record.get("pmid") or "") for record in records}
    missing = [pmid for pmid in pmids if pmid not in found]
    if missing and batch_size > 1:
        if log is not None:
            log(f"retry {len(missing)}")
        records.extend(
            europepmc_records(missing, fetch, batch_size=max(1, batch_size // 4), log=log)
        )
    return records


def methods_text(xml: str) -> str:
    """The methods section of a JATS full-text document, with nucleotide runs removed.

    The longest matching section is kept. Nested method sections are shorter than their parent.
    """
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return ""
    texts = ["".join(element.itertext()) for element in root.iter() if _is_methods(element)]
    if not texts:
        return ""
    return clean_text(max(texts, key=len))


def _is_methods(element: ET.Element) -> bool:
    if not str(element.tag).endswith("sec"):
        return False
    title = ""
    for child in list(element):
        if str(child.tag).endswith("title"):
            title = "".join(child.itertext())
            break
    label = f"{element.attrib.get('sec-type', '')} {title}".casefold()
    return "method" in label or "material" in label


def fulltext_url(pmcid: str) -> str:
    """Europe PMC full-text XML URL for a PMC identifier."""
    identifier = pmcid if pmcid.startswith("PMC") else f"PMC{pmcid}"
    return f"{EUROPEPMC_FULLTEXT}/{identifier}/fullTextXML"


def live_fetch(url: str, body: bytes | None = None) -> dict[str, object]:
    """GET or POST JSON, with one retry. Sleeps so a full corpus stays under NCBI's rate."""
    payload = _request(url, body)
    time.sleep(0.35)
    parsed = json.loads(payload)
    if not isinstance(parsed, dict):
        raise InputError("The literature service returned a response that is not an object.")
    return parsed


def live_fetch_text(url: str) -> str:
    """GET a full-text document. Sleeps between calls."""
    text = _request(url, None).decode()
    time.sleep(0.35)
    return text


def clean_text(value: str) -> str:
    """Drop markup and long nucleotide runs from text that will be stored or quoted."""
    stripped = re.sub(r"<[^>]+>", " ", html.unescape(value))
    collapsed = re.sub(r"\s+", " ", stripped).strip()
    return redact(collapsed)


def _request(url: str, body: bytes | None) -> bytes:
    headers = {"User-Agent": _USER_AGENT}
    if body is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    request = urllib.request.Request(url, data=body, headers=headers)
    try:
        return _read(request)
    except urllib.error.URLError:
        pass
    try:
        return _read(request)
    except urllib.error.URLError as error:
        raise InputError(f"Literature request failed: {error}") from error


def _context() -> ssl.SSLContext:
    if certifi is not None:
        return ssl.create_default_context(cafile=certifi.where())
    return ssl.create_default_context()


def _read(request: urllib.request.Request) -> bytes:
    with urllib.request.urlopen(request, timeout=90, context=_context()) as response:
        payload = response.read()
    if not isinstance(payload, bytes):
        raise InputError("The literature service returned a response that is not bytes.")
    return payload


def _esearch(payload: Mapping[str, object]) -> tuple[int, list[str]]:
    result = payload.get("esearchresult")
    if not isinstance(result, dict):
        raise InputError("PubMed returned a response without an esearch result.")
    count = result.get("count")
    raw_ids = result.get("idlist", [])
    if not isinstance(count, str) or not count.isdigit():
        raise InputError("PubMed returned a response without a count.")
    if not isinstance(raw_ids, list) or not all(isinstance(item, str) for item in raw_ids):
        raise InputError("PubMed returned a response without a PMID list.")
    return int(count), list(raw_ids)


def _europepmc_page(
    payload: Mapping[str, object], cursor: str
) -> tuple[list[dict[str, object]], str, bool]:
    result_list = payload.get("resultList")
    results: list[dict[str, object]] = []
    if isinstance(result_list, dict) and isinstance(result_list.get("result"), list):
        for item in result_list["result"]:
            if isinstance(item, dict):
                results.append(item)
    nxt = payload.get("nextCursorMark")
    if not isinstance(nxt, str) or nxt == cursor or not results:
        return results, cursor, True
    return results, nxt, False


def _batches(pmids: Sequence[str], size: int) -> list[tuple[str, ...]]:
    return [tuple(pmids[start : start + size]) for start in range(0, len(pmids), size)]
