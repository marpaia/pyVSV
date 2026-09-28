"""Build the construction-paper index from review reference lists and PubMed."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from urllib.parse import quote

from vsv.errors import InputError
from vsv.literature.catalog import CensusIndex, CitedPaper
from vsv.literature.fetch import fulltext_url
from vsv.literature.records import Paper
from vsv.literature.reviews import REVIEWS, Review
from vsv.literature.xmlbody import local_name

FetchJson = Callable[[str], Mapping[str, object]]
FetchText = Callable[[str], str]

EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
IDCONV = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
OAI = "https://pmc.ncbi.nlm.nih.gov/api/oai/v1/mh/"
_PAGE = 1000
_ID_BATCH = 200
_PUBMED_BATCH = 100


def references_url(review: Review, page: int) -> str:
    """Europe PMC reference-list URL for one index review."""
    if review.pmcid:
        kind, identifier = "PMC", review.pmcid
    else:
        kind, identifier = "MED", review.pmid
    return f"{EUROPEPMC}/{kind}/{identifier}/references?format=json&pageSize={_PAGE}&page={page}"


def idconv_url(pmids: Sequence[str], email: str) -> str:
    """PMC id-converter URL for a batch of PMIDs."""
    joined = ",".join(pmids)
    return f"{IDCONV}?ids={quote(joined)}&format=json&tool=vsv&email={quote(email)}"


def pubmed_url(pmids: Sequence[str], email: str) -> str:
    """E-utilities efetch URL for PubMed front matter."""
    joined = ",".join(pmids)
    return (
        f"{EFETCH}?db=pubmed&id={quote(joined)}&retmode=xml&rettype=abstract"
        f"&tool=vsv&email={quote(email)}"
    )


def efetch_pmc_url(pmcid: str, email: str) -> str:
    """E-utilities efetch URL for one PMC article."""
    numeric = pmcid.removeprefix("PMC")
    return f"{EFETCH}?db=pmc&id={numeric}&retmode=xml&tool=vsv&email={quote(email)}"


def oai_url(pmcid: str) -> str:
    """PMC OAI-PMH URL for one article. Used when efetch does not return XML."""
    numeric = pmcid.removeprefix("PMC")
    identifier = f"oai:pubmedcentral.nih.gov:{numeric}"
    return f"{OAI}?verb=GetRecord&identifier={quote(identifier)}&metadataPrefix=pmc"


def europepmc_fulltext_url(pmcid: str) -> str:
    """Europe PMC open-access XML URL."""
    return fulltext_url(pmcid)


def build_index(
    fetch_json: FetchJson,
    fetch_text: FetchText,
    papers: Mapping[str, Paper],
    email: str,
    retrieved: str,
) -> CensusIndex:
    """Cited PMIDs, confirmed against the PMC id converter.

    A PMID already in ``papers`` keeps that record's title, year, DOI, and open-access
    flag. The PMCID comes from the id converter. A converter error means the PMID has
    no PMC article.
    """
    if not email:
        raise InputError("Set VSV_NCBI_EMAIL to a contact address for NCBI requests.")
    cited = _cited(fetch_json)
    order = sorted(cited, key=int)
    pmcids = _pmcids(fetch_json, order, email)
    missing = [pmid for pmid in order if pmid not in papers]
    front = _pubmed_front(fetch_text, missing, email) if missing else {}
    records: list[CitedPaper] = []
    for pmid in order:
        known = papers.get(pmid)
        extra = front.get(pmid)
        title = known.title if known is not None else (extra.title if extra else "")
        year = known.year if known is not None else (extra.year if extra else "")
        doi = known.doi if known is not None else (extra.doi if extra else "")
        records.append(
            CitedPaper(
                pmid,
                pmcids.get(pmid, ""),
                doi,
                title,
                year,
                known.open_access if known is not None else False,
                known is not None,
                tuple(sorted(cited[pmid], key=int)),
            )
        )
    return CensusIndex(retrieved, tuple(records))


def _cited(fetch_json: FetchJson) -> dict[str, set[str]]:
    cited: dict[str, set[str]] = {}
    for review in REVIEWS:
        page = 1
        while True:
            payload = fetch_json(references_url(review, page))
            refs = _references(payload)
            for ref in refs:
                if ref.get("source") != "MED":
                    continue
                raw = ref.get("id")
                if not isinstance(raw, str | int):
                    continue
                pmid = str(raw)
                if not pmid.isdigit():
                    continue
                cited.setdefault(pmid, set()).add(review.pmid)
            hit = _hit_count(payload)
            if not refs or page * _PAGE >= hit or page >= 20:
                break
            page += 1
    return cited


def _references(payload: Mapping[str, object]) -> list[Mapping[str, object]]:
    container = payload.get("referenceList")
    if not isinstance(container, dict):
        return []
    return [item for item in _as_list(container.get("reference")) if isinstance(item, dict)]


def _hit_count(payload: Mapping[str, object]) -> int:
    raw = payload.get("hitCount", 0)
    if isinstance(raw, int):
        return raw
    if isinstance(raw, str) and raw.isdigit():
        return int(raw)
    return 0


def _pmcids(fetch_json: FetchJson, pmids: Sequence[str], email: str) -> dict[str, str]:
    found: dict[str, str] = {}
    for start in range(0, len(pmids), _ID_BATCH):
        batch = pmids[start : start + _ID_BATCH]
        payload = fetch_json(idconv_url(batch, email))
        records = payload.get("records")
        if not isinstance(records, list):
            raise InputError("PMC id converter returned a response without records.")
        for item in records:
            if not isinstance(item, dict):
                continue
            requested = item.get("requested-id")
            pmid = item.get("pmid")
            key = requested if isinstance(requested, str) else pmid
            if not isinstance(key, str):
                continue
            pmcid = item.get("pmcid")
            if item.get("status") == "error" or not isinstance(pmcid, str):
                found[key] = ""
            else:
                found[key] = pmcid if pmcid.startswith("PMC") else f"PMC{pmcid}"
    return found


@dataclass(frozen=True, slots=True)
class _Front:
    title: str
    year: str
    doi: str


def _pubmed_front(fetch_text: FetchText, pmids: Sequence[str], email: str) -> dict[str, _Front]:
    found: dict[str, _Front] = {}
    for start in range(0, len(pmids), _PUBMED_BATCH):
        batch = pmids[start : start + _PUBMED_BATCH]
        xml = fetch_text(pubmed_url(batch, email))
        found.update(pubmed_records(xml))
    return found


def pubmed_records(xml: str) -> dict[str, _Front]:
    """Title, year, and DOI from a PubMed efetch article set."""
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as error:
        raise InputError(f"PubMed returned a response that is not XML: {error}") from error
    found: dict[str, _Front] = {}
    for article in root.iter():
        if local_name(article.tag) != "PubmedArticle":
            continue
        pmid = _child_text(article, "PMID")
        if not pmid.isdigit():
            continue
        title = _child_text(article, "ArticleTitle")
        year = _year(article)
        doi = _article_id(article, "doi")
        found[pmid] = _Front(title, year, doi)
    return found


def _year(article: ET.Element) -> str:
    for element in article.iter():
        if local_name(element.tag) != "PubDate":
            continue
        for child in element:
            if local_name(child.tag) == "Year" and (child.text or "").isdigit():
                return child.text or ""
            if local_name(child.tag) == "MedlineDate":
                digits = "".join(ch for ch in (child.text or "") if ch.isdigit())
                if len(digits) >= 4:
                    return digits[:4]
    return ""


def _article_id(article: ET.Element, kind: str) -> str:
    for element in article.iter():
        if local_name(element.tag) != "ArticleId":
            continue
        if element.attrib.get("IdType") == kind:
            return (element.text or "").strip()
    return ""


def _child_text(article: ET.Element, name: str) -> str:
    for element in article.iter():
        if local_name(element.tag) == name and element.text:
            return " ".join(element.text.split())
    return ""


def _as_list(value: object) -> list[object]:
    if isinstance(value, list):
        return value
    if isinstance(value, dict):
        return [value]
    return []
