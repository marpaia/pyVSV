"""Bibliography construction and article XML download, with no network."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pytest

from vsv.cli import main
from vsv.errors import InputError
from vsv.literature.acquire import acquire
from vsv.literature.bibliography import build_index
from vsv.literature.catalog import CitedPaper
from vsv.literature.index import read_index
from vsv.literature.records import Paper
from vsv.literature.run import download_xml, make_index
from vsv.literature.xmlbody import BODY_MIN, body_characters

_BODY = "<article><body><p>" + ("word " * 500) + "</p></body></article>"
_FRONT = (
    "<pmc-articleset><article><front><article-title>Title</article-title></front>"
    "</article></pmc-articleset>"
)
_LICENSED = (
    "<article><front><article-meta><license>Text mining is allowed.</license>"
    "</article-meta></front><body><p>" + ("word " * 500) + "</p></body></article>"
)


def test_body_characters_accepts_an_article_and_rejects_other_xml() -> None:
    count = body_characters(_BODY)
    assert count is not None
    assert count >= BODY_MIN
    assert body_characters(_FRONT) == 0
    assert body_characters('{"status": 500}') is None
    assert body_characters("<html><body>words</body></html>") is None


def test_an_index_keeps_pubmed_ids_and_confirms_pmcids(tmp_path: Path) -> None:
    def fetch_json(url: str) -> Mapping[str, object]:
        if "/references" in url and "page=1" in url and "15121047" in url:
            return {
                "hitCount": 1001,
                "referenceList": {"reference": [{"source": "MED", "id": "111"}]},
            }
        if "/references" in url and "page=2" in url and "15121047" in url:
            return {
                "hitCount": 1001,
                "referenceList": {"reference": [{"source": "MED", "id": "444"}]},
            }
        if "/references" in url and "PMC8470541" in url:
            return {
                "hitCount": 2,
                "referenceList": {
                    "reference": [
                        {"source": "MED", "id": "222"},
                        {"source": "MED", "id": "333"},
                        {"source": "DOI", "id": "10.1/skip"},
                    ]
                },
            }
        if "/references" in url:
            return {"hitCount": 0, "referenceList": {"reference": []}}
        if "idconv" in url:
            return {
                "records": [
                    {"requested-id": "111", "status": "error", "errmsg": "not in PMC"},
                    {"requested-id": "222", "pmcid": "PMC222"},
                    {"requested-id": "333", "pmcid": "333"},
                    {"requested-id": "444", "status": "error"},
                ]
            }
        raise AssertionError(url)

    def fetch_text(url: str) -> str:
        assert "db=pubmed" in url
        return """
            <PubmedArticleSet>
              <PubmedArticle>
                <MedlineCitation>
                  <PMID>111</PMID>
                  <Article>
                    <ArticleTitle>Not in PMC</ArticleTitle>
                    <Journal>
                      <JournalIssue>
                        <PubDate><MedlineDate>2004 Jan</MedlineDate></PubDate>
                      </JournalIssue>
                    </Journal>
                  </Article>
                </MedlineCitation>
                <PubmedData><ArticleIdList>
                  <ArticleId IdType="doi">10.1000/one</ArticleId>
                </ArticleIdList></PubmedData>
              </PubmedArticle>
              <PubmedArticle>
                <MedlineCitation>
                  <PMID>333</PMID>
                  <Article>
                    <ArticleTitle>Outside the corpus</ArticleTitle>
                    <Journal><JournalIssue><PubDate><Year>2019</Year></PubDate></JournalIssue></Journal>
                  </Article>
                </MedlineCitation>
                <PubmedData><ArticleIdList>
                  <ArticleId IdType="doi">10.1000/three</ArticleId>
                </ArticleIdList></PubmedData>
              </PubmedArticle>
              <PubmedArticle>
                <MedlineCitation>
                  <PMID>444</PMID>
                  <Article>
                    <ArticleTitle>Second page</ArticleTitle>
                    <Journal><JournalIssue><PubDate><Year>2001</Year></PubDate></JournalIssue></Journal>
                  </Article>
                </MedlineCitation>
              </PubmedArticle>
            </PubmedArticleSet>
            """

    known = _paper("222", title="In the corpus", pmcid="PMC-ignored", open_access=True)
    index = make_index(
        tmp_path,
        {"222": known},
        fetch_json,
        fetch_text,
        "person@example.edu",
        "2026-09-26",
    )
    by_pmid = {paper.pmid: paper for paper in index.papers}
    assert list(by_pmid) == ["111", "222", "333", "444"]
    assert by_pmid["111"].pmcid == ""
    assert by_pmid["111"].title == "Not in PMC"
    assert by_pmid["111"].year == "2004"
    assert by_pmid["111"].doi == "10.1000/one"
    assert by_pmid["111"].in_corpus is False
    assert by_pmid["111"].cited_by == ("15121047",)
    assert by_pmid["222"].pmcid == "PMC222"
    assert by_pmid["222"].title == "In the corpus"
    assert by_pmid["222"].in_corpus is True
    assert by_pmid["222"].open_access is True
    assert by_pmid["222"].cited_by == ("34578125",)
    assert by_pmid["333"].pmcid == "PMC333"
    assert by_pmid["333"].title == "Outside the corpus"
    assert by_pmid["444"].title == "Second page"
    loaded = read_index(tmp_path)
    assert loaded.papers == index.papers
    assert '"label": "Liu 2021"' in (tmp_path / "index.json").read_text(encoding="utf-8")


def test_index_requires_a_contact_email() -> None:
    with pytest.raises(InputError, match="VSV_NCBI_EMAIL"):
        build_index(lambda _url: {}, lambda _url: "", {}, "", "2026-09-26")


def test_a_cached_body_is_not_downloaded_again(tmp_path: Path) -> None:
    fulltext = tmp_path / "fulltext"
    fulltext.mkdir()
    path = fulltext / "PMC222.xml"
    path.write_text(_LICENSED, encoding="utf-8")
    calls: list[str] = []

    def fetch_text(url: str) -> str:
        calls.append(url)
        raise AssertionError(url)

    row = acquire(_cited("222", "PMC222"), fulltext, fetch_text, "person@example.edu", "2026-09-26")
    assert calls == []
    assert row.status == "xml-body"
    assert row.source == "cache"
    assert row.license == "Text mining is allowed."
    assert row.sha256
    assert row.nbytes == path.stat().st_size


def test_europepmc_xml_is_saved_when_the_body_is_long_enough(tmp_path: Path) -> None:
    fulltext = tmp_path / "fulltext"

    def fetch_text(url: str) -> str:
        assert "europepmc" in url
        return _BODY

    row = acquire(_cited("222", "PMC222"), fulltext, fetch_text, "person@example.edu", "2026-09-26")
    assert row.status == "xml-body"
    assert "europepmc" in row.source
    saved = (fulltext / "PMC222.xml").read_text(encoding="utf-8")
    assert saved == _BODY


def test_efetch_supplies_the_body_when_europepmc_fails(tmp_path: Path) -> None:
    fulltext = tmp_path / "fulltext"
    seen: list[str] = []

    def fetch_text(url: str) -> str:
        seen.append(url)
        if "europepmc" in url:
            raise InputError("HTTP 500")
        if "efetch" in url:
            return _LICENSED
        raise AssertionError(url)

    row = acquire(_cited("222", "PMC222"), fulltext, fetch_text, "a@b.edu", "2026-09-26")
    assert row.status == "xml-body"
    assert "db=pmc" in row.source
    assert row.license == "Text mining is allowed."
    assert not any("oai" in url for url in seen)


def test_front_matter_is_not_saved_and_does_not_try_oai(tmp_path: Path) -> None:
    fulltext = tmp_path / "fulltext"
    seen: list[str] = []

    def fetch_text(url: str) -> str:
        seen.append(url)
        if "europepmc" in url:
            raise InputError("no open xml")
        return _FRONT

    row = acquire(_cited("222", "PMC222"), fulltext, fetch_text, "a@b.edu", "2026-09-26")
    assert row.status == "front-matter-only"
    assert row.sha256 == ""
    assert list(fulltext.glob("*.xml")) == []
    assert not any("oai" in url for url in seen)


def test_oai_is_used_when_efetch_is_not_xml(tmp_path: Path) -> None:
    fulltext = tmp_path / "fulltext"

    def fetch_text(url: str) -> str:
        if "europepmc" in url or "efetch" in url:
            return '{"error": "not xml"}'
        assert "oai" in url
        return _BODY

    row = acquire(_cited("222", "PMC222"), fulltext, fetch_text, "a@b.edu", "2026-09-26")
    assert row.status == "xml-body"
    assert "oai" in row.source
    assert (fulltext / "PMC222.xml").is_file()


def test_a_pmid_with_no_pmcid_is_left_for_the_browser(tmp_path: Path) -> None:
    row = acquire(_cited("111", ""), tmp_path, lambda _url: _BODY, "a@b.edu", "2026-09-26")
    assert row.status == "no-pmc"
    assert list(tmp_path.glob("*.xml")) == []


def test_fetch_writes_a_manifest_and_counts_statuses(tmp_path: Path) -> None:
    census = tmp_path / "census"
    fulltext = tmp_path / "fulltext"

    def fetch_json(url: str) -> Mapping[str, object]:
        if "/references" in url and "PMC8470541" in url:
            return {
                "hitCount": 2,
                "referenceList": {
                    "reference": [
                        {"source": "MED", "id": "222"},
                        {"source": "MED", "id": "111"},
                    ]
                },
            }
        if "/references" in url:
            return {"hitCount": 0, "referenceList": {"reference": []}}
        if "idconv" in url:
            return {
                "records": [
                    {"requested-id": "222", "pmcid": "PMC222"},
                    {"requested-id": "111", "status": "error"},
                ]
            }
        raise AssertionError(url)

    def fetch_text(url: str) -> str:
        if "db=pubmed" in url:
            return """
                <PubmedArticleSet><PubmedArticle><MedlineCitation>
                  <PMID>111</PMID>
                  <Article><ArticleTitle>Paywalled</ArticleTitle>
                    <Journal><JournalIssue><PubDate><Year>2004</Year></PubDate></JournalIssue></Journal>
                  </Article>
                </MedlineCitation></PubmedArticle></PubmedArticleSet>
                """
        if "europepmc" in url:
            raise InputError("no xml")
        return _BODY

    make_index(census, {}, fetch_json, fetch_text, "a@b.edu", "2026-09-26")
    counts = download_xml(census, fulltext, fetch_text, "a@b.edu", "2026-09-26")
    assert counts["xml-body"] == 1
    assert counts["no-pmc"] == 1
    manifest = (census / "manifest.jsonl").read_text(encoding="utf-8")
    assert "Paywalled" in manifest
    assert '"status": "xml-body"' in manifest
    assert (fulltext / "PMC222.xml").is_file()


def test_census_help() -> None:
    assert main(["census", "--help"]) == 0


def _paper(pmid: str, title: str, pmcid: str, open_access: bool) -> Paper:
    return Paper(
        pmid,
        "10.1000/corpus",
        pmcid,
        "2011",
        (),
        (),
        title,
        "",
        open_access,
        "pubmed",
        "2026-09-21",
        True,
    )


def _cited(pmid: str, pmcid: str) -> CitedPaper:
    return CitedPaper(pmid, pmcid, "10.1000/x", "Title", "2011", False, True, ("34578125",))
