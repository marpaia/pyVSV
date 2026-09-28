"""Catalog identity, document formats, and workspace commands. No network."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

import pytest

from vsv.cli import main
from vsv.literature.catalog import CensusIndex, CitedPaper, ManifestRow
from vsv.literature.documents import download_xml, fetch_collection
from vsv.literature.index import read_manifest, write_index, write_manifest
from vsv.literature.papers import digest_pmids, load_ledger, write_ledger, write_papers
from vsv.literature.records import PUBMED, REVIEW_REFERENCES, Ledger, Paper
from vsv.literature.store import (
    census_dir,
    corpus_dir,
    document_from_row,
    fulltext_dir,
    load_catalog,
)
from vsv.literature.text import parse_collection

_BODY = "<article><body><p>" + ("word " * 500) + "</p></body></article>"
_TEXT = "Results\n\nMice were protected.\n"


def test_a_paper_outside_pubmed_keeps_its_record_and_its_citation(tmp_path: Path) -> None:
    _write_overlap(tmp_path)
    catalog = load_catalog(tmp_path)
    by_pmid = {paper.pmid: paper for paper in catalog.papers}
    assert by_pmid["111"].title == "Outside the pubmed snapshot"
    assert by_pmid["111"].year == "2024"
    assert catalog.collections_of("111") == (REVIEW_REFERENCES,)
    assert by_pmid["222"].title == "In the corpus"
    assert by_pmid["222"].pmcid == "PMC-corpus"
    assert catalog.collections_of("222") == (PUBMED, REVIEW_REFERENCES)
    assert catalog.collections_of("333") == (PUBMED,)
    cited = {(item.source, item.target) for item in catalog.citations}
    assert ("15121047", "111") in cited
    assert ("34578125", "222") in cited
    review = catalog.named(REVIEW_REFERENCES)
    pubmed = catalog.named(PUBMED)
    assert review is not None and pubmed is not None
    assert review.retrieved == "2026-09-26"
    assert pubmed.retrieved == "2026-09-21"
    assert "111" in review.pmids
    assert "111" not in pubmed.pmids


def test_documents_split_format_method_and_availability(tmp_path: Path) -> None:
    _write_manifest(tmp_path)
    catalog = load_catalog(tmp_path)
    by_format = {
        (document.pmid, document.format, document.availability, document.method)
        for document in catalog.documents
    }
    assert ("111", "pdf", "available", "browser") in by_format
    assert ("222", "html", "available", "browser") in by_format
    assert ("444", "", "unavailable", "") in by_format
    assert ("555", "xml", "partial", "download") in by_format
    cached = document_from_row(_row("666", "PMC666", "xml-body", source="cache", sha256="abc"))
    assert cached.format == "xml"
    assert cached.availability == "available"
    assert cached.method == "cache"
    assert cached.location == "corpus/fulltext/PMC666.xml"


def test_one_paper_keeps_an_xml_file_and_an_html_file(tmp_path: Path) -> None:
    _write_overlap(tmp_path)
    xml = fulltext_dir(tmp_path) / "PMC222.xml"
    xml.parent.mkdir(parents=True)
    xml.write_bytes(_BODY.encode())
    html = census_dir(tmp_path) / "acquired" / "222.html"
    html.parent.mkdir(parents=True)
    html.write_bytes(b"<html>publisher</html>")
    digest = hashlib.sha256(xml.read_bytes()).hexdigest()
    write_manifest(
        census_dir(tmp_path),
        (_row("222", "PMC222", "xml-body", sha256=digest, nbytes=xml.stat().st_size),),
    )
    formats = sorted(
        document.format
        for document in load_catalog(tmp_path).documents
        if document.pmid == "222" and document.availability == "available"
    )
    assert formats == ["html", "xml"]


def test_download_xml_keeps_rows_it_does_not_replace(tmp_path: Path) -> None:
    census = _write_manifest(tmp_path)
    original = (census / "manifest.jsonl").read_text(encoding="utf-8")
    fulltext = fulltext_dir(tmp_path)

    def fetch_text(_url: str) -> str:
        return _BODY

    download_xml(census, fulltext, fetch_text, "a@b.edu", "2026-09-26", limit=1)
    assert (census / "manifest.jsonl").read_text(encoding="utf-8") == original

    download_xml(census, fulltext, fetch_text, "a@b.edu", "2026-09-26")
    rows = {row.pmid: row for row in read_manifest(census)}
    assert rows["111"].status == "browser-pdf"
    assert rows["222"].status == "xml-body"
    assert rows["222"].sha256
    assert rows["333"].status == "browser-html"
    assert rows["444"].status == "no-pmc"
    assert rows["555"].status == "front-matter-only"
    assert (fulltext / "PMC222.xml").is_file()

    download_xml(census, fulltext, fetch_text, "a@b.edu", "2026-09-26")
    stable = (census / "manifest.jsonl").read_text(encoding="utf-8")
    download_xml(census, fulltext, fetch_text, "a@b.edu", "2026-09-26")
    assert (census / "manifest.jsonl").read_text(encoding="utf-8") == stable
    kept = {row.pmid: row.status for row in read_manifest(census)}
    assert kept == {
        "111": "browser-pdf",
        "222": "xml-body",
        "333": "browser-html",
        "444": "no-pmc",
        "555": "front-matter-only",
    }


def test_pubmed_fetch_leaves_files_outside_the_run(tmp_path: Path) -> None:
    corpus = corpus_dir(tmp_path)
    write_papers(
        (
            _paper("1", "One", "2020", pmcid="PMC1", open_access=True),
            _paper("2", "Two", "2021", pmcid="PMC2", open_access=True),
        ),
        corpus,
    )
    cache = fulltext_dir(tmp_path)
    cache.mkdir(parents=True)
    (cache / "PMC2.xml").write_text(_BODY, encoding="utf-8")
    calls: list[str] = []

    def fetch_text(url: str) -> str:
        calls.append(url)
        return _BODY

    fetch_collection(
        tmp_path,
        PUBMED,
        fetch_text=fetch_text,
        email="",
        retrieved="2026-09-26",
        limit=1,
    )
    assert calls
    assert "PMC1" in calls[0]
    assert (cache / "PMC2.xml").read_text(encoding="utf-8") == _BODY
    calls.clear()
    fetch_collection(
        tmp_path,
        PUBMED,
        fetch_text=fetch_text,
        email="",
        retrieved="2026-09-26",
    )
    assert calls == []
    assert (cache / "PMC1.xml").is_file()


def test_parse_reads_each_format_and_can_be_repeated(tmp_path: Path) -> None:
    _write_overlap(tmp_path)
    xml_path = fulltext_dir(tmp_path) / "PMC333.xml"
    xml_path.parent.mkdir(parents=True)
    xml_path.write_text(_BODY, encoding="utf-8")
    acquired = census_dir(tmp_path) / "acquired"
    acquired.mkdir()
    (acquired / "111.txt").write_text(_TEXT, encoding="utf-8")
    write_manifest(
        census_dir(tmp_path),
        (
            _row("111", "", "browser-html", sha256="html"),
            _row("222", "PMC222", "no-pmc"),
        ),
    )
    assert parse_collection(tmp_path, REVIEW_REFERENCES) == 1
    assert parse_collection(tmp_path, PUBMED) == 1
    html_reading = census_dir(tmp_path) / "sections" / "111.json"
    xml_reading = census_dir(tmp_path) / "sections" / "333.json"
    first_html = html_reading.read_text(encoding="utf-8")
    first_xml = xml_reading.read_text(encoding="utf-8")
    assert "Mice were protected." in first_html
    assert "word word" in first_xml
    assert parse_collection(tmp_path, REVIEW_REFERENCES) == 1
    assert html_reading.read_text(encoding="utf-8") == first_html
    assert xml_reading.read_text(encoding="utf-8") == first_xml
    catalog = load_catalog(tmp_path)
    assert {reading.pmid for reading in catalog.readings} == {"111", "333"}


def test_literature_commands_use_the_workspace(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _write_overlap(tmp_path)
    assert main(["literature"]) == 0
    assert "Usage: vsv literature" in capsys.readouterr().out
    assert main(["literature", "sync"]) == 0
    assert "review-references" in capsys.readouterr().out
    assert main(["literature", "--root", str(tmp_path), "collections"]) == 0
    assert "the pubmed query" in capsys.readouterr().out
    assert (
        main(
            [
                "literature",
                "--root",
                str(tmp_path),
                "papers",
                "--collection",
                "review-references",
                "--year",
                "2024",
                "--limit",
                "0",
            ]
        )
        == 0
    )
    listed = capsys.readouterr().out
    assert "111" in listed
    assert "Outside the pubmed snapshot" in listed
    assert "In the corpus" not in listed
    assert main(["literature", "--json", "--root", str(tmp_path), "status"]) == 0
    status = _captured_json(capsys)
    assert status["papers"] == 3
    assert status["in_both"] == 1
    assert main(["literature", "--root", str(tmp_path), "show", "222"]) == 0
    shown = capsys.readouterr().out
    assert "In the corpus" in shown
    assert "pubmed, review-references" in shown
    assert main(["literature", "--root", str(tmp_path), "show", "999"]) == 1
    assert main(["literature", "--root", str(tmp_path), "check"]) == 0
    assert capsys.readouterr().out.strip() == "ok"
    assert main(["corpus", "--help"]) == 0
    assert main(["census", "--help"]) == 0
    assert main(["papers", "--help"]) == 0


def test_check_fails_when_an_available_file_does_not_match(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _write_overlap(tmp_path)
    raw = b"not the hash we recorded"
    path = census_dir(tmp_path) / "acquired" / "111.pdf"
    path.parent.mkdir(parents=True)
    path.write_bytes(raw)
    digest = hashlib.sha256(raw).hexdigest()
    write_manifest(
        census_dir(tmp_path),
        (_row("111", "", "browser-pdf", sha256="0" * 64, nbytes=len(raw)),),
    )
    assert main(["literature", "--root", str(tmp_path), "check"]) == 1
    assert "checksum mismatch" in capsys.readouterr().out
    write_manifest(
        census_dir(tmp_path),
        (_row("111", "", "browser-pdf", sha256=digest, nbytes=len(raw)),),
    )
    assert main(["literature", "--root", str(tmp_path), "check"]) == 0
    ledger = load_ledger(corpus_dir(tmp_path))
    assert ledger.pmids == ("222", "333")


def _write_overlap(root: Path) -> None:
    corpus = corpus_dir(root)
    write_papers(
        (
            _paper("222", "In the corpus", "2011", pmcid="PMC-corpus"),
            _paper("333", "Only pubmed", "2020", pmcid="PMC333"),
        ),
        corpus,
    )
    pmids = ("222", "333")
    write_ledger(
        Ledger("the pubmed query", "residual", "2026-09-21", 2, 0, pmids, (), digest_pmids(pmids)),
        corpus,
    )
    write_index(
        census_dir(root),
        _index(
            (
                CitedPaper(
                    "111",
                    "",
                    "10.1000/outside",
                    "Outside the pubmed snapshot",
                    "2024",
                    False,
                    True,
                    ("15121047",),
                ),
                CitedPaper(
                    "222",
                    "PMC222",
                    "10.1000/inside",
                    "Index title is not the paper",
                    "2011",
                    True,
                    False,
                    ("34578125",),
                ),
            )
        ),
    )


def _write_manifest(root: Path) -> Path:
    census = census_dir(root)
    write_index(
        census,
        _index(
            (
                CitedPaper("111", "", "10.1/a", "Pdf", "2001", False, False, ("15121047",)),
                CitedPaper("222", "PMC222", "10.1/b", "Html", "2002", False, True, ("34578125",)),
            )
        ),
    )
    write_manifest(
        census,
        (
            _row("111", "", "browser-pdf", sha256="a" * 64, nbytes=10),
            _row("222", "PMC222", "browser-html", sha256="b" * 64, nbytes=20),
            _row("333", "", "browser-html", sha256="c" * 64, nbytes=30),
            _row("444", "", "no-pmc"),
            _row("555", "PMC555", "front-matter-only"),
        ),
    )
    return census


def _index(papers: tuple[CitedPaper, ...]) -> CensusIndex:
    return CensusIndex("2026-09-26", papers)


def _paper(
    pmid: str,
    title: str,
    year: str,
    *,
    pmcid: str = "",
    open_access: bool = False,
) -> Paper:
    return Paper(
        pmid,
        "",
        pmcid,
        year,
        (),
        (),
        title,
        "",
        open_access,
        f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        "2026-09-21",
        True,
    )


def _row(
    pmid: str,
    pmcid: str,
    status: str,
    *,
    source: str = "https://example.test/article",
    sha256: str = "",
    nbytes: int = 0,
) -> ManifestRow:
    return ManifestRow(
        pmid,
        pmcid,
        "10.1000/x",
        f"Title {pmid}",
        "2001",
        source,
        "2026-09-26",
        sha256,
        nbytes,
        0,
        "",
        status,
        "",
        ("15121047",),
    )


def _captured_json(capsys: pytest.CaptureFixture[str]) -> Mapping[str, object]:
    out = capsys.readouterr().out
    parsed = json.loads(out)
    assert isinstance(parsed, dict)
    return parsed
