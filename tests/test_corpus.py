"""Corpus fetch, classification, extraction, and promotion stay offline."""

import urllib.parse
from pathlib import Path

import pytest

from vsv.cli import main
from vsv.design.library import load_design
from vsv.design.literature import clusters_from, write_design
from vsv.literature.analysis.classify import classify
from vsv.literature.analysis.evaluate import audit_messages, stratified_pmids
from vsv.literature.analysis.extract import describe
from vsv.literature.analysis.records import PaperKind
from vsv.literature.analysis.workflow import extract_corpus
from vsv.literature.checks import coverage_report
from vsv.literature.discovery.pubmed import europepmc_records, ingest_corpus
from vsv.literature.documents.jats import methods_text
from vsv.literature.documents.text import clean_text
from vsv.literature.records import Paper
from vsv.literature.search import mentions
from vsv.literature.storage.corpus import load_extractions, load_ledger, load_papers, write_papers


def paper(
    *,
    pmid: str = "1",
    year: str = "2020",
    pub_types: tuple[str, ...] = (),
    mesh: tuple[str, ...] = (),
    title: str = "",
    abstract: str = "",
) -> Paper:
    return Paper(
        pmid,
        "",
        "",
        year,
        pub_types,
        mesh,
        title,
        abstract,
        False,
        f"https://europepmc.org/article/MED/{pmid}",
        "2026-09-21",
        True,
    )


def test_classifier_keeps_vector_reviews_and_pseudotypes() -> None:
    vector = classify(
        paper(
            title="Recombinant VSV expressing Ebola virus glycoprotein",
            pub_types=("Review", "Journal Article"),
        )
    )
    assert (vector.kind, vector.review, vector.status) == (
        PaperKind.VACCINE_VECTOR,
        True,
        "classified",
    )
    pseudo = classify(paper(title="A lentiviral vector pseudotyped with VSV-G"))
    assert pseudo.kind is PaperKind.PSEUDOTYPE
    blank = classify(paper(title="Notes"))
    assert (blank.kind, blank.status) == (PaperKind.OTHER, "unresolved")
    baculovirus = classify(
        paper(
            title="The VSV nucleocapsid protein expressed by recombinant baculovirus",
            abstract="A vaccine preparation expressing the nucleoprotein was used in an ELISA.",
        )
    )
    assert baculovirus.kind is not PaperKind.VACCINE_VECTOR


def test_extraction_quotes_a_replacement_and_redacts_nucleotides() -> None:
    sentence = (
        "The rVSVΔG-ZEBOV-GP vaccine " + ("ACGT" * 15) + " protected mice from lethal challenge."
    )
    item = describe(paper(title="rVSV Ebola vaccine", abstract=sentence))
    assert item.kind == "vaccine-vector"
    assert item.edit_activity.value == "replace-module"
    assert item.edit_target.value == "Glycoprotein"
    assert item.antigen_slug.value == "filovirus-glycoprotein"
    assert item.locus.value == "unknown"
    assert item.stage.value == "preclinical"
    assert item.claims[0].host == "mouse"
    assert item.claims[0].phenotype == "protection"
    assert "ACGTACGT" not in item.edit_activity.quote
    assert "[sequence omitted]" in item.edit_activity.quote


def test_a_g_l_expression_site_inserts_zika_envelope_after_glycoprotein() -> None:
    item = describe(
        paper(
            title="rVSV Zika vaccine",
            abstract=(
                "Controls were included. rVSV constructs were successfully generated using "
                "the G-L foreign gene expression site. The viruses express ZIKV envelope "
                "protein and are fully "
                "replication competent. Mice were protected from lethal challenge. "
                "It warrants further evaluation in phase I trials in human subjects."
            ),
        )
    )
    assert item.edit_activity.value == "insert-module"
    assert item.edit_target.value == "Glycoprotein"
    assert item.antigen_slug.value == "zika-envelope"
    assert item.replication.value == "replication-competent"
    assert item.stage.value == "preclinical"
    assert item.edit_activity.quote.startswith("rVSV constructs")


def test_a_saved_article_supplies_the_architecture_fields(tmp_path: Path) -> None:
    record = paper(
        pmid="28289159",
        title="rVSV Zika vaccine",
        abstract="Recombinant VSV expressing a Zika protein was tested in mice.",
    )
    corpus = tmp_path / "corpus"
    write_papers([record], corpus)
    acquired = tmp_path / "census" / "acquired"
    acquired.mkdir(parents=True)
    body = (
        "rVSV constructs were successfully generated using the G-L foreign gene "
        "expression site.\n\n"
        "The viruses express ZIKV envelope protein and are fully replication competent.\n\n"
    )
    (acquired / "28289159.txt").write_text(body + ("Noted. " * 400), encoding="utf-8")
    report = extract_corpus(corpus)
    stored = load_extractions(corpus)
    assert report.fulltext == 1
    assert stored[0].edit_activity.value == "insert-module"
    assert stored[0].edit_target.value == "Glycoprotein"
    assert stored[0].antigen_slug.value == "zika-envelope"
    assert stored[0].text_source == "fulltext"
    assert stored[0].edit_activity.quote.startswith("rVSV constructs")


def test_species_glycoproteins_keep_their_own_classes() -> None:
    sudan = describe(
        paper(
            pmid="2",
            abstract=(
                "Recombinant VSV expressing Sudan virus glycoprotein replaced the glycoprotein."
            ),
        )
    )
    ebola = describe(
        paper(
            pmid="3",
            abstract=(
                "Recombinant VSV expressing Ebola virus glycoprotein replaced the glycoprotein."
            ),
        )
    )
    assert sudan.antigen_slug.value == "sudan-ebolavirus-glycoprotein"
    assert ebola.antigen_slug.value == "filovirus-glycoprotein"


def test_conflicting_backbones_stay_unknown() -> None:
    item = describe(
        paper(
            title="Recombinant VSV vector",
            abstract=(
                "The Indiana construct is replication-competent. "
                "The New Jersey construct is replication-competent."
            ),
        )
    )
    assert item.backbone.value == "unknown"
    assert "backbone" in item.conflicts
    assert item.replication.value == "replication-competent"


def test_clusters_keep_stated_loci_apart() -> None:
    shared = (
        "Recombinant VSV expressing Ebola virus glycoprotein replaced the glycoprotein. "
        "The construct is replication-competent."
    )
    placed = shared + " The insert is at gene position 1."
    first = describe(paper(pmid="10", abstract=shared))
    second = describe(paper(pmid="11", abstract=shared))
    third = describe(paper(pmid="12", abstract=placed))
    clusters = clusters_from([first, second, third], curated=[])
    assert len(clusters) == 2
    sizes = sorted(len(cluster.pmids) for cluster in clusters)
    assert sizes == [1, 2]
    unknown = next(cluster for cluster in clusters if cluster.locus == "unknown")
    assert unknown.pmids == ("10", "11")
    assert "locus-unknown" in unknown.gaps
    assert "parent-unknown" in unknown.gaps


def test_an_alias_matches_ignoring_case() -> None:
    assert mentions("the rVSV-ZEBOV vaccine", ("rVSV-ZEBOV",))
    assert not mentions("vesicular stomatitis virus", ("Ervebo",))


def test_an_unstated_backbone_joins_one_stated_backbone() -> None:
    shared = "Recombinant VSV expressing Ebola virus glycoprotein replaced the glycoprotein. "
    indiana = describe(
        paper(pmid="1", abstract=shared + "The Indiana construct is replication-competent.")
    )
    unnamed = describe(paper(pmid="2", abstract=shared + "The construct is replication-competent."))
    jersey = describe(
        paper(pmid="3", abstract=shared + "The New Jersey construct is replication-competent.")
    )
    (cluster,) = clusters_from([indiana, unnamed], curated=[])
    assert cluster.backbone == "indiana"
    assert "backbone-unknown" in cluster.gaps
    apart = clusters_from([indiana, jersey], curated=[])
    assert {item.backbone for item in apart} == {"indiana", "new-jersey"}


def test_promote_writes_a_loadable_design_without_sequences(tmp_path: Path) -> None:
    item = describe(
        paper(
            pmid="15121047",
            abstract=(
                "Recombinant VSV expressing Ebola virus glycoprotein replaced the glycoprotein "
                "and is replication-competent. Mice were protected from lethal challenge."
            ),
        )
    )
    (cluster,) = clusters_from([item], curated=[])
    dest = tmp_path / cluster.id
    write_design(cluster, [item], [paper(pmid="15121047", abstract=item.edit_activity.quote)], dest)
    loaded = load_design(dest)
    assert loaded.architecture.edit is not None
    assert loaded.architecture.antigen_classes == ("filovirus-glycoprotein",)
    blob = (dest / "design.ttl").read_text()
    assert "sequence" not in blob
    assert "ACGT" not in blob
    assert loaded.observations[0].host.value == "mouse"


def test_europepmc_retries_an_empty_batch() -> None:
    def fetch(url: str, body: bytes | None = None) -> dict[str, object]:
        del body
        query = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["query"][0]
        if " OR " in query:
            return {"resultList": {"result": []}, "nextCursorMark": "*"}
        pmid = query.split("EXT_ID:", 1)[1].split(" ", 1)[0].rstrip(")")
        return {"resultList": {"result": [{"pmid": pmid}]}, "nextCursorMark": "*"}

    records = europepmc_records(["10", "20"], fetch, batch_size=2)
    assert {str(record["pmid"]) for record in records} == {"10", "20"}


def test_methods_section_is_taken_from_jats() -> None:
    xml = """
    <article>
      <body>
        <sec sec-type="intro"><title>Introduction</title><p>Background.</p></sec>
        <sec sec-type="methods"><title>Methods</title><p>Inserted after the matrix gene.</p></sec>
      </body>
    </article>
    """
    assert "matrix gene" in methods_text(xml)
    assert "Background" not in methods_text(xml)


def test_ingest_stores_every_pmid_and_records_misses(tmp_path: Path) -> None:
    def fetch(url: str, body: bytes | None = None) -> dict[str, object]:
        del body
        parsed = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
        if "eutils.ncbi" in url:
            term = parsed["term"][0]
            if " NOT " in term:
                return {"esearchresult": {"count": "4", "idlist": []}}
            start = int(parsed["retstart"][0])
            size = int(parsed["retmax"][0])
            page = ["10", "20", "30"][start : start + size]
            return {"esearchresult": {"count": "3", "idlist": page}}
        cursor = parsed["cursorMark"][0]
        query = parsed["query"][0]
        if cursor != "*" or "EXT_ID:10" not in query:
            return {"resultList": {"result": []}, "nextCursorMark": cursor}
        return {
            "resultList": {
                "result": [
                    {
                        "pmid": "10",
                        "doi": "10.1000/vsv",
                        "pmcid": "PMC10",
                        "pubYear": 2004,
                        "title": "Recombinant VSV ebola vaccine",
                        "abstractText": "ACGT" * 15,
                        "isOpenAccess": "Y",
                        "pubTypeList": {"pubType": ["Journal Article"]},
                        "meshHeadingList": {
                            "meshHeading": [
                                {"descriptorName": "Vesicular stomatitis Indiana virus"}
                            ]
                        },
                    },
                    {
                        "pmid": "20",
                        "doi": "",
                        "title": "Field notes",
                        "abstractText": "Livestock outbreak.",
                        "isOpenAccess": "N",
                        "pubTypeList": {"pubType": "Journal Article"},
                    },
                ]
            },
            "nextCursorMark": "*",
        }

    report = ingest_corpus(tmp_path, fetch=fetch, page_size=2, batch_size=2)
    assert (report.pubmed_count, report.stored, report.new, report.europepmc_misses) == (3, 3, 3, 1)
    assert report.residual_count == 4
    papers = {item.pmid: item for item in load_papers(tmp_path)}
    assert papers["30"].in_europepmc is False
    assert "[sequence omitted]" in papers["10"].abstract
    assert "ACGTACGT" not in papers["10"].abstract
    assert papers["10"].mesh == ("Vesicular stomatitis Indiana virus",)
    ledger = load_ledger(tmp_path)
    assert ledger.pmids == ("10", "20", "30")
    text, ok = coverage_report(tmp_path)
    assert ok
    assert "extractions  absent" in text
    again = ingest_corpus(tmp_path, fetch=fetch, page_size=2, batch_size=2)
    assert (again.new, again.removed) == (0, 0)


def test_audit_reports_a_missed_vector() -> None:
    item = describe(paper(pmid="5", title="Notes on cattle"))
    messages = audit_messages(
        [item],
        [{"pmid": "5", "kind": "vaccine-vector", "vector": True}],
    )
    assert any("labeled vaccine-vector" in message for message in messages)
    assert any("labeled vector" in message for message in messages)


def test_stratified_sample_is_stable() -> None:
    papers = [paper(pmid=str(number), year=str(1970 + number)) for number in range(1, 30)]
    assert stratified_pmids(papers, size=6, seed=3) == stratified_pmids(papers, size=6, seed=3)
    assert len(stratified_pmids(papers, size=6, seed=3)) == 6


def test_corpus_command_fails_when_the_ledger_is_missing(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    status = main(["corpus", "--dest", str(tmp_path)])
    captured = capsys.readouterr()
    assert status == 1
    assert "ledger" in captured.err


def test_ingest_all_rejects_a_date_window(capsys: pytest.CaptureFixture[str]) -> None:
    status = main(["ingest", "--all", "--since", "2020-01-01"])
    captured = capsys.readouterr()
    assert status == 1
    assert "frozen" in captured.err


def test_clean_text_drops_markup() -> None:
    assert clean_text("<p>Indiana</p>") == "Indiana"
