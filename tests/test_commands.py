"""Command-line behavior, run in process."""

from pathlib import Path

import pytest
import sbol3

from tests.synthetic import SyntheticGenome, build_genome
from vsv import sbol
from vsv.cli import main
from vsv.formats import OUTPUT_FORMATS


@pytest.fixture
def pair(tmp_path: Path, synthetic: SyntheticGenome) -> list[Path]:
    codons = {name: "GCT" * (10 + 5 * i) for i, name in enumerate(synthetic.genes)}
    codons["G"] = "GCTGCA" + "GCT" * 12
    other = build_genome("SYN0002.1", codons)
    first, second = tmp_path / "a.fasta", tmp_path / "b.fasta"
    first.write_text(synthetic.fasta)
    second.write_text(other.fasta)
    return [first, second]


def run(capsys: pytest.CaptureFixture[str], *argv: str) -> tuple[int, str, str]:
    status = main(list(argv))
    captured = capsys.readouterr()
    return status, captured.out, captured.err


def test_summary_reports_each_record(capsys: pytest.CaptureFixture[str], pair: list[Path]) -> None:
    status, out, err = run(capsys, "summary", *map(str, pair))
    assert (status, err) == (0, "")
    assert out.splitlines()[0].split()[:2] == ["id", "length"]
    assert [line.split()[0] for line in out.splitlines()[1:]] == ["SYN0001.1", "SYN0002.1"]


def test_annotate_table(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path, synthetic: SyntheticGenome
) -> None:
    status, out, err = run(capsys, "annotate", str(synthetic_fasta))
    assert (status, err) == (0, "")
    rows = [line.split() for line in out.splitlines()[1:]]
    start, end = synthetic.cds["N"]
    assert [
        "SYN0001.1",
        "CDS",
        "N",
        str(start + 1),
        str(end),
        str(end - start),
        str(len(synthetic.proteins["N"])),
    ] in rows
    assert [row[1] for row in rows].count("gene") == 5


def test_annotate_gff3(capsys: pytest.CaptureFixture[str], synthetic_fasta: Path) -> None:
    status, out, _ = run(capsys, "annotate", "--format", "gff3", str(synthetic_fasta))
    assert status == 0
    assert out.startswith("##gff-version 3\n")
    assert out.count("\tgene\t") == 5


@pytest.mark.parametrize("output_format", ["turtle", "ntriples", "rdfxml", "jsonld"])
def test_annotate_sbol(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path, output_format: str
) -> None:
    status, out, err = run(capsys, "annotate", "--format", output_format, str(synthetic_fasta))
    assert (status, err) == (0, "")
    document = sbol3.Document()
    document.read_string(out, OUTPUT_FORMATS[output_format])
    assert not document.validate()
    (genome,) = sbol.genomes(document)
    assert len(genome.features(sbol.SO_GENE)) == 5


def test_curated_annotations_are_kept(
    capsys: pytest.CaptureFixture[str], synthetic_genbank: Path
) -> None:
    status, out, _ = run(capsys, "annotate", str(synthetic_genbank))
    kinds = [line.split()[1] for line in out.splitlines()[1:]]
    assert status == 0
    # The converter's whole-genome source feature and the five curated coding sequences;
    # nothing from the heuristic annotator.
    assert kinds == ["region", "CDS", "CDS", "CDS", "CDS", "CDS"]


def test_compare(capsys: pytest.CaptureFixture[str], pair: list[Path]) -> None:
    status, out, err = run(capsys, "compare", "-k", "8", *map(str, pair))
    assert (status, err) == (0, "")
    assert "Genome k-mer similarity (k=8)" in out
    genes = out.split("Protein identity by gene")[1].splitlines()
    assert [line.split()[2] for line in genes[2:]] == ["N", "P", "M", "G", "L"]


def test_compare_accepts_mixed_formats(
    capsys: pytest.CaptureFixture[str], synthetic_genbank: Path, pair: list[Path]
) -> None:
    status, out, _ = run(capsys, "compare", str(synthetic_genbank), str(pair[1]))
    assert status == 0
    assert "SYN0001" in out
    assert "100.00%" in out


def test_compare_needs_two_records(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path
) -> None:
    status, out, err = run(capsys, "compare", str(synthetic_fasta))
    assert (status, out) == (1, "")
    assert "requires at least 2 records" in err


@pytest.mark.parametrize(
    ("name", "text"),
    [("bad.fasta", "ACGT\n"), ("bad.gb", "nope\n"), ("bad.txt", ">A\nAC\n")],
)
def test_bad_input_exits_with_an_error(
    capsys: pytest.CaptureFixture[str], tmp_path: Path, name: str, text: str
) -> None:
    path = tmp_path / name
    path.write_text(text)
    status, out, err = run(capsys, "summary", str(path))
    assert (status, out) == (1, "")
    assert err.startswith(f"vsv: error: {path}:")


def test_no_command_prints_help(capsys: pytest.CaptureFixture[str]) -> None:
    status, out, _ = run(capsys)
    assert status == 0
    assert "usage: vsv" in out


def test_genes_lists_the_five_genes(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path, synthetic: SyntheticGenome
) -> None:
    status, out, err = run(capsys, "genes", str(synthetic_fasta))
    assert (status, err) == (0, "")
    rows = [line.split() for line in out.splitlines()]
    assert rows[0][:3] == ["id", "gene", "start"]
    assert [row[1] for row in rows[1:]] == ["N", "P", "M", "G", "L"]
    start, end = synthetic.genes["N"]
    cds_start, cds_end = synthetic.cds["N"]
    assert rows[1] == [
        "SYN0001.1",
        "N",
        str(start + 1),
        str(end),
        str(end - start),
        str(cds_start + 1),
        str(cds_end),
        str(cds_end - cds_start),
        str(len(synthetic.proteins["N"])),
    ]


def test_genes_leaves_missing_spans_blank(
    capsys: pytest.CaptureFixture[str], synthetic_genbank: Path
) -> None:
    status, out, _ = run(capsys, "genes", str(synthetic_genbank))
    assert status == 0
    # The record annotates coding sequences only, so the gene columns are blank.
    assert out.splitlines()[1].split()[:2] == ["SYN0001", "N"]
    assert len(out.splitlines()[1].split()) == 6


def test_extract_writes_cds_fasta_by_default(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path, synthetic: SyntheticGenome
) -> None:
    status, out, err = run(capsys, "extract", "-g", "P", "-g", "N", str(synthetic_fasta))
    assert (status, err) == (0, "")
    headers = [line for line in out.splitlines() if line.startswith(">")]
    assert headers == [">SYN0001.1 P cds", ">SYN0001.1 N cds"]
    start, end = synthetic.cds["P"]
    assert out.splitlines()[1] == synthetic.sequence[start:end]


def test_extract_writes_protein_fasta(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path, synthetic: SyntheticGenome
) -> None:
    status, out, _ = run(capsys, "extract", "--kind", "protein", "-g", "M", str(synthetic_fasta))
    assert status == 0
    assert out.splitlines() == [">SYN0001.1 M protein", synthetic.proteins["M"]]


def test_extract_reports_an_unannotated_product(
    capsys: pytest.CaptureFixture[str], synthetic_genbank: Path
) -> None:
    status, out, err = run(capsys, "extract", "--kind", "gene", "-g", "N", str(synthetic_genbank))
    assert (status, out) == (1, "")
    assert "SYN0001: no gene annotated for gene N." in err


def test_incomplete_genomes_warn_on_stderr(
    capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    short = tmp_path / "short.fasta"
    short.write_text(">SHORT\nACGT\n")
    status, _, err = run(capsys, "genes", str(short))
    assert status == 0
    assert "vsv: warning: SHORT: genes not identified: N, P, M, G, L" in err


def test_extract_rejects_unknown_genes(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path
) -> None:
    with pytest.raises(SystemExit) as raised:
        main(["extract", "-g", "X", str(synthetic_fasta)])
    assert raised.value.code == 2
    assert "invalid Gene value" in capsys.readouterr().err
