"""Command-line behavior, run in process."""

import json
from pathlib import Path

import pytest
import sbol3

from tests.synthetic import SyntheticGenome, build_genome
from vsv import sbol
from vsv.cli import main
from vsv.formats import OUTPUT_FORMATS
from vsv.genes import Gene


@pytest.fixture
def pair(tmp_path: Path, synthetic: SyntheticGenome) -> list[Path]:
    codons = {name: "GCT" * (10 + 5 * i) for i, name in enumerate(synthetic.genes)}
    codons[Gene.GLYCOPROTEIN] = "GCTGCA" + "GCT" * 12
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
    start, end = synthetic.cds[Gene.NUCLEOPROTEIN]
    assert any(
        row[1] == "CDS" and "n CDS" in line and str(start + 1) in line and str(end) in line
        for line, row in zip(out.splitlines()[1:], rows, strict=True)
    )
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
    assert [line.split()[2] for line in genes[2:]] == [gene.short for gene in Gene]


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
    assert [row[1] for row in rows[1:]] == [gene.short for gene in Gene]
    start, end = synthetic.genes[Gene.NUCLEOPROTEIN]
    cds_start, cds_end = synthetic.cds[Gene.NUCLEOPROTEIN]
    assert rows[1] == [
        "SYN0001.1",
        "n",
        str(start + 1),
        str(end),
        str(end - start),
        str(cds_start + 1),
        str(cds_end),
        str(cds_end - cds_start),
        str(len(synthetic.proteins[Gene.NUCLEOPROTEIN])),
    ]


def test_genes_leaves_missing_spans_blank(
    capsys: pytest.CaptureFixture[str], synthetic_genbank: Path
) -> None:
    status, out, _ = run(capsys, "genes", str(synthetic_genbank))
    assert status == 0
    # The record annotates coding sequences only, so the gene columns are blank.
    assert out.splitlines()[1].split()[:2] == ["SYN0001", "n"]
    assert len(out.splitlines()[1].split()) == 6


def test_extract_writes_cds_fasta_by_default(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path, synthetic: SyntheticGenome
) -> None:
    status, out, err = run(capsys, "extract", "-g", "p", "-g", "n", str(synthetic_fasta))
    assert (status, err) == (0, "")
    headers = [line for line in out.splitlines() if line.startswith(">")]
    assert headers == [">SYN0001.1 p cds", ">SYN0001.1 n cds"]
    start, end = synthetic.cds[Gene.PHOSPHOPROTEIN]
    assert out.splitlines()[1] == synthetic.sequence[start:end]


def test_extract_writes_protein_fasta(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path, synthetic: SyntheticGenome
) -> None:
    status, out, _ = run(capsys, "extract", "--kind", "protein", "-g", "m", str(synthetic_fasta))
    assert status == 0
    assert out.splitlines() == [
        ">SYN0001.1 m protein",
        synthetic.proteins[Gene.MATRIX_PROTEIN],
    ]


def test_extract_reports_an_unannotated_product(
    capsys: pytest.CaptureFixture[str], synthetic_genbank: Path
) -> None:
    status, out, err = run(capsys, "extract", "--kind", "gene", "-g", "n", str(synthetic_genbank))
    assert (status, out) == (1, "")
    assert "SYN0001: no gene annotated for gene n." in err


def test_incomplete_genomes_warn_on_stderr(
    capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    short = tmp_path / "short.fasta"
    short.write_text(">SHORT\nACGT\n")
    status, _, err = run(capsys, "genes", str(short))
    assert status == 0
    assert ("vsv: warning: SHORT: genes not identified: n, p, m, g, l") in err


def test_designs_can_be_filtered_by_antigen_class(capsys: pytest.CaptureFixture[str]) -> None:
    status, out, err = run(capsys, "designs", "--antigen-class", "henipavirus-glycoprotein")
    assert (status, err) == (0, "")
    assert "rvsv_niv_g" in out
    assert "rvsv_niv_f" not in out
    status, out, err = run(capsys, "designs")
    assert (status, err) == (0, "")
    assert "gaps" in out.splitlines()[0]


def test_show_separates_lassa_implementations(capsys: pytest.CaptureFixture[str]) -> None:
    status, out, err = run(capsys, "show", "rvsv_lasv_gpc")
    assert (status, err) == (0, "")
    assert "derivation   replace g" in out
    assert "  3      g        arenavirus-gpc" in out
    assert "preclinical  garbutt-2004" in out
    assert "  no observations" in out
    assert "clinical  malkin-2025" in out
    assert "phenotype       host   assay           outcome" in out
    assert "immunogenicity  human  immunogenicity  reported" in out
    assert "reactogenicity" in out
    assert "hearing loss" in out
    assert out.count("https://doi.org/10.1056/NEJMoa2501073") == 1
    assert "Zaire" not in out


def test_show_records_the_nipah_locus_as_unknown(capsys: pytest.CaptureFixture[str]) -> None:
    status, out, err = run(capsys, "show", "data/designs/rvsv_niv_g")
    assert (status, err) == (0, "")
    assert "parent       rvsv_dg_filovirus_gp" in out
    assert not any(line.startswith("derivation") for line in out.splitlines())
    assert "  unknown  niv_g    henipavirus-glycoprotein" in out
    assert "locus unknown for niv_g" in out


def test_check_accepts_the_canonical_chassis(capsys: pytest.CaptureFixture[str]) -> None:
    status, out, err = run(capsys, "check", "data/designs/vsv_indiana")
    assert (status, err) == (0, "")
    assert out.strip() == "vsv_indiana: ok"


def test_check_reports_a_knowledge_gap_without_failing(capsys: pytest.CaptureFixture[str]) -> None:
    status, out, err = run(capsys, "check", "data/designs/rvsv_niv_g")
    assert (status, err) == (0, "")
    assert out.splitlines()[0].split()[0] == "gap"
    assert "locus-unknown" in out


def test_check_fails_a_structural_violation(
    capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    design = tmp_path / "bad"
    design.mkdir()
    (design / "design.json").write_text(
        json.dumps(
            {
                "id": "bad",
                "name": "bad",
                "backbone": "VSV",
                "replication": "replication-competent",
                "modules": [],
            }
        )
    )
    status, out, err = run(capsys, "check", str(design))
    assert (status, err) == (1, "")
    assert "replication" in out


def test_extract_rejects_unknown_genes(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path
) -> None:
    with pytest.raises(SystemExit) as raised:
        main(["extract", "-g", "X", str(synthetic_fasta)])
    assert raised.value.code == 2
    assert "invalid gene" in capsys.readouterr().err
