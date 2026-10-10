"""Command-line behavior, run in process."""

from pathlib import Path

import pytest
import sbol3

from tests.synthetic import SyntheticGenome, genbank_text
from vsv import Gene, GlycoproteinForm, Replication, Vaccine, Viral, indiana
from vsv.cli import main
from vsv.sbol import DNA_ENCODING, PROTEIN_ENCODING

RAW = Path(__file__).resolve().parents[1] / "data" / "raw"


def run(capsys: pytest.CaptureFixture[str], *argv: str) -> tuple[int, str, str]:
    status = main(list(argv))
    captured = capsys.readouterr()
    return status, captured.out, captured.err


def _sbol_records(text: str) -> dict[str, tuple[str, str]]:
    """Component name to sequence elements and encoding."""
    sbol3.set_namespace("https://vsv.bio")
    document = sbol3.Document()
    document.read_string(text, "ttl")
    found: dict[str, tuple[str, str]] = {}
    for component in document.objects:
        if not isinstance(component, sbol3.Component):
            continue
        (uri,) = list(component.sequences)
        sequence = document.find(str(uri))
        assert isinstance(sequence, sbol3.Sequence)
        found[str(component.name)] = (str(sequence.elements), str(sequence.encoding))
    return found


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


def test_coding_sequences_keep_gene_spans(
    capsys: pytest.CaptureFixture[str], tmp_path: Path, synthetic: SyntheticGenome
) -> None:
    path = tmp_path / "cds.gb"
    path.write_text(genbank_text(synthetic))
    status, out, err = run(capsys, "annotate", str(path))
    assert (status, err) == (0, "")
    kinds = [line.split()[1] for line in out.splitlines()[1:]]
    assert kinds == ["gene", "CDS"] * 5
    status, genes_out, _ = run(capsys, "genes", str(path))
    assert status == 0
    row = genes_out.splitlines()[1].split()
    start, end = synthetic.cds[Gene.NUCLEOPROTEIN]
    assert row[:2] == ["SYN0001.1", "n"]
    assert row[2:8] == [
        str(start + 1),
        str(end),
        str(end - start),
        str(start + 1),
        str(end),
        str(end - start),
    ]


@pytest.mark.parametrize(
    ("name", "text"),
    [("bad.fasta", "ACGT\n"), ("bad.gb", "nope\n"), ("bad.txt", ">A\nAC\n")],
)
def test_bad_input_exits_with_an_error(
    capsys: pytest.CaptureFixture[str], tmp_path: Path, name: str, text: str
) -> None:
    path = tmp_path / name
    path.write_text(text)
    status, out, err = run(capsys, "genes", str(path))
    assert (status, out) == (1, "")
    assert err.startswith(f"vsv: error: {path}:")


def test_no_command_prints_help(capsys: pytest.CaptureFixture[str]) -> None:
    status, out, _ = run(capsys)
    assert status == 0
    assert "Usage: vsv" in out


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


def test_extract_writes_cds_sbol_by_default(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path, synthetic: SyntheticGenome
) -> None:
    status, out, err = run(capsys, "extract", "-g", "p", "-g", "n", str(synthetic_fasta))
    assert (status, err) == (0, "")
    found = _sbol_records(out)
    start, end = synthetic.cds[Gene.PHOSPHOPROTEIN]
    n_start, n_end = synthetic.cds[Gene.NUCLEOPROTEIN]
    assert found["SYN0001.1 p cds"] == (synthetic.sequence[start:end], DNA_ENCODING)
    assert found["SYN0001.1 n cds"] == (synthetic.sequence[n_start:n_end], DNA_ENCODING)


def test_extract_writes_protein_sbol(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path, synthetic: SyntheticGenome
) -> None:
    status, out, _ = run(capsys, "extract", "--kind", "protein", "-g", "m", str(synthetic_fasta))
    assert status == 0
    found = _sbol_records(out)
    assert found == {
        "SYN0001.1 m protein": (synthetic.proteins[Gene.MATRIX_PROTEIN], PROTEIN_ENCODING)
    }


def test_extract_writes_the_coding_span_as_the_gene(
    capsys: pytest.CaptureFixture[str], tmp_path: Path, synthetic: SyntheticGenome
) -> None:
    path = tmp_path / "cds.gb"
    path.write_text(genbank_text(synthetic))
    status, out, err = run(capsys, "extract", "--kind", "gene", "-g", "n", str(path))
    assert (status, err) == (0, "")
    start, end = synthetic.cds[Gene.NUCLEOPROTEIN]
    elements, encoding = _sbol_records(out)["SYN0001.1 n gene"]
    assert encoding == DNA_ENCODING
    assert elements == synthetic.sequence[start:end]


def test_incomplete_genomes_warn_on_stderr(
    capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    short = tmp_path / "short.fasta"
    short.write_text(">SHORT\nACGT\n")
    status, _, err = run(capsys, "genes", str(short))
    assert status == 0
    assert "vsv: warning: SHORT: genes not identified: n, p, m, g, l" in err


def test_indiana_genes_and_extract_name_five_genes(capsys: pytest.CaptureFixture[str]) -> None:
    fasta, genbank = RAW / "NC_001560.fasta", RAW / "NC_001560.gb"
    if not fasta.exists() or not genbank.exists():
        pytest.skip("NC_001560 is not in data/raw")
    status, out, err = run(capsys, "genes", str(fasta))
    assert (status, err) == (0, "")
    assert [line.split()[1] for line in out.splitlines()[1:]] == ["n", "p", "m", "g", "l"]
    turtle = RAW / "NC_001560.ttl"
    if turtle.exists():
        status, out, err = run(capsys, "genes", str(turtle))
        assert (status, err) == (0, "")
        assert [line.split()[1] for line in out.splitlines()[1:]] == ["n", "p", "m", "g", "l"]
    status, out, err = run(capsys, "extract", str(genbank))
    assert (status, err) == (0, "")
    found = _sbol_records(out)
    assert set(found) == {f"NC_001560.1 {gene.short} cds" for gene in Gene}
    assert {encoding for _, encoding in found.values()} == {DNA_ENCODING}


def test_designs_can_be_filtered(capsys: pytest.CaptureFixture[str]) -> None:
    status, out, err = run(capsys, "designs", "--antigen-class", "zaire-ebolavirus-glycoprotein")
    assert (status, err) == (0, "")
    assert "rvsv_dg_zebov_kikwit_gp" in out
    assert "rvsv_dg_zebov_mayinga_gp" in out
    assert "rvsv_gl_zebov_sgp" in out
    assert "rvsv_dg_sudv_gp" not in out
    status, out, err = run(capsys, "designs", "--antigen-class", "filovirus-glycoprotein")
    assert (status, err) == (0, "")
    assert "vsv_indiana" not in out
    assert out.count("rvsv_") == 7
    status, out, err = run(
        capsys, "designs", "--phenotype", "protection", "--host", "nonhuman-primate"
    )
    assert (status, err) == (0, "")
    assert {line.split()[0] for line in out.splitlines()[1:]} == {
        "rvsv_dg_bdbv_gp",
        "rvsv_dg_marv_gp",
        "rvsv_dg_sudv_gp",
        "rvsv_dg_zebov_mayinga_gp",
    }
    status, out, err = run(capsys, "designs", "--evidence", "human")
    assert (status, err) == (0, "")
    assert [line.split()[0] for line in out.splitlines()[1:]] == ["rvsv_dg_zebov_kikwit_gp"]
    status, out, err = run(capsys, "designs")
    assert (status, err) == (0, "")
    assert out.splitlines()[0].split()[:6] == [
        "id",
        "edit",
        "antigens",
        "replication",
        "evidence",
        "gaps",
    ]


def test_show_separates_kikwit_implementations(capsys: pytest.CaptureFixture[str]) -> None:
    status, out, err = run(capsys, "show", "rvsv_dg_zebov_kikwit_gp")
    assert (status, err) == (0, "")
    assert "derivation   replace g, replication-competent" in out
    assert "  3      g       zaire-ebolavirus-glycoprotein  Kikwit-1995  full-length" in out
    assert "clinical  regules-2017" in out
    assert "clinical  agnandji-2016" in out
    assert "clinical  henao-restrepo-2017" in out
    assert "protection  human  field-efficacy  reported" in out
    assert "The most common adverse events were injection-site pain" in out
    assert out.count("https://doi.org/10.1056/NEJMoa1414216") == 1
    assert (
        "A gene closer to the 3' end is transcribed more, and position 1 attenuates the virus."
    ) in out


def test_precedents_locus_separates_a_leader_insert_from_a_g_replacement(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status, out, err = run(
        capsys, "precedents", "data/designs/rvsv_dg_zebov_kikwit_gp", "--locus", "0"
    )
    assert (status, out, err) == (0, "no precedents\n", "")
    status, out, err = run(
        capsys, "precedents", "data/designs/rvsv_dg_zebov_kikwit_gp", "--locus", "3"
    )
    assert (status, err) == (0, "")
    assert "rvsv_dg_zebov_mayinga_gp" in out
    status, out, err = run(capsys, "precedents", "data/designs/rvsv_gl_zebov_sgp", "--locus", "0")
    assert (status, err) == (0, "")
    assert "rvsv_dg_zebov_mayinga_gp" not in out
    assert "rvsv_dg_zebov_kikwit_gp" not in out
    status, out, err = run(capsys, "precedents", "data/designs/rvsv_gl_zebov_sgp", "--locus", "4")
    assert (status, out, err) == (0, "no precedents\n", "")
    status, out, err = run(capsys, "precedents", "rvsv_dg_zebov_mayinga_gp")
    assert (status, err) == (0, "")
    assert "rvsv_dg_zebov_kikwit_gp" in out
    assert "rvsv_dg_bombali_gp" not in out
    status, out, err = run(capsys, "precedents", "vsv_indiana")
    assert (status, out, err) == (0, "no precedents\n", "")


def test_show_places_soluble_glycoprotein_after_vsv_g(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status, out, err = run(capsys, "show", "data/designs/rvsv_gl_zebov_sgp")
    assert (status, err) == (0, "")
    assert "parent       vsv_indiana" in out
    assert "derivation   insert after g with zebov_sgp, replication-competent" in out
    assert "  4      zebov_sgp  zaire-ebolavirus-glycoprotein  Mayinga  soluble" in out
    assert "  plasmid  VSVXN2" in out
    assert "  no observations" in out
    assert "locus unknown" not in out


def test_check_accepts_the_indiana_chassis(capsys: pytest.CaptureFixture[str]) -> None:
    status, out, err = run(capsys, "check", "data/designs/vsv_indiana")
    assert (status, err) == (0, "")
    assert out.strip() == "vsv_indiana: ok"


def test_check_reports_a_knowledge_gap_without_failing(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status, out, err = run(capsys, "check", "rvsv_dg_bombali_gp")
    assert (status, err) == (0, "")
    assert out.splitlines()[0].split()[0] == "gap"
    assert "replication-unknown" in out


def test_check_fails_a_structural_violation(
    capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    vaccine = Vaccine(
        "single_cycle",
        "single cycle",
        Replication.SINGLE_CYCLE,
        tuple(Viral(gene) for gene in Gene),
    )
    path = tmp_path / "design.ttl"
    vaccine.write(path)
    status, out, err = run(capsys, "check", str(path))
    assert (status, err) == (1, "")
    assert "replication" in out


def test_precedents_and_antigenome_commands(
    capsys: pytest.CaptureFixture[str], tmp_path: Path
) -> None:
    coding = "ATG" + "GCT" * 40 + "CTA" * 20 + "TAA"
    vaccine = indiana().replace(
        coding,
        identifier="rvsv_dg_outbreak_gp",
        name="rVSV outbreak glycoprotein",
        antigen_class="filovirus-glycoprotein",
        replication=Replication.COMPETENT,
        form=GlycoproteinForm.FULL_LENGTH,
        strain="Kikwit-1995",
    )
    path = tmp_path / "design.ttl"
    vaccine.write(path)
    status, out, err = run(
        capsys,
        "precedents",
        str(path),
        "--antigen-class",
        "filovirus-glycoprotein",
        "--phenotype",
        "protection",
    )
    assert (status, err) == (0, "")
    assert out.splitlines()[0].split() == ["id", "phenotype", "host", "assay", "citation"]
    assert "rvsv_dg_filovirus_gp" not in out
    assert "rvsv_dg_marv_gp" in out
    assert "rvsv_dg_zebov_mayinga_gp" in out
    assert "rvsv_dg_zebov_kikwit_gp" in out
    assert "rvsv_dg_bombali_gp" not in out
    assert "https://doi.org/10.1038/nm1258" in out
    status, out, err = run(capsys, "antigenome", str(path))
    assert (status, err) == (0, "")
    elements, encoding = _sbol_records(out)["rvsv_dg_outbreak_gp"]
    assert encoding == DNA_ENCODING
    assert coding in elements
    status, out, err = run(capsys, "antigenome", "rvsv_dg_zebov_kikwit_gp")
    assert status == 1
    assert out == ""
    assert "no coding sequence" in err


def test_extract_rejects_unknown_genes(
    capsys: pytest.CaptureFixture[str], synthetic_fasta: Path
) -> None:
    status, _out, err = run(capsys, "extract", "-g", "X", str(synthetic_fasta))
    assert status == 2
    assert "invalid gene" in err
