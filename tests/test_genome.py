"""FASTA, GenBank, and SBOL genomes, and the annotations read from them."""

from pathlib import Path

import pytest

from tests.synthetic import SyntheticGenome, genbank_text
from vsv import GeneWarning, Genome, InputError, indiana, load, load_all, read_record
from vsv.genes import Gene, Product
from vsv.genome.annotate import annotate
from vsv.genome.formats import read_records
from vsv.genome.sbol import write_genome

RAW = Path(__file__).resolve().parents[1] / "data" / "raw"
ACCESSIONS = ["NC_001560", "MF196237", "MF196238"]


def test_fasta_annotation_recovers_the_synthetic_genome(synthetic: SyntheticGenome) -> None:
    genome = annotate(
        synthetic.sequence, identifier=synthetic.identifier, description=synthetic.description
    )
    assert genome.leader is not None
    assert (genome.leader.start, genome.leader.end) == synthetic.leader
    assert genome.trailer is not None
    assert (genome.trailer.start, genome.trailer.end) == synthetic.trailer
    for gene in Gene:
        found = genome.gene(gene)
        assert (found.unit.start, found.unit.end) == synthetic.genes[gene]
        assert found.cds is not None
        assert (found.cds.start, found.cds.end) == synthetic.cds[gene]
        assert genome.protein(gene) == synthetic.proteins[gene]


def test_five_genes_are_named_in_file_order(tmp_path: Path, synthetic: SyntheticGenome) -> None:
    path = tmp_path / "genes.gb"
    path.write_text(genbank_text(synthetic, genes=True))
    genome = read_record(path)
    assert [gene.gene for gene in genome.genes] == list(Gene)
    assert genome.id == "SYN0001.1"


def test_five_coding_sequences_name_the_genes(tmp_path: Path, synthetic: SyntheticGenome) -> None:
    path = tmp_path / "cds.gb"
    path.write_text(genbank_text(synthetic))
    genome = read_record(path)
    for gene in Gene:
        found = genome.gene(gene)
        start, end = synthetic.cds[gene]
        assert found.cds is not None
        assert (found.unit.start, found.unit.end) == (start, end)
        assert (found.cds.start, found.cds.end) == (start, end)


def test_a_shorter_cds_inside_phosphoprotein_is_not_p(
    tmp_path: Path, synthetic: SyntheticGenome
) -> None:
    start, _end = synthetic.cds[Gene.PHOSPHOPROTEIN]
    partial = (
        f"     CDS             <{start + 5}..>{start + 9}\n"
        '                     /product="protein C\'"\n'
    )
    path = tmp_path / "partial.gb"
    path.write_text(genbank_text(synthetic, genes=True).replace("ORIGIN", partial + "ORIGIN"))
    genome = read_record(path)
    coding = genome.gene(Gene.PHOSPHOPROTEIN)
    assert coding.cds is not None
    assert (coding.cds.start, coding.cds.end) == synthetic.cds[Gene.PHOSPHOPROTEIN]
    assert coding.product == "phosphoprotein"


def test_load_warns_and_read_record_does_not(tmp_path: Path, synthetic: SyntheticGenome) -> None:
    path = tmp_path / "short.fasta"
    path.write_text(">SHORT\nACGT\n")
    with pytest.warns(GeneWarning, match="genes not identified: Nucleoprotein"):
        assert not load(path).genes
    assert not read_record(path).genes


def test_load_rejects_more_than_one_record(tmp_path: Path, synthetic: SyntheticGenome) -> None:
    path = tmp_path / "two.fasta"
    path.write_text(synthetic.fasta + synthetic.fasta.replace("SYN0001.1", "SYN0002.1"))
    with pytest.raises(InputError, match="found 2; use load_all"):
        load(path)
    genomes = load_all([path])
    assert [genome.id for genome in genomes] == ["SYN0001.1", "SYN0002.1"]


def test_duplicate_and_unknown_files_name_the_problem(
    tmp_path: Path, synthetic: SyntheticGenome
) -> None:
    fasta = tmp_path / "a.fasta"
    fasta.write_text(synthetic.fasta)
    with pytest.raises(InputError, match="Duplicate record identifier"):
        read_records([fasta, fasta])
    unknown = tmp_path / "genome.txt"
    unknown.write_text(">A\nACGT\n")
    with pytest.raises(InputError, match=rf"{unknown}: Unrecognized file type"):
        read_records([unknown])
    bad = tmp_path / "bad.fasta"
    bad.write_text("ACGT\n")
    with pytest.raises(InputError, match=rf"{bad}:") as error:
        read_records([bad])
    assert "before the first '>' header" in str(error.value)


@pytest.mark.parametrize("accession", ACCESSIONS)
def test_annotated_spans_match_the_curated_record(accession: str) -> None:
    fasta = RAW / f"{accession}.fasta"
    genbank = RAW / f"{accession}.gb"
    if not fasta.exists() or not genbank.exists():
        pytest.skip(f"{accession} FASTA and GenBank files are not in data/raw")
    heuristic = read_record(fasta)
    curated = read_record(genbank)
    assert heuristic.bases == curated.bases
    for gene in Gene:
        found = heuristic.gene(gene)
        expected = curated.gene(gene)
        assert found.cds is not None and expected.cds is not None
        assert (found.cds.start, found.cds.end) == (expected.cds.start, expected.cds.end)
        assert found.unit.start <= found.cds.start < found.cds.end <= found.unit.end
        assert expected.unit.start <= expected.cds.start < expected.cds.end <= expected.unit.end
        assert heuristic.protein(gene) == curated.protein(gene)


@pytest.mark.parametrize("name", ["NC_001560.gb", "NC_001560.ttl"])
def test_the_partial_phosphoprotein_orf_is_not_the_p_cds(name: str) -> None:
    path = RAW / name
    if not path.exists():
        pytest.skip(f"{name} is not in data/raw")
    genome = read_record(path)
    coding = genome.gene(Gene.PHOSPHOPROTEIN)
    assert coding.cds is not None
    assert (coding.cds.start, coding.cds.end) == (1395, 2193)
    assert coding.product == "NS protein"
    assert genome.id == "NC_001560.1"


def test_the_indiana_document_matches_the_genbank_record() -> None:
    genbank = RAW / "NC_001560.gb"
    turtle = RAW / "NC_001560.ttl"
    if not genbank.exists() or not turtle.exists():
        pytest.skip("NC_001560 files are not in data/raw")
    curated = read_record(genbank)
    document = read_record(turtle)
    assert document.bases == curated.bases
    assert document.id == curated.id == "NC_001560.1"
    for gene in Gene:
        assert document.gene(gene) == curated.gene(gene)


def test_a_genome_document_round_trips(tmp_path: Path, synthetic: SyntheticGenome) -> None:
    genome = annotate(
        synthetic.sequence, identifier=synthetic.identifier, description=synthetic.description
    )
    path = tmp_path / "genome.ttl"
    path.write_text(write_genome(genome))
    assert read_record(path) == genome


def test_a_sequence_only_document_is_annotated(tmp_path: Path, synthetic: SyntheticGenome) -> None:
    path = tmp_path / "sequence.ttl"
    path.write_text(
        write_genome(Genome(synthetic.identifier, synthetic.description, synthetic.sequence))
    )
    assert read_record(path) == annotate(
        synthetic.sequence, identifier=synthetic.identifier, description=synthetic.description
    )


def test_a_design_file_has_no_genome_sequence(tmp_path: Path) -> None:
    path = tmp_path / "design.ttl"
    indiana().write(path)
    with pytest.raises(InputError, match="no genome sequence"):
        read_record(path)


def test_unreadable_sbol_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "bad.ttl"
    path.write_text("this is not turtle\n")
    with pytest.raises(InputError, match="not a readable SBOL file"):
        read_record(path)


def test_sequences_are_read_by_product(synthetic_fasta: Path, synthetic: SyntheticGenome) -> None:
    genome = read_record(synthetic_fasta)
    gene = genome.gene(Gene.MATRIX_PROTEIN)
    start, end = synthetic.cds[Gene.MATRIX_PROTEIN]
    assert gene.sequence(genome.bases, Product.CDS) == synthetic.sequence[start:end]
    assert genome.protein(Gene.MATRIX_PROTEIN) == synthetic.proteins[Gene.MATRIX_PROTEIN]
