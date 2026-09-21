"""Loading genomes with their genes identified, and typed access to each gene."""

import warnings
from pathlib import Path

import pytest

import vsv
from tests.synthetic import SyntheticGenome
from vsv import sbol
from vsv.api import gene_problems
from vsv.genes import Gene, Product


def test_gene_is_a_string_enum_in_genome_order() -> None:
    assert list(Gene) == ["N", "P", "M", "G", "L"]
    assert Gene("P") is Gene.P
    assert Gene.N == "N"


def test_load_annotates_a_fasta_file(synthetic_fasta: Path, synthetic: SyntheticGenome) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        genome = vsv.load(synthetic_fasta)
    assert list(genome.genes) == list(Gene)
    for gene, annotation in genome.genes.items():
        assert annotation.gene_span == synthetic.genes[gene]
        assert annotation.cds_span == synthetic.cds[gene]
        assert annotation.protein == synthetic.proteins[gene]


def test_load_accepts_a_string_path(synthetic_fasta: Path) -> None:
    assert vsv.load(str(synthetic_fasta)).genes[Gene.N].cds_span is not None


def test_load_names_a_curated_genbank_record(
    synthetic_genbank: Path, synthetic: SyntheticGenome
) -> None:
    genome = vsv.load(synthetic_genbank)
    assert list(genome.genes) == list(Gene)
    assert genome.genes[Gene.G].cds_span == synthetic.cds["G"]
    assert genome.genes[Gene.G].gene_span is None


def test_sequences_are_read_by_product(synthetic_fasta: Path, synthetic: SyntheticGenome) -> None:
    annotation = vsv.load(synthetic_fasta).genes[Gene.M]
    start, end = synthetic.cds["M"]
    assert annotation.sequence(Product.CDS) == synthetic.sequence[start:end]
    assert annotation.sequence(Product.PROTEIN) == synthetic.proteins["M"]
    gene_start, gene_end = synthetic.genes["M"]
    assert annotation.sequence(Product.GENE) == synthetic.sequence[gene_start:gene_end]


def test_incomplete_genomes_warn(tmp_path: Path, synthetic: SyntheticGenome) -> None:
    truncated = tmp_path / "short.fasta"
    end = synthetic.genes["G"][1] + 7
    truncated.write_text(f">SHORT\n{synthetic.sequence[:end]}\n")
    with pytest.warns(vsv.GeneWarning, match=r"SHORT: genes not identified: N, P, M, G, L"):
        genome = vsv.load(truncated)
    assert not genome.genes


def test_gene_problems_names_genes_without_a_coding_sequence(synthetic: SyntheticGenome) -> None:
    genome = sbol.add_genome(sbol.new_document(), "ONLY", "", synthetic.sequence)
    for gene, (start, end) in synthetic.genes.items():
        genome.add_region(gene, start, end, [sbol.SO_GENE])
    assert gene_problems(genome) == ["genes without a coding sequence: N, P, M, G, L"]


def test_load_rejects_multi_record_files(tmp_path: Path, synthetic: SyntheticGenome) -> None:
    path = tmp_path / "two.fasta"
    path.write_text(synthetic.fasta + synthetic.fasta.replace("SYN0001.1", "SYN0002.1"))
    with pytest.raises(vsv.InputError, match="found 2; use load_all"):
        vsv.load(path)
    _, genomes = vsv.load_all([path])
    assert [g.id for g in genomes] == ["SYN0001.1", "SYN0002.1"]
