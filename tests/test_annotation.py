"""Annotation of VSV genomes as SBOL3 features."""

import pytest
import sbol3

import vsv.genome.sbol as sbol
from tests.synthetic import SyntheticGenome, build_genome
from vsv.genes import Gene
from vsv.genome.annotation import (
    GENE_END_SIGNAL,
    TRANSCRIPT_END,
    annotate,
    find_gene_ends,
    find_gene_start,
    find_orf,
    is_annotated,
    kind,
    name_curated_genes,
    to_gff3,
)


def make(synthetic: SyntheticGenome) -> sbol.Genome:
    return sbol.add_genome(
        sbol.new_document(), synthetic.identifier, synthetic.description, synthetic.sequence
    )


def spans(genome: sbol.Genome, role: str) -> dict[str, tuple[int, int]]:
    return {str(f.name): genome.span(f) for f in genome.features(role)}


def test_find_gene_ends_reports_every_signal() -> None:
    signal = "TATG" + "A" * 7
    assert find_gene_ends(f"CC{signal}GG{signal}") == [2, 15]
    assert find_gene_ends("TATGAAAAAA") == []


def test_find_gene_start_tolerates_one_mismatch_and_prefers_exact() -> None:
    assert find_gene_start("GGAACAGTTATCGG", 0, 14) == 2
    assert find_gene_start("GGAACAGTTTTCGG", 0, 14) == 2
    assert find_gene_start("GGAACTTTTTTCGG", 0, 14) is None
    assert find_gene_start("AACTGTTATCAACAGTTATC", 0, 20) == 10


def test_find_orf_needs_a_start_and_a_closed_frame() -> None:
    assert find_orf("CCATGGCTTAACC", 0, 13) == (2, 11)
    assert find_orf("CCGCTTAACC", 0, 10) is None
    assert find_orf("CCATGGCTGCT", 0, 11) is None


def test_annotate_finds_every_feature(synthetic: SyntheticGenome) -> None:
    genome = make(synthetic)
    assert not is_annotated(genome)
    annotate(genome)
    assert is_annotated(genome)
    assert spans(genome, sbol.SO_GENE) == synthetic.genes
    assert {name: genome.span(genome.cds(name)) for name in Gene} == synthetic.cds
    assert list(spans(genome, sbol.SO_GENE)) == list(Gene)
    coding = genome.cds(Gene.NUCLEOPROTEIN)
    assert coding is not None
    assert coding.name == "Nucleoprotein CDS"
    assert coding.display_id == "Nucleoprotein_cds"
    gene = genome.gene(Gene.NUCLEOPROTEIN)
    assert gene is not None
    assert any(
        constraint.restriction == sbol3.SBOL_CONTAINS
        and constraint.subject == gene.identity
        and constraint.object == coding.identity
        for constraint in genome.component.constraints
    )
    assert spans(genome, sbol.ROLE_LEADER) == {"leader": synthetic.leader}
    assert spans(genome, sbol.ROLE_TRAILER) == {"trailer": synthetic.trailer}


def test_annotate_records_signals_and_proteins(synthetic: SyntheticGenome) -> None:
    genome = make(synthetic)
    annotate(genome)
    starts = spans(genome, sbol.ROLE_GENE_START)
    ends = spans(genome, sbol.ROLE_GENE_END)
    for name, (start, end) in synthetic.genes.items():
        assert starts[f"{name} gene start"] == (start, start + 10)
        assert ends[f"{name} gene end"] == (end - 4, end + 7)
        assert genome.protein(name) == synthetic.proteins[name]
    assert len(genome.component.interactions) == len(Gene)
    assert not genome.document.validate()


def test_signals_carry_standard_and_vsv_roles(synthetic: SyntheticGenome) -> None:
    genome = make(synthetic)
    annotate(genome)
    (leader,) = genome.features(sbol.ROLE_LEADER)
    assert {sbol.role_key(r) for r in leader.roles} == {"SO:0000001", sbol.ROLE_LEADER}
    assert genome.features(sbol.SO_TERMINATOR) == genome.features(sbol.ROLE_GENE_END)
    assert genome.features(sbol.SO_REGULATORY_REGION) == genome.features(sbol.ROLE_GENE_START)
    assert [kind(f) for f in genome.features()][:4] == ["leader", "gene", "gene_start", "CDS"]


def test_unexpected_gene_counts_get_positional_names() -> None:
    synthetic = build_genome()
    genome = make(synthetic)
    truncated = sbol.add_genome(
        genome.document,
        "SHORT",
        "",
        synthetic.sequence[
            : synthetic.genes[Gene.GLYCOPROTEIN][1] + len(GENE_END_SIGNAL) - len(TRANSCRIPT_END)
        ],
    )
    annotate(truncated)
    assert [f.name for f in truncated.features(sbol.SO_GENE)] == [f"gene{i}" for i in range(1, 5)]


def test_sequence_without_signals_gets_no_features() -> None:
    genome = sbol.add_genome(sbol.new_document(), "PLAIN", "", "ACGT" * 30)
    annotate(genome)
    assert genome.features() == []
    assert not is_annotated(genome)


def curated(genome: sbol.Genome, synthetic: SyntheticGenome, *, genes: bool) -> None:
    for i, name in enumerate(Gene):
        if genes:
            genome.add_region(f"gene {i}", *synthetic.genes[name], [sbol.SO_GENE])
        genome.add_region(f"product {i}", *synthetic.cds[name], [sbol.SO_CDS])


def test_curated_cds_are_named_by_order(synthetic: SyntheticGenome) -> None:
    genome = make(synthetic)
    curated(genome, synthetic, genes=False)
    name_curated_genes(genome)
    assert {name: genome.span(genome.cds(name)) for name in Gene} == synthetic.cds
    assert genome.protein(Gene.GLYCOPROTEIN) == synthetic.proteins[Gene.GLYCOPROTEIN]
    assert [f.name for f in genome.features(sbol.SO_CDS)] == [f"product {i}" for i in range(5)]
    nucleocapsid = genome.cds(Gene.NUCLEOPROTEIN)
    assert nucleocapsid is not None and nucleocapsid.display_id == "Nucleoprotein_cds"


def test_curated_genes_take_their_longest_cds(synthetic: SyntheticGenome) -> None:
    genome = make(synthetic)
    curated(genome, synthetic, genes=True)
    start, _ = synthetic.cds[Gene.PHOSPHOPROTEIN]
    genome.add_region("short orf", start + 3, start + 12, [sbol.SO_CDS])
    name_curated_genes(genome)
    assert spans(genome, sbol.SO_GENE) == synthetic.genes
    phosphoprotein = genome.cds(Gene.PHOSPHOPROTEIN)
    assert phosphoprotein is not None
    assert genome.span(phosphoprotein) == synthetic.cds[Gene.PHOSPHOPROTEIN]
    assert [f.name for f in genome.features(sbol.SO_CDS)].count("short orf") == 1


def test_naming_curated_genes_is_idempotent(synthetic: SyntheticGenome) -> None:
    genome = make(synthetic)
    curated(genome, synthetic, genes=True)
    name_curated_genes(genome)
    before = [(f.name, f.description) for f in genome.features()]
    name_curated_genes(genome)
    assert [(f.name, f.description) for f in genome.features()] == before


def test_naming_leaves_other_gene_counts_alone(synthetic: SyntheticGenome) -> None:
    genome = make(synthetic)
    genome.add_region("only", *synthetic.cds[Gene.NUCLEOPROTEIN], [sbol.SO_CDS])
    name_curated_genes(genome)
    assert [f.name for f in genome.features()] == ["only"]


def test_gff3_uses_one_based_inclusive_coordinates(synthetic: SyntheticGenome) -> None:
    genome = make(synthetic)
    annotate(genome)
    lines = to_gff3([genome]).splitlines()
    assert lines[:2] == ["##gff-version 3", f"##sequence-region SYN0001.1 1 {len(genome)}"]
    start, end = synthetic.cds[Gene.NUCLEOPROTEIN]
    assert (
        "SYN0001.1\tvsv\tCDS\t"
        f"{start + 1}\t{end}\t.\t+\t0\t"
        "ID=SYN0001.1:CDS:Nucleoprotein CDS;Name=Nucleoprotein CDS"
    ) in lines
    assert (
        f"SYN0001.1\tvsv\tleader\t1\t{synthetic.leader[1]}\t.\t+\t.\tID=SYN0001.1:leader:leader;Name=leader"
        in lines
    )


def test_gff3_reports_the_reverse_strand(synthetic: SyntheticGenome) -> None:
    genome = make(synthetic)
    location = sbol3.Range(genome.sequence, 1, 9, orientation=sbol3.SBOL_REVERSE_COMPLEMENT)
    genome.component.features.append(
        sbol3.SequenceFeature([location], roles=[sbol.SO_GENE], name="R")
    )
    (line,) = (ln for ln in to_gff3([genome]).splitlines() if "\tgene\t" in ln)
    assert line.split("\t")[6] == "-"


@pytest.mark.parametrize(
    ("role", "expected"),
    [
        ("https://identifiers.org/SO:0000234", "mRNA"),
        ("http://identifiers.org/so/SO:0000419", "mature_protein_region"),
        ("https://identifiers.org/SO:0009999", "region"),
    ],
)
def test_kind_maps_roles_to_labels(role: str, expected: str) -> None:
    genome = sbol.add_genome(sbol.new_document(), "K", "", "ACGT" * 4)
    assert kind(genome.add_region("x", 0, 4, [role])) == expected
