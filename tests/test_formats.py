"""Reading genomes from FASTA, GenBank, and SBOL3 files."""

from pathlib import Path

import pytest
import sbol3

from tests.synthetic import SyntheticGenome, genbank_text
from vsv import sbol
from vsv.annotation import annotate
from vsv.formats import OUTPUT_FORMATS, SBOL_SUFFIXES, read_genomes, write_document
from vsv.sequence import InputError


def test_fasta_becomes_an_sbol_genome(synthetic_fasta: Path, synthetic: SyntheticGenome) -> None:
    document, genomes = read_genomes([synthetic_fasta])
    (genome,) = genomes
    assert genome.id == "SYN0001.1"
    assert genome.component.identity == "https://vsv.bio/SYN0001"
    assert genome.description == synthetic.description
    assert genome.bases == synthetic.sequence
    assert genome.features() == []
    assert not document.validate()


def test_genbank_becomes_an_sbol_genome_with_named_genes(
    synthetic_genbank: Path, synthetic: SyntheticGenome
) -> None:
    document, (genome,) = read_genomes([synthetic_genbank])
    assert genome.component.identity == "https://vsv.bio/SYN0001"
    assert genome.bases == synthetic.sequence
    assert genome.description.startswith("synthetic vesicular stomatitis virus")
    assert {str(f.name): genome.span(f) for f in genome.features(sbol.SO_CDS)} == synthetic.cds
    assert genome.protein("L") == synthetic.proteins["L"]
    assert not document.validate()


def test_genbank_genes_are_named_by_position(tmp_path: Path, synthetic: SyntheticGenome) -> None:
    path = tmp_path / "genes.gb"
    path.write_text(genbank_text(synthetic, genes=True))
    _, (genome,) = read_genomes([path])
    assert {str(f.name): genome.span(f) for f in genome.features(sbol.SO_GENE)} == synthetic.genes


def test_fasta_and_genbank_agree_on_the_genome(
    synthetic_fasta: Path, synthetic_genbank: Path, tmp_path: Path
) -> None:
    _, (from_fasta,) = read_genomes([synthetic_fasta])
    _, (from_genbank,) = read_genomes([synthetic_genbank])
    assert from_fasta.bases == from_genbank.bases
    assert from_fasta.component.identity == from_genbank.component.identity


@pytest.mark.parametrize("output_format", OUTPUT_FORMATS)
def test_sbol_files_round_trip(synthetic_fasta: Path, tmp_path: Path, output_format: str) -> None:
    document, (genome,) = read_genomes([synthetic_fasta])
    annotate(genome)
    suffix = next(s for s, f in SBOL_SUFFIXES.items() if f == OUTPUT_FORMATS[output_format])
    path = tmp_path / f"genome{suffix}"
    path.write_text(write_document(document, output_format))
    reread, (loaded,) = read_genomes([path])
    assert loaded.id == genome.id
    assert loaded.bases == genome.bases
    assert [(loaded.span(f), f.name) for f in loaded.features()] == [
        (genome.span(f), f.name) for f in genome.features()
    ]
    assert not reread.validate()


def test_inputs_are_merged_in_order(tmp_path: Path, synthetic_fasta: Path) -> None:
    other = tmp_path / "other.fasta"
    other.write_text(">AAA\nACGT\n>ZZZ\nGGCC\n")
    document, genomes = read_genomes([other, synthetic_fasta])
    assert [g.id for g in genomes] == ["AAA", "ZZZ", "SYN0001.1"]
    assert {g.document for g in genomes} == {document}
    assert len(sbol.genomes(document)) == 3


def test_mixed_formats_share_one_document(synthetic_genbank: Path, tmp_path: Path) -> None:
    other = tmp_path / "other.fasta"
    other.write_text(">OTHER\nACGT\n")
    document, genomes = read_genomes([synthetic_genbank, other])
    assert [g.id for g in genomes] == ["SYN0001", "OTHER"]
    assert not document.validate()


@pytest.mark.parametrize("second", ["fasta", "gb"])
def test_duplicate_genomes_are_rejected(
    synthetic_fasta: Path, synthetic_genbank: Path, second: str
) -> None:
    other = synthetic_fasta if second == "fasta" else synthetic_genbank
    with pytest.raises(InputError, match="Duplicate record identifier"):
        read_genomes([synthetic_fasta, other])


def test_unknown_file_types_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "genome.txt"
    path.write_text(">A\nACGT\n")
    with pytest.raises(InputError, match=r"genome\.txt: Unrecognized file type \.txt"):
        read_genomes([path])


def test_missing_files_are_rejected(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="No such file or directory"):
        read_genomes([tmp_path / "missing.fasta"])


def test_files_without_records_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "empty.fasta"
    path.write_text("")
    with pytest.raises(InputError, match="No sequence records"):
        read_genomes([path])


@pytest.mark.parametrize(
    ("name", "text", "message"),
    [
        ("bad.fasta", "ACGT\n", "before the first '>' header"),
        ("bad.gb", "this is not GenBank\n", "not a readable GenBank file"),
        ("empty.gb", "", "not a readable GenBank file"),
        ("bad.ttl", "@prefix : <x> .\nnot turtle", "not a readable SBOL3 file"),
        ("bad.xml", "<?xml version='1.0'?><unclosed", "not a readable SBOL3 file"),
        ("bad.jsonld", "{not json", "not a readable SBOL3 file"),
    ],
)
def test_malformed_files_name_the_file(tmp_path: Path, name: str, text: str, message: str) -> None:
    path = tmp_path / name
    path.write_text(text)
    with pytest.raises(InputError, match=message) as error:
        read_genomes([path])
    assert str(error.value).startswith(f"{path}:")


def test_sbol_without_dna_sequences_has_no_records(tmp_path: Path) -> None:
    document = sbol.new_document()
    document.add(sbol3.Component("nothing", sbol3.SBO_DNA))
    path = tmp_path / "nothing.ttl"
    path.write_text(document.write_string(sbol3.TURTLE))
    with pytest.raises(InputError, match="No sequence records"):
        read_genomes([path])
