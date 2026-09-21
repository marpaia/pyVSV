"""The SBOL3 genome model: identifiers, roles, coordinates, strands, and serialization."""

import pytest
import sbol3

from vsv import sbol
from vsv.sequence import InputError, reverse_complement


@pytest.fixture
def genome() -> sbol.Genome:
    return sbol.add_genome(sbol.new_document(), "X_1.2", "test genome", "ATGAAATTTGGGTAACCCATGCCC")


@pytest.mark.parametrize(
    ("identifier", "expected"),
    [
        ("NC_001560.1", "NC_001560"),
        ("MF196237.1", "MF196237"),
        ("plain", "plain"),
        ("with-dash.and space", "with_dash_and_space"),
        ("1abc", "_1abc"),
    ],
)
def test_display_id(identifier: str, expected: str) -> None:
    assert sbol.display_id(identifier) == expected
    assert sbol3.Component(expected, sbol3.SBO_DNA) is not None


@pytest.mark.parametrize(
    "uri",
    [
        "https://identifiers.org/SO:0000704",
        "http://identifiers.org/so/SO:0000704",
        "http://purl.obolibrary.org/obo/SO_0000704",
    ],
)
def test_role_key_normalizes_sequence_ontology_uris(uri: str) -> None:
    assert sbol.role_key(uri) == "SO:0000704"


def test_role_key_leaves_other_uris_alone() -> None:
    assert sbol.role_key(sbol.ROLE_LEADER) == sbol.ROLE_LEADER


def test_genome_identity(genome: sbol.Genome) -> None:
    assert genome.id == "X_1.2"
    assert genome.component.identity == "https://vsv.bio/X_1"
    assert genome.description == "test genome"
    assert genome.bases == "ATGAAATTTGGGTAACCCATGCCC"
    assert len(genome) == 24
    assert str(genome.sequence.elements) == genome.bases.lower()
    assert genome.sequence.encoding == sbol3.IUPAC_DNA_ENCODING


def test_add_region_uses_one_based_inclusive_ranges(genome: sbol.Genome) -> None:
    feature = genome.add_region("N", 0, 12, [sbol.SO_CDS])
    location = feature.locations[0]
    assert (location.start, location.end) == (1, 12)
    assert genome.span(feature) == (0, 12)
    assert genome.feature_bases(feature) == "ATGAAATTTGGG"
    assert genome.translate(feature) == "MKFG"


@pytest.mark.parametrize(("start", "end"), [(-1, 4), (4, 4), (5, 4), (0, 25)])
def test_add_region_rejects_spans_outside_the_sequence(
    genome: sbol.Genome, start: int, end: int
) -> None:
    with pytest.raises(ValueError, match="outside a sequence"):
        genome.add_region("bad", start, end, [sbol.SO_REGION])


def test_reverse_strand_features_read_the_reverse_complement(genome: sbol.Genome) -> None:
    location = sbol3.Range(genome.sequence, 13, 24, orientation=sbol3.SBOL_REVERSE_COMPLEMENT)
    feature = sbol3.SequenceFeature([location], roles=[sbol.SO_CDS], name="R")
    genome.component.features.append(feature)
    assert genome.is_reverse(feature)
    assert genome.feature_bases(feature) == reverse_complement(genome.bases[12:24])
    assert not genome.is_reverse(genome.add_region("F", 0, 3, [sbol.SO_CDS]))


def test_features_filter_by_role_across_uri_spellings(genome: sbol.Genome) -> None:
    genome.add_region("gene", 0, 12, [sbol.SO_GENE])
    legacy = sbol3.SequenceFeature(
        [sbol3.Range(genome.sequence, 13, 24)], roles=["http://identifiers.org/so/SO:0000704"]
    )
    genome.component.features.append(legacy)
    assert len(genome.features(sbol.SO_GENE)) == 2
    assert genome.features(sbol.SO_CDS) == []


def test_features_are_ordered_by_position_then_length(genome: sbol.Genome) -> None:
    genome.add_region("late", 12, 24, [sbol.SO_REGION])
    genome.add_region("short", 0, 6, [sbol.SO_REGION])
    genome.add_region("long", 0, 12, [sbol.SO_REGION])
    assert [f.name for f in genome.features()] == ["long", "short", "late"]


def test_named_lookups(genome: sbol.Genome) -> None:
    genome.add_region("N", 0, 12, [sbol.SO_GENE])
    genome.add_region("N", 0, 12, [sbol.SO_CDS])
    assert genome.gene("N") is not None
    assert genome.cds("N") is not None
    assert genome.gene("P") is None
    assert genome.protein("N") == "MKFG"
    assert genome.protein("P") is None


def test_add_protein_links_the_coding_sequence_to_its_product(genome: sbol.Genome) -> None:
    cds = genome.add_region("N", 0, 12, [sbol.SO_CDS])
    protein = genome.add_protein(cds, "N")
    assert protein.identity == "https://vsv.bio/X_1_N_protein"
    assert protein.types == [sbol3.SBO_PROTEIN]
    (sequence,) = (
        o for o in genome.document.objects if o.identity == "https://vsv.bio/X_1_N_protein_seq"
    )
    assert sequence.elements == "MKFG"
    assert sequence.encoding == sbol3.IUPAC_PROTEIN_ENCODING
    (interaction,) = genome.component.interactions
    assert interaction.types == [sbol3.SBO_GENETIC_PRODUCTION]
    by_role = {p.roles[0]: p.participant for p in interaction.participations}
    assert by_role[sbol3.SBO_TEMPLATE] == cds.identity
    (product,) = (f for f in genome.component.features if isinstance(f, sbol3.SubComponent))
    assert by_role[sbol3.SBO_PRODUCT] == product.identity
    assert product.instance_of == protein.identity
    assert not genome.document.validate()


def test_protein_components_are_not_genomes(genome: sbol.Genome) -> None:
    genome.add_protein(genome.add_region("N", 0, 12, [sbol.SO_CDS]), "N")
    assert [g.id for g in sbol.genomes(genome.document)] == ["X_1.2"]


def test_duplicate_identifiers_are_rejected(genome: sbol.Genome) -> None:
    with pytest.raises(InputError, match=r"Duplicate record identifier: X_1\.9"):
        sbol.add_genome(genome.document, "X_1.9", "", "ACGT")


def test_component_without_a_dna_sequence_is_not_a_genome() -> None:
    document = sbol.new_document()
    component = sbol3.Component("empty", sbol3.SBO_DNA)
    document.add(component)
    with pytest.raises(InputError, match="no DNA sequence"):
        sbol.Genome(document, component)


@pytest.mark.parametrize(
    "output_format", [sbol3.TURTLE, sbol3.NTRIPLES, sbol3.RDF_XML, sbol3.JSONLD]
)
def test_annotations_survive_serialization(genome: sbol.Genome, output_format: str) -> None:
    genome.add_region("N", 0, 12, [sbol.SO_GENE, sbol.ROLE_GENE_START])
    genome.add_protein(genome.add_region("N", 0, 12, [sbol.SO_CDS]), "N")
    assert not genome.document.validate()
    reread = sbol3.Document()
    reread.read_string(genome.document.write_string(output_format), output_format)
    (loaded,) = sbol.genomes(reread)
    assert loaded.id == "X_1.2"
    assert loaded.bases == genome.bases
    assert [(loaded.span(f), f.name) for f in loaded.features()] == [
        ((0, 12), "N"),
        ((0, 12), "N"),
    ]
    assert loaded.features(sbol.ROLE_GENE_START)
    assert loaded.protein("N") == "MKFG"
