"""Architectures: canonical projection, order constraints, and SBOL round-trip."""

from dataclasses import replace
from pathlib import Path

import pytest
import sbol3

from tests.synthetic import SyntheticGenome
from vsv import NAMESPACE, Architecture, sbol
from vsv.annotation import annotate
from vsv.architecture import (
    Derivation,
    Module,
    ModuleKind,
    Replication,
    add_order_constraints,
    architecture_of,
    canonical_chassis,
    project,
    to_document,
)
from vsv.formats import read_genomes
from vsv.sbol import display_id
from vsv.sequence import InputError

RAW = Path(__file__).resolve().parents[1] / "data" / "raw"
_DOI = "https://doi.org/10.1371/journal.pcbi.1011912"


def _reload(architecture: Architecture) -> Architecture:
    document = to_document(architecture)
    reread = sbol3.Document()
    reread.read_string(document.write_string(sbol3.TURTLE), sbol3.TURTLE)
    return architecture_of(reread)


def _logical(architecture: Architecture) -> tuple[object, ...]:
    return (
        architecture.id,
        architecture.name,
        architecture.backbone,
        architecture.replication,
        architecture.parent_id,
        architecture.citations,
        architecture.derivation,
        tuple(
            (module.id, module.kind, module.name, module.locus, module.antigen_class)
            for module in architecture.modules
        ),
    )


def test_canonical_chassis_is_leader_five_genes_and_trailer() -> None:
    architecture = canonical_chassis()
    assert [(module.id, module.locus) for module in architecture.modules] == [
        ("leader", None),
        ("Nucleoprotein", 0),
        ("Phosphoprotein", 1),
        ("MatrixProtein", 2),
        ("Glycoprotein", 3),
        ("LargePolymerase", 4),
        ("trailer", None),
    ]
    assert architecture.replication is Replication.COMPETENT
    assert architecture.parent_id is None


def test_ids_are_sbol_display_ids() -> None:
    with pytest.raises(InputError, match="architecture id"):
        replace(canonical_chassis(), id="vsv-indiana")
    with pytest.raises(InputError, match="module id"):
        Module("antigen-1", ModuleKind.ANTIGEN, "antigen", antigen_class="filovirus-glycoprotein")


def test_an_antigen_needs_a_class_slug() -> None:
    with pytest.raises(InputError, match="antigen class"):
        Module("antigen_1", ModuleKind.ANTIGEN, "antigen")
    with pytest.raises(InputError, match="not an antigen"):
        Module(
            "Nucleoprotein",
            ModuleKind.VIRAL_GENE,
            "Nucleoprotein",
            locus=0,
            antigen_class="filovirus-glycoprotein",
        )


def test_module_ids_are_unique() -> None:
    architecture = canonical_chassis()
    duplicated = architecture.modules[:1] + architecture.modules[:1]
    with pytest.raises(InputError, match="unique"):
        replace(architecture, modules=duplicated)


def test_project_reads_an_annotated_synthetic_genome(synthetic: SyntheticGenome) -> None:
    genome = sbol.add_genome(
        sbol.new_document(), synthetic.identifier, synthetic.description, synthetic.sequence
    )
    annotate(genome)
    architecture = project(genome, backbone="VSV")
    assert architecture.id == display_id(synthetic.identifier)
    assert architecture.name == synthetic.identifier
    assert [(module.id, module.locus) for module in architecture.modules] == [
        ("leader", None),
        ("Nucleoprotein", 0),
        ("Phosphoprotein", 1),
        ("MatrixProtein", 2),
        ("Glycoprotein", 3),
        ("LargePolymerase", 4),
        ("trailer", None),
    ]
    nucleoprotein = genome.gene("Nucleoprotein")
    assert nucleoprotein is not None
    assert architecture.modules[1].feature_identity == nucleoprotein.identity


def test_project_requires_the_five_coding_sequences(synthetic: SyntheticGenome) -> None:
    genome = sbol.add_genome(sbol.new_document(), synthetic.identifier, "", synthetic.sequence)
    with pytest.raises(InputError, match="coding sequences"):
        project(genome)


def test_order_constraints_follow_the_canonical_features(synthetic: SyntheticGenome) -> None:
    genome = sbol.add_genome(sbol.new_document(), synthetic.identifier, "", synthetic.sequence)
    annotate(genome)
    add_order_constraints(genome)
    add_order_constraints(genome)
    names = []
    for constraint in genome.component.constraints:
        if constraint.restriction != sbol3.SBOL_STRICTLY_PRECEDES:
            continue
        earlier = genome.document.find(str(constraint.subject))
        later = genome.document.find(str(constraint.object))
        names.append((earlier.name, later.name))
    assert names == [
        ("leader", "Nucleoprotein"),
        ("Nucleoprotein", "Phosphoprotein"),
        ("Phosphoprotein", "MatrixProtein"),
        ("MatrixProtein", "Glycoprotein"),
        ("Glycoprotein", "LargePolymerase"),
        ("LargePolymerase", "trailer"),
    ]
    assert not genome.document.validate()


def test_architecture_round_trip_keeps_order_lineage_and_citations() -> None:
    antigen = Module(
        "antigen_1",
        ModuleKind.ANTIGEN,
        "filovirus glycoprotein",
        locus=3,
        antigen_class="filovirus-glycoprotein",
    )
    modules = tuple(
        module if module.id != "Glycoprotein" else antigen for module in canonical_chassis().modules
    )
    architecture = Architecture(
        "rvsv_dg",
        "rVSV filovirus GP",
        "VSV-Indiana",
        Replication.SINGLE_CYCLE,
        modules,
        parent_id="vsv_indiana",
        citations=(_DOI, "https://pubmed.ncbi.nlm.nih.gov/36508878"),
        derivation=Derivation(
            "replace-module", "Glycoprotein", "antigen_1", Replication.SINGLE_CYCLE
        ),
    )
    document = to_document(architecture)
    assert not document.validate()
    loaded = _reload(architecture)
    assert _logical(loaded) == _logical(architecture)
    assert loaded.modules[0].feature_identity == f"{NAMESPACE}/rvsv_dg/leader"
    glycoprotein = next(module for module in loaded.modules if module.id == "antigen_1")
    assert glycoprotein.feature_identity == f"{NAMESPACE}/rvsv_dg/antigen_1"


def test_a_genome_is_not_an_architecture(synthetic: SyntheticGenome) -> None:
    genome = sbol.add_genome(sbol.new_document(), synthetic.identifier, "", synthetic.sequence)
    with pytest.raises(InputError, match="found 0"):
        architecture_of(genome.document)


def test_unknown_locus_and_replication_are_not_serialized() -> None:
    architecture = canonical_chassis()
    unknown = replace(
        next(module for module in architecture.modules if module.id == "LargePolymerase"),
        locus=None,
    )
    modules = tuple(
        unknown if module.id == "LargePolymerase" else module for module in architecture.modules
    )
    with pytest.raises(InputError, match="needs a locus"):
        to_document(replace(architecture, modules=modules))
    with pytest.raises(InputError, match="replication is unknown"):
        to_document(replace(architecture, replication=Replication.UNKNOWN))


def test_serialization_requires_locus_order() -> None:
    architecture = canonical_chassis()
    gap = [
        replace(module, locus=5) if module.id == "LargePolymerase" else module
        for module in architecture.modules
    ]
    with pytest.raises(InputError, match="locus order"):
        to_document(replace(architecture, modules=tuple(gap)))


def test_reference_indiana_projects_to_the_five_genes() -> None:
    path = RAW / "NC_001560.gb"
    if not path.exists():
        pytest.skip("NC_001560 GenBank file is not in data/raw")
    _, (genome,) = read_genomes([path])
    architecture = project(genome, backbone="VSV-Indiana")
    transcribed = [
        module for module in architecture.modules if module.kind is ModuleKind.VIRAL_GENE
    ]
    assert [(module.id, module.locus) for module in transcribed] == [
        ("Nucleoprotein", 0),
        ("Phosphoprotein", 1),
        ("MatrixProtein", 2),
        ("Glycoprotein", 3),
        ("LargePolymerase", 4),
    ]
    assert architecture.replication is Replication.COMPETENT
