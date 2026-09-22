"""The curated design library, its queries, and the linked intracellular model."""

from pathlib import Path

import pytest
import sbol3

from vsv.design import (
    Antigen,
    Architecture,
    Design,
    EvidenceLevel,
    Host,
    Implementation,
    Leader,
    LocusStatus,
    Observation,
    Outcome,
    ParentStatus,
    Phenotype,
    Replication,
    Stage,
    Trailer,
    Unmodeled,
    canonical_chassis,
    check,
    design_to_document,
    load_design,
    load_library,
    parameters_for,
    select,
    to_document,
)
from vsv.design.evidence import Assay
from vsv.design.model import attach_model
from vsv.errors import InputError
from vsv.genes import Gene


def test_seed_designs_match_their_recorded_transformations() -> None:
    records = {record.architecture.id: record for record in load_library()}
    assert set(records) == {
        "vsv_indiana",
        "rvsv_dg_filovirus_gp",
        "rvsv_dg_marv_gp",
        "rvsv_niv_g",
        "rvsv_niv_f",
        "rvsv_lasv_gpc",
        "rvsv_dg_sars_cov_2_s",
        "rvsv_hiv_env",
        "rvsv_dg_hiv_env",
        "rvsv_influenza_ha",
        "rvsv_gfp",
        "rvsv_luc",
        "rvsv_zika_e",
        "rvsv_cchf_gpc",
        "rvsv_rsv_f",
        "rvsv_hantavirus_gp",
        "rvsv_chikv_e",
        "rvsv_rabies_g",
    }
    indiana = records["vsv_indiana"]
    assert indiana.realization == "NC_001560"
    assert indiana.model_source == "https://doi.org/10.1371/journal.pcbi.1011373"
    assert indiana.architecture.parent is ParentStatus.ROOT
    assert indiana.implementations == ()
    assert check(indiana.architecture, indiana.observations) == []
    assert indiana.card.evidence_level is EvidenceLevel.NONE

    filo = records["rvsv_dg_filovirus_gp"].architecture
    applied = indiana.architecture.replace(
        Gene.GLYCOPROTEIN,
        Antigen(Gene.GLYCOPROTEIN, "filovirus glycoprotein", "filovirus-glycoprotein"),
        filo.id,
        name=filo.name,
    )
    assert [
        (module.id, module.kind, module.name, module.locus, module.antigen_class)
        for module in filo.modules
    ] == [
        (module.id, module.kind, module.name, module.locus, module.antigen_class)
        for module in applied.modules
    ]
    assert filo.edit == applied.edit
    assert filo.replication == applied.replication
    assert records["rvsv_dg_filovirus_gp"].card.evidence_level is EvidenceLevel.NONHUMAN_PRIMATE
    assert "Ervebo" in records["rvsv_dg_filovirus_gp"].aliases

    nipah = records["rvsv_niv_g"].architecture
    added = next(module for module in nipah.modules if module.id == "niv_g")
    assert nipah.parent_id == filo.id
    assert nipah.edit is None
    assert added.locus is None
    assert added.locus_status is LocusStatus.UNKNOWN
    assert next(module.locus for module in nipah.modules if module.id == "LargePolymerase") == 4
    fusion = records["rvsv_niv_f"].architecture
    assert fusion.edit is None
    assert next(module for module in fusion.modules if module.id == "niv_f").locus is None
    assert records["rvsv_niv_f"].card.evidence_level is EvidenceLevel.NONHUMAN_PRIMATE

    lassa = records["rvsv_lasv_gpc"]
    assert [(item.id, item.stage) for item in lassa.implementations] == [
        ("garbutt-2004", Stage.PRECLINICAL),
        ("malkin-2025", Stage.CLINICAL),
    ]
    assert {item.implementation_id for item in lassa.observations} == {"malkin-2025"}
    assert lassa.card.evidence_level is EvidenceLevel.HUMAN

    marburg = records["rvsv_dg_marv_gp"]
    assert marburg.card.antigen_classes == ("marburgvirus-glycoprotein",)
    assert marburg.card.evidence_level is EvidenceLevel.NONHUMAN_PRIMATE
    assert marburg.observations[0].host is Host.OTHER
    assert "Zaire" not in marburg.observations[0].note
    ebola = records["rvsv_dg_filovirus_gp"].observations[0]
    assert ebola.implementation_id == "garbutt-2004"
    assert "Zaire" in ebola.note


def test_queries_separate_antigen_class_and_evidence() -> None:
    records = load_library()
    glycoprotein = select(records, antigen_class="henipavirus-glycoprotein")
    assert [record.architecture.id for record in glycoprotein] == ["rvsv_niv_g"]
    protected = select(records, phenotype=Phenotype.PROTECTION, host=Host.NONHUMAN_PRIMATE)
    assert [record.architecture.id for record in protected] == [
        "rvsv_dg_filovirus_gp",
        "rvsv_dg_marv_gp",
        "rvsv_niv_f",
        "rvsv_niv_g",
    ]
    assert protected[0].observations[0].phenotype is Phenotype.PROTECTION
    clinical = select(records, min_evidence=EvidenceLevel.HUMAN)
    assert [record.architecture.id for record in clinical] == [
        "rvsv_dg_sars_cov_2_s",
        "rvsv_lasv_gpc",
    ]
    marburg = select(records, antigen_class="marburgvirus-glycoprotein")
    assert [record.architecture.id for record in marburg] == ["rvsv_dg_marv_gp"]
    filovirus = select(records, antigen_class="filovirus-glycoprotein")
    assert "rvsv_dg_marv_gp" not in [record.architecture.id for record in filovirus]


_DOI = "https://doi.org/10.1128/jvi.78.10.5458-5465.2004"


def _partial() -> Design:
    architecture = Architecture(
        "partial",
        "partial",
        "VSV",
        Replication.UNKNOWN,
        (Leader(), Antigen("antigen_1", "antigen", "filovirus-glycoprotein"), Trailer()),
        parent=ParentStatus.UNKNOWN,
        citations=(_DOI,),
    )
    return Design(
        architecture,
        (
            Observation(
                "safety",
                "paper-2004",
                Phenotype.SAFETY,
                Host.OTHER,
                Assay.SAFETY,
                Outcome.REPORTED,
                _DOI,
            ),
        ),
        (Implementation("paper-2004", Stage.PRECLINICAL, _DOI),),
        unmodeled=(
            Unmodeled("dose", "not quantified", ""),
            Unmodeled("route", "not stated", ""),
            Unmodeled("species", "unspecified", "The abstract names no species."),
        ),
    )


def test_unknown_fields_are_kept(tmp_path: Path) -> None:
    design = tmp_path / "partial"
    design.mkdir()
    document = design_to_document(_partial())
    (design / "design.ttl").write_text(document.write_string(sbol3.TURTLE), encoding="utf-8")
    record = load_design(design)
    kept = {item.key: item for item in record.unmodeled}
    assert kept["dose"].value == "not quantified"
    assert kept["species"].source == "The abstract names no species."
    assert kept["route"].value == "not stated"
    assert record.architecture.replication.value == "unknown"
    assert record.architecture.parent is ParentStatus.UNKNOWN


def test_a_missing_document_is_rejected(tmp_path: Path) -> None:
    design = tmp_path / "bad"
    design.mkdir()
    with pytest.raises(InputError, match=r"missing design\.ttl"):
        load_design(design)


def test_a_sequence_field_is_rejected(tmp_path: Path) -> None:
    design = tmp_path / "bad"
    design.mkdir()
    document = design_to_document(_partial())
    document.add(sbol3.Sequence("bad_seq", elements="atgc", encoding=sbol3.IUPAC_DNA_ENCODING))
    (design / "design.ttl").write_text(document.write_string(sbol3.TURTLE), encoding="utf-8")
    with pytest.raises(InputError, match="sequence payload"):
        load_design(design)


def test_a_sequence_key_inside_unmodeled_is_rejected() -> None:
    record = _partial()
    rejected = Design(
        record.architecture,
        record.observations,
        record.implementations,
        unmodeled=(Unmodeled("sequence", "ATGC", "omitted"),),
    )
    with pytest.raises(InputError, match="sequence payload"):
        design_to_document(rejected)


def test_parameters_describe_the_chassis_and_the_model_validates() -> None:
    architecture = canonical_chassis()
    assert parameters_for(architecture) == {
        "replication": "replication-competent",
        "gene_count": 5,
        "loci": [gene.value for gene in Gene],
    }
    document = to_document(architecture)
    model = attach_model(document, architecture.id)
    assert model.source == "https://doi.org/10.1371/journal.pcbi.1011373"
    assert not document.validate()
