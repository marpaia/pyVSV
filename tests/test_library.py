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
    Insert,
    Leader,
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
        "rvsv_dg_bdbv_gp",
        "rvsv_dg_bombali_gp",
        "rvsv_dg_filovirus_gp",
        "rvsv_dg_hcv_e",
        "rvsv_dg_marv_gp",
        "rvsv_dg_mers_s",
        "rvsv_dg_niv_f",
        "rvsv_dg_niv_g",
        "rvsv_dg_sars_s",
        "rvsv_dg_sftsv_gp",
        "rvsv_dg_sudv_gp",
        "rvsv_niv_g",
        "rvsv_niv_f",
        "rvsv_niv_n",
        "rvsv_lasv_gpc",
        "rvsv_dg_sars_cov_2_s",
        "rvsv_hiv_env",
        "rvsv_dg_hiv_env",
        "rvsv_influenza_ha",
        "rvsv_gl_hev_g",
        "rvsv_gl_hiv_gag",
        "rvsv_gl_niv_f",
        "rvsv_gl_niv_g",
        "rvsv_gl_sars_s",
        "rvsv_gfp",
        "rvsv_luc",
        "rvsv_zika_e",
        "rvsv_cchf_gpc",
        "rvsv_rsv_f",
        "rvsv_hantavirus_gp",
        "rvsv_chikv_e",
        "rvsv_rabies_g",
        "rvsv_cat",
        "rvsv_ch_hiv_gag",
        "rvsv_ch_hsv2_gd",
        "rvsv_ch_wnv_e",
        "rvsv_gl_chikv_e",
        "rvsv_gl_crpv_l1",
        "rvsv_gl_denv2_pre",
        "rvsv_gl_hsv2_gd",
        "rvsv_gl_ifnb",
        "rvsv_gl_il4",
        "rvsv_gl_lcrv",
        "rvsv_gl_tk",
        "rvsv_gmcsf",
        "rvsv_hcv_ce1e2",
        "rvsv_leader_gfp",
        "rvsv_leader_lcrv",
        "rvsv_n_zika_prme",
        "rvsv_nj_lcrv",
        "rvsv_wnv_e",
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
    assert records["rvsv_dg_filovirus_gp"].card.evidence_level is EvidenceLevel.HUMAN
    assert "Ervebo" in records["rvsv_dg_filovirus_gp"].aliases

    nipah = records["rvsv_niv_g"].architecture
    added = next(module for module in nipah.modules if module.id == "niv_g")
    assert nipah.parent_id == filo.id
    assert isinstance(nipah.edit, Insert)
    assert nipah.edit.after == "Glycoprotein"
    assert added.locus == 4
    assert next(module.locus for module in nipah.modules if module.id == "LargePolymerase") == 5
    fusion = records["rvsv_niv_f"].architecture
    assert isinstance(fusion.edit, Insert)
    assert next(module for module in fusion.modules if module.id == "niv_f").locus == 4
    assert records["rvsv_niv_f"].card.evidence_level is EvidenceLevel.NONHUMAN_PRIMATE

    lassa = records["rvsv_lasv_gpc"]
    assert [(item.id, item.stage) for item in lassa.implementations] == [
        ("garbutt-2004", Stage.PRECLINICAL),
        ("malkin-2025", Stage.CLINICAL),
        ("pmid-30774999", Stage.PRECLINICAL),
    ]
    assert {item.implementation_id for item in lassa.observations} == {
        "malkin-2025",
        "pmid-30774999",
    }
    assert lassa.card.evidence_level is EvidenceLevel.HUMAN

    marburg = records["rvsv_dg_marv_gp"]
    assert marburg.card.antigen_classes == ("marburgvirus-glycoprotein",)
    assert marburg.card.evidence_level is EvidenceLevel.NONHUMAN_PRIMATE
    assert marburg.observations[0].host is Host.OTHER
    assert "Zaire" not in marburg.observations[0].note
    ebola = records["rvsv_dg_filovirus_gp"].observations[0]
    assert ebola.implementation_id == "garbutt-2004"
    assert "Zaire" in ebola.note
    sudan = records["rvsv_dg_sudv_gp"]
    assert sudan.architecture.antigen_classes == ("sudan-ebolavirus-glycoprotein",)
    assert sudan.observations[0].challenge_virus == "SEBOV"
    assert sudan.observations[0].route.value == "intramuscular"
    hendra = records["rvsv_gl_hev_g"]
    assert hendra.architecture.edit is not None
    assert hendra.observations[0].route.value == "intranasal"
    assert records["rvsv_dg_mers_s"].architecture.replication.value == "replication-competent"
    hiv = records["rvsv_hiv_env"].architecture
    assert isinstance(hiv.edit, Insert)
    assert next(module for module in hiv.modules if module.id == "hiv_env").locus == 4
    single = records["rvsv_dg_hiv_env"]
    assert single.observations[0].host is Host.MOUSE
    assert single.observations[0].route.value == "intramuscular"
    assert records["rvsv_luc"].observations == ()


def test_queries_separate_antigen_class_and_evidence() -> None:
    records = load_library()
    glycoprotein = select(records, antigen_class="henipavirus-glycoprotein")
    assert [record.architecture.id for record in glycoprotein] == [
        "rvsv_dg_niv_g",
        "rvsv_gl_hev_g",
        "rvsv_gl_niv_g",
        "rvsv_niv_g",
    ]
    protected = select(records, phenotype=Phenotype.PROTECTION, host=Host.NONHUMAN_PRIMATE)
    assert [record.architecture.id for record in protected] == [
        "rvsv_dg_bdbv_gp",
        "rvsv_dg_filovirus_gp",
        "rvsv_dg_marv_gp",
        "rvsv_dg_sudv_gp",
        "rvsv_niv_f",
        "rvsv_niv_g",
    ]
    assert protected[0].observations[0].phenotype is Phenotype.PROTECTION
    clinical = select(records, min_evidence=EvidenceLevel.HUMAN)
    assert [record.architecture.id for record in clinical] == [
        "rvsv_dg_filovirus_gp",
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
