"""The curated design library, its queries, and the linked intracellular model."""

import json
from pathlib import Path

import pytest

from vsv.architecture import (
    LocusStatus,
    Module,
    ModuleKind,
    ParentStatus,
    canonical_chassis,
    to_document,
)
from vsv.evidence import EvidenceLevel, Phenotype, Stage
from vsv.genes import Gene
from vsv.library import load_design, load_library
from vsv.model import attach_model, parameters_for
from vsv.query import select
from vsv.rules import Host, check
from vsv.sequence import InputError
from vsv.transform import ReplaceModule


def test_seed_designs_match_their_recorded_transformations() -> None:
    records = {record.architecture.id: record for record in load_library()}
    assert set(records) == {
        "vsv_indiana",
        "rvsv_dg_filovirus_gp",
        "rvsv_dg_marv_gp",
        "rvsv_niv_g",
        "rvsv_niv_f",
        "rvsv_lasv_gpc",
    }
    indiana = records["vsv_indiana"]
    assert indiana.realization == "NC_001560"
    assert indiana.model_source == "https://doi.org/10.1371/journal.pcbi.1011373"
    assert indiana.architecture.parent_status is ParentStatus.ROOT
    assert indiana.implementations == ()
    assert check(indiana.architecture, indiana.observations) == []
    assert indiana.card.evidence_level is EvidenceLevel.NONE

    filo = records["rvsv_dg_filovirus_gp"].architecture
    applied = ReplaceModule(
        Gene.GLYCOPROTEIN,
        Module(
            Gene.GLYCOPROTEIN,
            ModuleKind.ANTIGEN,
            "filovirus glycoprotein",
            antigen_class="filovirus-glycoprotein",
        ),
        filo.id,
        name=filo.name,
    ).apply(indiana.architecture)
    assert filo.modules == applied.modules
    assert filo.derivation == applied.derivation
    assert filo.replication == applied.replication
    assert records["rvsv_dg_filovirus_gp"].card.evidence_level is EvidenceLevel.RODENT

    nipah = records["rvsv_niv_g"].architecture
    added = next(module for module in nipah.modules if module.id == "niv_g")
    assert nipah.parent_id == filo.id
    assert nipah.derivation is None
    assert added.locus is None
    assert added.locus_status is LocusStatus.UNKNOWN
    assert next(module.locus for module in nipah.modules if module.id == "LargePolymerase") == 4
    fusion = records["rvsv_niv_f"].architecture
    assert fusion.derivation is None
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
    assert marburg.card.evidence_level is EvidenceLevel.NONE
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
    assert [record.architecture.id for record in protected] == ["rvsv_niv_f", "rvsv_niv_g"]
    assert protected[0].observations[0].phenotype is Phenotype.PROTECTION
    clinical = select(records, min_evidence=EvidenceLevel.HUMAN)
    assert [record.architecture.id for record in clinical] == ["rvsv_lasv_gpc"]
    marburg = select(records, antigen_class="marburgvirus-glycoprotein")
    assert [record.architecture.id for record in marburg] == ["rvsv_dg_marv_gp"]
    filovirus = select(records, antigen_class="filovirus-glycoprotein")
    assert "rvsv_dg_marv_gp" not in [record.architecture.id for record in filovirus]


def _partial_design(extra: dict[str, object]) -> dict[str, object]:
    return {
        "id": "partial",
        "name": "partial",
        "backbone": "VSV",
        "replication": "unknown",
        "parent_status": "unknown",
        "citations": ["https://doi.org/10.1128/jvi.78.10.5458-5465.2004"],
        "modules": [
            {"id": "leader", "kind": "leader", "name": "leader"},
            {
                "id": "antigen_1",
                "kind": "antigen",
                "name": "antigen",
                "antigen_class": "filovirus-glycoprotein",
                "locus_status": "unknown",
            },
            {"id": "trailer", "kind": "trailer", "name": "trailer"},
        ],
        "implementations": [
            {
                "id": "paper-2004",
                "stage": "preclinical",
                "citation": "https://doi.org/10.1128/jvi.78.10.5458-5465.2004",
            }
        ],
        **extra,
    }


def test_unknown_fields_are_kept(tmp_path: Path) -> None:
    design = tmp_path / "partial"
    design.mkdir()
    (design / "design.json").write_text(
        json.dumps(
            _partial_design(
                {
                    "dose": "not quantified",
                    "unmodeled": [
                        {
                            "key": "species",
                            "value": "unspecified",
                            "source": "The abstract names no species.",
                        }
                    ],
                }
            )
        )
    )
    (design / "observations.json").write_text(
        json.dumps(
            [
                {
                    "id": "safety",
                    "implementation_id": "paper-2004",
                    "phenotype": "safety",
                    "host": "other",
                    "assay": "safety",
                    "outcome": "reported",
                    "citation": "https://doi.org/10.1128/jvi.78.10.5458-5465.2004",
                    "route": "not stated",
                }
            ]
        )
    )
    record = load_design(design)
    kept = {item.key: item for item in record.unmodeled}
    assert kept["dose"].value == "not quantified"
    assert kept["species"].source == "The abstract names no species."
    assert kept["route"].value == "not stated"
    assert record.architecture.replication.value == "unknown"
    assert record.architecture.parent_status is ParentStatus.UNKNOWN


def test_a_sequence_field_is_rejected(tmp_path: Path) -> None:
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
                "sequence": "ATGC",
            }
        )
    )
    with pytest.raises(InputError, match="sequence payload"):
        load_design(design)


def test_a_sequence_key_inside_unmodeled_is_rejected(tmp_path: Path) -> None:
    design = tmp_path / "bad"
    design.mkdir()
    (design / "design.json").write_text(
        json.dumps(
            _partial_design(
                {"unmodeled": [{"key": "sequence", "value": "ATGC", "source": "omitted"}]}
            )
        )
    )
    with pytest.raises(InputError, match="sequence payload"):
        load_design(design)


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
