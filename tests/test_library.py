"""The curated filovirus shelf."""

from pathlib import Path

from vsv import (
    DoseUnit,
    EvidenceLevel,
    GlycoproteinForm,
    Host,
    Insert,
    Phenotype,
    Replication,
    Route,
    Timing,
    Vaccine,
    library,
    read,
)
from vsv.design.library import select
from vsv.genes import Gene

_SHELF = (
    "rvsv_dg_bdbv_gp",
    "rvsv_dg_bombali_gp",
    "rvsv_dg_marv_gp",
    "rvsv_dg_sudv_gp",
    "rvsv_dg_zebov_kikwit_gp",
    "rvsv_dg_zebov_mayinga_gp",
    "rvsv_gl_zebov_sgp",
    "vsv_indiana",
)


def _by_id(identifier: str) -> Vaccine:
    return next(item for item in library() if item.id == identifier)


def test_the_shelf_is_the_indiana_chassis_and_seven_filovirus_vaccines() -> None:
    shelf = library()
    assert [item.id for item in shelf] == list(_SHELF)
    indiana = _by_id("vsv_indiana")
    assert indiana.parent is None
    assert indiana.edit is None
    assert indiana.realization == "NC_001560"
    assert [module.id for module in indiana.modules] == [gene.value for gene in Gene]


def test_kikwit_is_ervevo_and_records_the_ring_trial() -> None:
    kikwit = _by_id("rvsv_dg_zebov_kikwit_gp")
    assert kikwit.aliases == ("Ervebo",)
    assert kikwit.citations == (
        "https://doi.org/10.1056/NEJMoa1414216",
        "https://doi.org/10.1056/NEJMoa1502924",
        "https://doi.org/10.1016/S1473-3099(15)00154-1",
        "https://doi.org/10.1016/S0140-6736(16)32621-6",
    )
    antigen = kikwit.antigens[0]
    assert antigen.name == "Zaire ebolavirus glycoprotein"
    assert antigen.strain == "Kikwit-1995"
    assert antigen.form is GlycoproteinForm.FULL_LENGTH
    assert [item.id for item in kikwit.implementations] == [
        "regules-2017",
        "agnandji-2016",
        "huttner-2015",
        "henao-restrepo-2017",
    ]
    ring = next(item for item in kikwit.observations if item.id == "henao-restrepo-2017-efficacy")
    assert ring.host is Host.HUMAN
    assert ring.phenotype is Phenotype.PROTECTION
    assert ring.assay.value == "field-efficacy"
    assert ring.route is Route.UNKNOWN
    assert ring.timing is Timing.UNKNOWN
    regules = next(item for item in kikwit.observations if item.id == "regules-2017-adverse")
    assert regules.note.startswith("The most common adverse events were injection-site pain")
    assert "Agnandji" not in regules.note


def test_mayinga_keeps_the_mouse_challenge_and_the_primate_result() -> None:
    mayinga = _by_id("rvsv_dg_zebov_mayinga_gp")
    assert mayinga.antigens[0].strain == "Mayinga"
    assert mayinga.antigens[0].name == "Zaire ebolavirus glycoprotein"
    mouse = next(
        item for item in mayinga.observations if item.id == "garbutt-2004-mouse-protection"
    )
    assert mouse.host is Host.MOUSE
    assert mouse.challenge is not None
    assert mouse.challenge.virus == "mouse-adapted Zaire Ebola virus"
    assert mouse.challenge.route is Route.INTRAPERITONEAL
    assert mouse.challenge.amount == "1,000 LD50"
    assert mouse.challenge.unit is DoseUnit.UNSPECIFIED
    primate = next(item for item in mayinga.observations if item.id == "nm-2005-nhp-protection")
    assert primate.host is Host.NONHUMAN_PRIMATE
    assert primate.challenge is not None
    assert primate.challenge.virus == "EBOV"
    marburg = _by_id("rvsv_dg_marv_gp")
    assert all(item.host is not Host.MOUSE for item in marburg.observations)


def test_soluble_glycoprotein_is_inserted_after_vsv_g() -> None:
    soluble = _by_id("rvsv_gl_zebov_sgp")
    assert isinstance(soluble.edit, Insert)
    assert soluble.edit == Insert("g", "zebov_sgp")
    assert soluble.replication is Replication.COMPETENT
    assert soluble.modules[3].id == "Glycoprotein"
    antigen = soluble.antigens[0]
    assert antigen.id == "zebov_sgp"
    assert antigen.name == "Zaire ebolavirus soluble glycoprotein"
    assert antigen.form is GlycoproteinForm.SOLUBLE
    assert antigen.strain == "Mayinga"
    assert soluble.implementations[0].plasmid == "VSVXN2"
    assert soluble.observations == ()


def test_sudan_is_a_post_exposure_rhesus_result() -> None:
    sudan = _by_id("rvsv_dg_sudv_gp")
    assert sudan.antigens[0].strain == "Boniface"
    assert sudan.antigens[0].name == "Sudan ebolavirus glycoprotein"
    result = next(item for item in sudan.observations if item.id == "geisbert-2008-postexposure")
    assert result.host is Host.RHESUS
    assert result.timing is Timing.POST_EXPOSURE
    assert result.route is Route.INTRAMUSCULAR
    assert result.day is None
    assert result.vaccine_amount == "2e7"
    assert result.vaccine_unit is DoseUnit.PFU
    assert result.challenge is not None
    assert result.challenge.virus == "SEBOV"
    assert result.challenge.route is Route.INTRAMUSCULAR
    assert result.challenge.amount == "1000"
    assert result.challenge.unit is DoseUnit.PFU


def test_blend_and_prime_boost_results_stay_off_the_vaccine_cards() -> None:
    for vaccine in library():
        notes = " ".join(item.note for item in vaccine.observations)
        viruses = [
            item.challenge.virus for item in vaccine.observations if item.challenge is not None
        ]
        assert "CIEBOV" not in viruses
        assert "Group 4" not in notes


def test_queries_follow_antigen_class_host_and_evidence() -> None:
    shelf = library()
    zaire = {item.id for item in select(shelf, antigen_class="zaire-ebolavirus-glycoprotein")}
    assert zaire == {
        "rvsv_dg_zebov_kikwit_gp",
        "rvsv_dg_zebov_mayinga_gp",
        "rvsv_gl_zebov_sgp",
    }
    filovirus = {item.id for item in select(shelf, antigen_class="filovirus-glycoprotein")}
    assert filovirus == set(_SHELF) - {"vsv_indiana"}
    protected = {
        item.id
        for item in select(shelf, phenotype=Phenotype.PROTECTION, host=Host.NONHUMAN_PRIMATE)
    }
    assert protected == {
        "rvsv_dg_bdbv_gp",
        "rvsv_dg_marv_gp",
        "rvsv_dg_sudv_gp",
        "rvsv_dg_zebov_mayinga_gp",
    }
    human = {item.id for item in select(shelf, min_evidence=EvidenceLevel.HUMAN)}
    assert human == {"rvsv_dg_zebov_kikwit_gp"}


def test_bombali_only_gap_is_unknown_replication() -> None:
    bombali = _by_id("rvsv_dg_bombali_gp")
    assert bombali.replication is Replication.UNKNOWN
    assert {item.rule for item in bombali.check().violations} == {"replication-unknown"}
    assert not bombali.check()


def test_shelf_cards_round_trip(tmp_path: Path) -> None:
    for vaccine in library():
        path = tmp_path / f"{vaccine.id}.ttl"
        vaccine.write(path)
        assert read(path) == vaccine
