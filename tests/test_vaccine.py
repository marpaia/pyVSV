"""Vaccine edits, sequence rules, and the Indiana antigenome cassettes."""

from pathlib import Path

import pytest

from vsv import (
    Antigen,
    Challenge,
    DoseUnit,
    Gene,
    GlycoproteinForm,
    Host,
    Implementation,
    InputError,
    Observation,
    Outcome,
    Phenotype,
    Replace,
    Replication,
    Route,
    Stage,
    Timing,
    Vaccine,
    Viral,
    indiana,
    library,
    read,
)
from vsv.design.antigenome import cassette
from vsv.design.evidence import Assay
from vsv.genome.annotate import GENE_END_SIGNAL
from vsv.genome.sequence import read_fasta

_DOI = "https://doi.org/10.1128/jvi.78.10.5458-5465.2004"
_CHARGED = "ATG" + "GAAAAA" * 30 + "TAA"
_GENE_END = "ATG" + "GCT" * 4 + "TATGAAAAAAA" + "C" + "GCT" * 40 + "TAA"
_HYDROPHOBIC = "CTA" * 20
_EIGHT_A = "ATG" + "GCT" * 4 + "AAAAAAAAC" + "GCT" * 4 + _HYDROPHOBIC + "TAA"
_SEVEN_A = "ATG" + "GCT" * 4 + "AAAAAAACC" + "GCT" * 4 + _HYDROPHOBIC + "TAA"
_OPEN = "ATG" + "GCT" * 20 + "TAA"
_RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "NC_001560.fasta"


def _rules(vaccine: Vaccine) -> set[str]:
    return {item.rule for item in vaccine.check().violations}


def _replacement(
    coding: str,
    *,
    form: GlycoproteinForm = GlycoproteinForm.FULL_LENGTH,
    replication: Replication = Replication.COMPETENT,
    antigen_class: str = "filovirus-glycoprotein",
    after: str | None = None,
) -> Vaccine:
    chassis = indiana()
    if after == "leader":
        return chassis.insert(
            coding,
            after="leader",
            identifier="rvsv_leader_gp",
            name="leader insert",
            antigen_class=antigen_class,
            module="leader_gp",
            form=form,
            strain="Kikwit-1995",
        )
    if after == "g":
        return chassis.insert(
            coding,
            after="g",
            identifier="rvsv_gl_gp",
            name="insert after G",
            antigen_class=antigen_class,
            module="insert_gp",
            form=form,
            strain="Kikwit-1995",
        )
    return chassis.replace(
        coding,
        identifier="rvsv_dg_outbreak_gp",
        name="glycoprotein replacement",
        antigen_class=antigen_class,
        replication=replication,
        form=form,
        strain="Kikwit-1995",
    )


def test_replace_cites_the_method_and_insert_keeps_replication() -> None:
    replaced = _replacement(_EIGHT_A, antigen_class="zaire-ebolavirus-glycoprotein")
    assert replaced.citations == (_DOI,)
    assert replaced.parent == "vsv_indiana"
    assert isinstance(replaced.edit, Replace)
    assert replaced.modules[3].id == "Glycoprotein"
    inserted = indiana().insert(
        _OPEN,
        after="g",
        identifier="rvsv_gl_gp",
        name="insert",
        antigen_class="zaire-ebolavirus-glycoprotein",
        module="zebov_sgp",
        form=GlycoproteinForm.SOLUBLE,
        strain="Mayinga",
    )
    assert inserted.replication is Replication.COMPETENT
    assert inserted.modules[3].id == "Glycoprotein"
    assert inserted.modules[4].id == "zebov_sgp"


def test_a_charged_replacement_needs_a_transmembrane_anchor() -> None:
    charged = _replacement(_CHARGED)
    assert _rules(charged) == {"transmembrane-anchor"}
    chimera = _replacement(_CHARGED, form=GlycoproteinForm.TAIL_CHIMERA)
    assert "transmembrane-anchor" in _rules(chimera)
    soluble = _replacement(_CHARGED, form=GlycoproteinForm.SOLUBLE)
    assert "soluble-entry" in _rules(soluble)
    assert "replication" in _rules(soluble)
    single_cycle = _replacement(_CHARGED, replication=Replication.SINGLE_CYCLE)
    assert "transmembrane-anchor" not in _rules(single_cycle)
    assert "soluble-entry" not in _rules(single_cycle)
    inserted = _replacement(_CHARGED, after="g")
    assert "transmembrane-anchor" not in _rules(inserted)


def test_zaire_editing_site_follows_the_glycoprotein_form() -> None:
    zaire = "zaire-ebolavirus-glycoprotein"
    assert _rules(_replacement(_EIGHT_A, antigen_class=zaire)) == set()
    assert _rules(_replacement(_SEVEN_A, antigen_class=zaire)) == {"editing-site"}
    soluble = _replacement(_SEVEN_A, form=GlycoproteinForm.SOLUBLE, antigen_class=zaire)
    assert "editing-site" not in _rules(soluble)
    assert "soluble-entry" in _rules(soluble)
    poly_k = "ATG" + "AAA" * 70 + "TAA"
    anchor = _rules(_replacement(poly_k, antigen_class=zaire))
    assert {"editing-site", "transmembrane-anchor"} <= anchor
    influenza = _replacement(_SEVEN_A, antigen_class="influenza-hemagglutinin")
    assert "editing-site" not in _rules(influenza)


def test_sequence_rules_reject_an_unclosed_frame_and_an_internal_gene_end() -> None:
    internal_stop = _replacement("ATG" + "GCT" * 4 + "TAA" + "GCT" * 4 + "TAA")
    assert "open-reading-frame" in _rules(internal_stop)
    assert "transmembrane-anchor" not in _rules(internal_stop)
    assert _rules(_replacement("GCT" * 10 + "TAA")) == {"open-reading-frame"}
    assert _rules(_replacement("ATG" + "GCT" * 10)) == {"open-reading-frame"}
    assert _rules(_replacement(_GENE_END)) == {"gene-end"}


def test_cassettes_keep_ct_between_a_gene_end_and_the_next_gene_start() -> None:
    if not _RAW.exists():
        pytest.skip("NC_001560 FASTA is not in data/raw")
    coding = _OPEN
    unit = "AACAGAGATC" + coding + GENE_END_SIGNAL
    replacement = _replacement(coding).antigenome()
    inserted = _replacement(coding, after="g").antigenome()
    leading = _replacement(coding, after="leader").antigenome()
    native = read_fasta(_RAW)[0].sequence
    assert "TATGAAAAAAACT" + unit + "CTAACAGCAATC" in replacement.bases
    assert replacement.cds(Gene.GLYCOPROTEIN) == coding
    assert len(replacement.genes) == 5
    assert GENE_END_SIGNAL + "CT" + unit + "CTAACAGCAATC" in inserted.bases
    assert cassette("g-l", coding) == "CT" + unit
    assert len(inserted.genes) == 6
    assert leading.bases.startswith(native[:47])
    assert cassette("leader", coding) == unit + "CT"
    assert GENE_END_SIGNAL + "CTAACAGTAATC" in leading.bases
    assert len(leading.genes) == 6
    for placed in (replacement, inserted, leading):
        assert "TATGAAAAAAAAC" not in placed.bases


def test_an_unknown_slot_has_no_cassette() -> None:
    with pytest.raises(InputError, match="Unknown slot"):
        cassette("n-p", "ATG" + "GCT" * 4 + "TAA")


def test_precedents_follow_the_published_slot() -> None:
    shelf = library()
    kikwit = next(item for item in shelf if item.id == "rvsv_dg_zebov_kikwit_gp")
    mayinga = next(item for item in shelf if item.id == "rvsv_dg_zebov_mayinga_gp")
    assert "rvsv_dg_zebov_mayinga_gp" in {item.id for item in kikwit.precedents(shelf, locus=3)}
    assert kikwit.precedents(shelf, locus=0) == []
    found = {item.id for item in mayinga.precedents(shelf)}
    assert "rvsv_dg_bombali_gp" in found
    assert "rvsv_dg_zebov_kikwit_gp" in found
    leader = _replacement(_EIGHT_A, after="leader", antigen_class="zaire-ebolavirus-glycoprotein")
    assert "rvsv_dg_zebov_mayinga_gp" not in {item.id for item in leader.precedents(shelf)}


def test_a_vaccine_round_trips_evidence_and_coding(tmp_path: Path) -> None:
    coding = _EIGHT_A
    modules = tuple(
        Antigen(
            "Glycoprotein",
            "Zaire ebolavirus glycoprotein",
            "zaire-ebolavirus-glycoprotein",
            "Kikwit-1995",
            "",
            GlycoproteinForm.FULL_LENGTH,
            coding,
        )
        if gene is Gene.GLYCOPROTEIN
        else Viral(gene)
        for gene in Gene
    )
    observation = Observation(
        "result",
        "lab",
        Phenotype.PROTECTION,
        Host.MOUSE,
        Assay.LETHAL_CHALLENGE,
        Outcome.REPORTED,
        _DOI,
        route=Route.INTRAMUSCULAR,
        timing=Timing.PROPHYLACTIC,
        day=-28,
        vaccine_amount="1e5",
        vaccine_unit=DoseUnit.PFU,
        challenge=Challenge("EBOV", "Mayinga", Route.INTRAPERITONEAL, "1000", DoseUnit.PFU),
    )
    vaccine = Vaccine(
        "example",
        "example vaccine",
        Replication.COMPETENT,
        modules,
        parent="vsv_indiana",
        edit=Replace("Glycoprotein"),
        citations=(_DOI,),
        implementations=(Implementation("lab", Stage.PRECLINICAL, _DOI),),
        observations=(observation,),
    )
    path = tmp_path / "design.ttl"
    vaccine.write(path)
    assert read(path) == vaccine
