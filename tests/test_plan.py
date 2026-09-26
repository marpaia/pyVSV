"""Building a vaccine plan writes SBOL and leaves the curated library alone."""

from pathlib import Path

import pytest
import sbol3

from tests.synthetic import SyntheticGenome
from vsv import (
    AntigenClass,
    Gene,
    Phenotype,
    Replication,
    indiana,
    load,
    plan,
    plans,
    read,
    read_plasmid,
    read_record,
)
from vsv.corpus.extract import _REPORTER_SLUGS, _SLUGS
from vsv.design.document import architecture_of
from vsv.design.edits import Insert, Replace
from vsv.errors import InputError
from vsv.genome.payload import AntigenFeature
from vsv.genome.sequence import Nucleotides

_REPLACE = "https://doi.org/10.1128/jvi.78.10.5458-5465.2004"
_AFTER_G = "https://doi.org/10.1128/jvi.71.8.5982-5989.1997"
_AFTER_LEADER = "https://doi.org/10.1128/jvi.00019-14"
_REPLACEMENTS = (
    "rvsv_cchf_gpc",
    "rvsv_chikv_e",
    "rvsv_dg_filovirus_gp",
    "rvsv_dg_hiv_env",
    "rvsv_dg_marv_gp",
    "rvsv_dg_sars_cov_2_s",
    "rvsv_hantavirus_gp",
    "rvsv_lasv_gpc",
    "rvsv_rabies_g",
)
_AFTER_GLYCOPROTEIN = (
    "rvsv_gfp",
    "rvsv_influenza_ha",
    "rvsv_rsv_f",
    "rvsv_zika_e",
)
_UNPOSITIONED = ("rvsv_hiv_env", "rvsv_niv_f", "rvsv_niv_g", "vsv_indiana")


def _payload(synthetic: SyntheticGenome, fasta: Path) -> tuple[AntigenFeature, str]:
    virus = read_record(fasta)
    start, end = synthetic.cds[Gene.GLYCOPROTEIN]
    coding = virus.add_cds("GP", start, end)
    payload = virus.mark(coding, AntigenClass.filovirus_glycoprotein)
    return payload, virus.feature_bases(coding)


def test_nucleotides_require_uppercase_iupac_dna() -> None:
    assert Nucleotides("ATGC").bases == "ATGC"
    with pytest.raises(InputError, match="uppercase IUPAC"):
        Nucleotides("atgc")


def test_read_record_does_not_name_vsv_genes(
    synthetic_fasta: Path, synthetic_genbank: Path
) -> None:
    fasta = read_record(synthetic_fasta)
    assert fasta.genes == {}
    assert set(load(synthetic_fasta).genes) == set(Gene)

    genbank = read_record(synthetic_genbank)
    assert genbank.genes == {}
    assert set(load(synthetic_genbank).genes) == set(Gene)


def test_mark_records_the_antigen_class_on_the_coding_feature(
    synthetic: SyntheticGenome, synthetic_fasta: Path
) -> None:
    payload, bases = _payload(synthetic, synthetic_fasta)
    assert payload.sequence.bases == bases
    assert payload.antigen_class is AntigenClass.filovirus_glycoprotein
    assert payload.antigen_class.role in {str(role) for role in payload.feature.roles}
    fresh = AntigenClass("paramyxovirus-glycoprotein")
    assert fresh.role == "https://vsv.bio/antigen-class/paramyxovirus-glycoprotein"
    declared = {
        value.slug for value in vars(AntigenClass).values() if isinstance(value, AntigenClass)
    }
    assert declared == {slug for _pattern, slug, _name in _SLUGS if slug not in _REPORTER_SLUGS}


def test_replace_writes_a_sequence_free_plan(
    synthetic: SyntheticGenome, synthetic_fasta: Path
) -> None:
    payload, bases = _payload(synthetic, synthetic_fasta)
    chassis = indiana()
    vaccine = chassis.replace(
        chassis.feature(Gene.GLYCOPROTEIN),
        payload,
        replication=Replication.COMPETENT,
        identifier="rvsv_dg_outbreak_gp",
        name="rVSV outbreak glycoprotein",
    )
    architecture = architecture_of(vaccine.document)
    assert [(module.id, module.locus) for module in architecture.modules] == [
        ("leader", None),
        ("Nucleoprotein", 0),
        ("Phosphoprotein", 1),
        ("MatrixProtein", 2),
        ("GP", 3),
        ("LargePolymerase", 4),
        ("trailer", None),
    ]
    assert architecture.edit == Replace("Glycoprotein", "GP", Replication.COMPETENT)
    assert architecture.citations == (_REPLACE,)
    assert vaccine.check()
    assert vaccine.observations == ()
    assert not any(isinstance(obj, sbol3.Sequence) for obj in vaccine.document.objects)
    assert bases not in vaccine.document.write_string(sbol3.TURTLE)
    assert chassis.id == "vsv_indiana"
    assert Gene.GLYCOPROTEIN in {module.id for module in architecture_of(chassis.document).modules}

    single_cycle = chassis.replace(
        chassis.feature(Gene.GLYCOPROTEIN),
        payload,
        replication=Replication.SINGLE_CYCLE,
        identifier="rvsv_dg_outbreak_gp_sc",
        name="single-cycle rVSV outbreak glycoprotein",
    )
    assert single_cycle.check()
    assert architecture_of(single_cycle.document).replication is Replication.SINGLE_CYCLE


def test_insert_after_glycoprotein_keeps_viral_g(
    synthetic: SyntheticGenome, synthetic_fasta: Path
) -> None:
    payload, _bases = _payload(synthetic, synthetic_fasta)
    chassis = indiana()
    vaccine = chassis.insert(
        payload,
        after=chassis.feature(Gene.GLYCOPROTEIN),
        identifier="rvsv_outbreak_gp",
        name="rVSV outbreak glycoprotein",
    )
    architecture = architecture_of(vaccine.document)
    assert [(module.id, module.locus) for module in architecture.modules] == [
        ("leader", None),
        ("Nucleoprotein", 0),
        ("Phosphoprotein", 1),
        ("MatrixProtein", 2),
        ("Glycoprotein", 3),
        ("GP", 4),
        ("LargePolymerase", 5),
        ("trailer", None),
    ]
    assert architecture.edit == Insert("Glycoprotein", "GP")
    assert architecture.citations == (_AFTER_G,)
    assert vaccine.check()


def test_insert_after_the_leader(synthetic: SyntheticGenome, synthetic_fasta: Path) -> None:
    payload, _bases = _payload(synthetic, synthetic_fasta)
    chassis = indiana()
    vaccine = chassis.insert(
        payload,
        after=chassis.leader,
        identifier="rvsv_outbreak_gp",
        name="rVSV outbreak glycoprotein",
    )
    architecture = architecture_of(vaccine.document)
    assert architecture.modules[1].id == "GP"
    assert architecture.modules[1].locus == 0
    assert architecture.modules[2].id == "Nucleoprotein"
    assert architecture.citations == (_AFTER_LEADER,)
    assert vaccine.check()


def test_precedents_follow_the_precedence_chain(
    synthetic: SyntheticGenome, synthetic_fasta: Path
) -> None:
    payload, _bases = _payload(synthetic, synthetic_fasta)
    library = plans()
    chassis = indiana()
    replaced = chassis.replace(
        chassis.feature(Gene.GLYCOPROTEIN),
        payload,
        replication=Replication.COMPETENT,
        identifier="rvsv_dg_outbreak_gp",
        name="rVSV outbreak glycoprotein",
    )
    inserted = chassis.insert(
        payload,
        after=chassis.feature(Gene.GLYCOPROTEIN),
        identifier="rvsv_outbreak_gp",
        name="rVSV outbreak glycoprotein",
    )
    leading = chassis.insert(
        payload,
        after=chassis.leader,
        identifier="rvsv_leader_gp",
        name="rVSV outbreak glycoprotein",
    )
    assert tuple(item.id for item in replaced.precedents(library)) == _REPLACEMENTS
    assert tuple(item.id for item in inserted.precedents(library)) == _AFTER_GLYCOPROTEIN
    assert tuple(item.id for item in leading.precedents(library)) == ("rvsv_luc",)
    for omitted in _UNPOSITIONED:
        assert omitted not in {item.id for item in replaced.precedents(library)}
    filo = next(item for item in library if item.id == "rvsv_dg_filovirus_gp")
    assert filo.observations
    assert filo.check()
    nipah = next(item for item in library if item.id == "rvsv_niv_g")
    assert not nipah.check()
    assert nipah.check().violations


def test_replace_rejects_a_glycoprotein_that_is_already_an_antigen(
    synthetic: SyntheticGenome, synthetic_fasta: Path
) -> None:
    payload, _bases = _payload(synthetic, synthetic_fasta)
    filo = next(item for item in plans() if item.id == "rvsv_dg_filovirus_gp")
    with pytest.raises(InputError, match="not a viral gene"):
        filo.replace(
            filo.feature(Gene.GLYCOPROTEIN),
            payload,
            replication=Replication.COMPETENT,
            identifier="child",
            name="child",
        )


def test_plan_round_trips_through_an_sbol_file(
    synthetic: SyntheticGenome, synthetic_fasta: Path, tmp_path: Path
) -> None:
    payload, _bases = _payload(synthetic, synthetic_fasta)
    chassis = indiana()
    vaccine = chassis.replace(
        chassis.feature(Gene.GLYCOPROTEIN),
        payload,
        replication=Replication.COMPETENT,
        identifier="rvsv_dg_outbreak_gp",
        name="rVSV outbreak glycoprotein",
    )
    path = tmp_path / "rvsv_dg_outbreak_gp.ttl"
    vaccine.write(path)
    loaded = read(path)
    assert loaded.antigen(AntigenClass.filovirus_glycoprotein).display_id == "GP"
    assert loaded.check()
    with pytest.raises(InputError, match="expected one marburgvirus-glycoprotein"):
        loaded.antigen(AntigenClass.marburgvirus_glycoprotein)


def test_plasmid_uses_the_sequence_recorded_at_mark(
    synthetic: SyntheticGenome, synthetic_fasta: Path, tmp_path: Path
) -> None:
    payload, bases = _payload(synthetic, synthetic_fasta)
    chassis = indiana()
    vaccine = chassis.replace(
        chassis.feature(Gene.GLYCOPROTEIN),
        payload,
        replication=Replication.COMPETENT,
        identifier="rvsv_dg_outbreak_gp",
        name="rVSV outbreak glycoprotein",
    )
    plasmid = vaccine.plasmid()
    assert plasmid.unit.coding_sequence.bases == bases
    assert plasmid.unit.gene_start.bases == "AACAGAGATC"
    assert plasmid.unit.gene_end.bases == "TATGAAAAAAA"
    assert plasmid.sequence.bases == f"AACAGAGATC{bases}TATGAAAAAAA"
    assert plasmid.check()
    assert plasmid.slot == "g"
    assert plasmid.edit == "replace g"
    assert plasmid.insert.iri == "https://vsv.bio/rvsv_dg_outbreak_gp/GP"
    assert plasmid.backbone.iri == "https://vsv.bio/backbone/pvsv-dg"
    assert plasmid.product.iri == "https://vsv.bio/rvsv_dg_outbreak_gp/plasmid"
    plan_path = tmp_path / "rvsv_dg_outbreak_gp.ttl"
    vaccine.write(plan_path)
    assert bases not in plan_path.read_text(encoding="utf-8")
    loaded = read(plan_path)
    with pytest.raises(InputError, match="no coding sequence"):
        loaded.plasmid()
    plasmid_path = tmp_path / "rvsv_dg_outbreak_gp_plasmid.ttl"
    plasmid.write(plasmid_path)
    restored = read_plasmid(plasmid_path)
    assert restored == plasmid
    assert bases.lower() in plasmid_path.read_text(encoding="utf-8").lower()


def test_plasmid_follows_the_insertion_site(
    synthetic: SyntheticGenome, synthetic_fasta: Path
) -> None:
    payload, _bases = _payload(synthetic, synthetic_fasta)
    chassis = indiana()
    after_g = chassis.insert(
        payload,
        after=chassis.feature(Gene.GLYCOPROTEIN),
        identifier="rvsv_outbreak_gp",
        name="rVSV outbreak glycoprotein",
    )
    after_leader = chassis.insert(
        payload,
        after=chassis.leader,
        identifier="rvsv_leader_gp",
        name="rVSV outbreak glycoprotein",
    )
    assert after_g.plasmid().slot == "g-l"
    assert after_g.plasmid().backbone.iri == "https://vsv.bio/backbone/pvsv-g-l"
    assert after_leader.plasmid().slot == "leader"
    assert after_leader.plasmid().backbone.iri == "https://vsv.bio/backbone/pvsv-leader"
    with pytest.raises(InputError, match="no coding sequence"):
        plans()[0].plasmid()


def test_precedents_can_require_an_antigen_class_and_a_phenotype(
    synthetic: SyntheticGenome, synthetic_fasta: Path
) -> None:
    payload, _bases = _payload(synthetic, synthetic_fasta)
    library = plans()
    chassis = indiana()
    replaced = chassis.replace(
        chassis.feature(Gene.GLYCOPROTEIN),
        payload,
        replication=Replication.COMPETENT,
        identifier="rvsv_dg_outbreak_gp",
        name="rVSV outbreak glycoprotein",
    )
    matched = replaced.precedents(
        library,
        antigen_class=AntigenClass.filovirus_glycoprotein,
        phenotype=Phenotype.PROTECTION,
    )
    assert tuple(item.id for item in matched) == ("rvsv_dg_filovirus_gp",)
    sudan = replaced.precedents(library, antigen_class=AntigenClass.sudan_ebolavirus_glycoprotein)
    assert sudan == []


def test_plan_of_a_sequenced_genome_has_no_sequence(synthetic_fasta: Path) -> None:
    genome = load(synthetic_fasta)
    designed = plan(genome)
    assert designed.id == "SYN0001"
    assert not any(isinstance(obj, sbol3.Sequence) for obj in designed.document.objects)
    assert genome.bases.lower() not in designed.document.write_string(sbol3.TURTLE).lower()
