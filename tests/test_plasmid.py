"""A plasmid is a transcription unit on a cited Indiana backbone."""

from pathlib import Path

from lab.protocols import BSAI, Assembly, AssemblyRequest, Part, ProtocolCompiler
from lab.targets import Manual

from tests.synthetic import SyntheticGenome
from vsv import AntigenClass, Gene, Plan, Plasmid, Replication, indiana, read_record
from vsv.design.plasmid import backbone
from vsv.genome.annotation import GENE_END_SIGNAL

_RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "NC_001560.fasta"
_CHARGED = "ATG" + "GAAAAA" * 30 + "TAA"
_GENE_END_ORF = "ATG" + "GCT" * 4 + "TATGAAAAAAA" + "C" + "GCT" * 40 + "TAA"


def _vaccine(
    tmp_path: Path,
    coding: str,
    *,
    replication: Replication = Replication.COMPETENT,
    after: str | None = None,
) -> Plan:
    path = tmp_path / "virus.fasta"
    path.write_text(f">outbreak\n{'A' * 30}{coding}{'A' * 30}\n")
    virus = read_record(path)
    payload = virus.mark(
        virus.add_cds("GP", 30, 30 + len(coding)), AntigenClass.filovirus_glycoprotein
    )
    chassis = indiana()
    if after == "leader":
        return chassis.insert(
            payload, after=chassis.leader, identifier="rvsv_leader_gp", name="leader insert"
        )
    if after == "g":
        return chassis.insert(
            payload,
            after=chassis.feature(Gene.GLYCOPROTEIN),
            identifier="rvsv_outbreak_gp",
            name="insert after G",
        )
    return chassis.replace(
        chassis.feature(Gene.GLYCOPROTEIN),
        payload,
        replication=replication,
        identifier="rvsv_dg_outbreak_gp",
        name="glycoprotein replacement",
    )


def _rules(plasmid: Plasmid) -> set[str]:
    return {violation.rule for violation in plasmid.check().violations}


def test_cited_backbones_match_the_indiana_antigenome() -> None:
    sequence = backbone("pvsv-dg").sequence.bases
    assert len(sequence) == 11161
    replaced = backbone("pvsv-dg")
    assert replaced.accession == "NC_001560.1"
    assert replaced.citation == "https://doi.org/10.1128/jvi.78.10.5458-5465.2004"
    assert (replaced.slot.start, replaced.slot.end) == (3048, 4720)
    assert replaced.slot.upstream.bases == "TATGAAAAAAACT"
    assert replaced.slot.downstream.bases == "CTAACAGCAATC"
    between = backbone("pvsv-g-l")
    assert between.slot.start == between.slot.end == 4720
    assert between.slot.upstream.bases == GENE_END_SIGNAL
    assert between.slot.downstream.bases == "CTAACAGCAATC"
    assert between.citation == "https://doi.org/10.1128/jvi.71.8.5982-5989.1997"
    leader = backbone("pvsv-leader")
    assert leader.slot.start == leader.slot.end == 50
    assert leader.slot.upstream.bases == sequence[30:50]
    assert leader.slot.downstream.bases == "AACAGTAATC"
    assert leader.citation == "https://doi.org/10.1128/jvi.00019-14"


def test_indiana_glycoprotein_passes_the_sequence_rules() -> None:
    virus = read_record(_RAW)
    payload = virus.mark(virus.add_cds("GP", 3077, 4613), AntigenClass.filovirus_glycoprotein)
    chassis = indiana()
    plasmid = chassis.replace(
        chassis.feature(Gene.GLYCOPROTEIN),
        payload,
        replication=Replication.COMPETENT,
        identifier="rvsv_dg_indiana_gp",
        name="Indiana glycoprotein",
    ).plasmid()
    assert plasmid.unit.coding_sequence.bases.startswith("ATG")
    assert plasmid.unit.coding_sequence.bases.endswith("TAA")
    assert len(plasmid.unit.coding_sequence.bases) == 1536
    assert plasmid.check()


def test_a_replacement_antigenome_is_the_slot_with_the_unit_inside(
    synthetic: SyntheticGenome, synthetic_fasta: Path
) -> None:
    virus = read_record(synthetic_fasta)
    start, end = synthetic.cds[Gene.GLYCOPROTEIN]
    payload = virus.mark(virus.add_cds("GP", start, end), AntigenClass.filovirus_glycoprotein)
    chassis = indiana()
    plasmid = chassis.replace(
        chassis.feature(Gene.GLYCOPROTEIN),
        payload,
        replication=Replication.COMPETENT,
        identifier="rvsv_dg_outbreak_gp",
        name="outbreak",
    ).plasmid()
    removed = 4720 - 3048
    assert len(plasmid.antigenome.bases) == 11161 - removed + len(plasmid.sequence.bases)
    slot = plasmid.backbone.slot
    assert (
        plasmid.antigenome.bases[slot.start : slot.start + len(plasmid.sequence.bases)]
        == plasmid.sequence.bases
    )


def test_sequence_rules_reject_an_unclosed_frame_and_an_internal_gene_end(tmp_path: Path) -> None:
    internal_stop = _vaccine(tmp_path, "ATG" + "GCT" * 4 + "TAA" + "GCT" * 4 + "TAA").plasmid()
    assert "open-reading-frame" in _rules(internal_stop)
    assert "transmembrane-anchor" not in _rules(internal_stop)
    no_start = _vaccine(tmp_path, "GCT" * 10 + "TAA").plasmid()
    assert _rules(no_start) == {"open-reading-frame"}
    no_stop = _vaccine(tmp_path, "ATG" + "GCT" * 10).plasmid()
    assert _rules(no_stop) == {"open-reading-frame"}
    internal_end = _vaccine(tmp_path, _GENE_END_ORF).plasmid()
    assert _rules(internal_end) == {"gene-end"}


def test_a_competent_replacement_needs_a_transmembrane_anchor(tmp_path: Path) -> None:
    charged = _vaccine(tmp_path, _CHARGED).plasmid()
    assert _rules(charged) == {"transmembrane-anchor"}
    single_cycle = _vaccine(tmp_path, _CHARGED, replication=Replication.SINGLE_CYCLE).plasmid()
    assert single_cycle.check()
    inserted = _vaccine(tmp_path, _CHARGED, after="g").plasmid()
    assert inserted.check()
    assert inserted.backbone.iri == "https://vsv.bio/backbone/pvsv-g-l"


def test_the_compiler_names_the_plasmid_parts(
    synthetic: SyntheticGenome, synthetic_fasta: Path
) -> None:
    virus = read_record(synthetic_fasta)
    start, end = synthetic.cds[Gene.GLYCOPROTEIN]
    payload = virus.mark(virus.add_cds("GP", start, end), AntigenClass.filovirus_glycoprotein)
    chassis = indiana()
    vaccine = chassis.replace(
        chassis.feature(Gene.GLYCOPROTEIN),
        payload,
        replication=Replication.COMPETENT,
        identifier="rvsv_dg_outbreak_gp",
        name="outbreak",
    )
    plasmid = vaccine.plasmid()
    assert isinstance(plasmid.backbone, Part)
    assert isinstance(plasmid.insert, Part)
    compiled = ProtocolCompiler().compile(
        AssemblyRequest(
            id=vaccine.id,
            assemblies=(
                Assembly(
                    id=f"{vaccine.id}-assembly",
                    product=plasmid.product,
                    backbone=plasmid.backbone,
                    parts=[plasmid.insert],
                    restriction_enzyme=BSAI,
                ),
            ),
        ),
        hardware=Manual(),
    )
    text = "\n".join(compiled.files.values())
    assert "pvsv-dg" in text
    assert "BsaI" in text
    assert "GP" in text
    assert plasmid.product.iri.endswith("/plasmid")
    assert plasmid.sequence.bases not in text
    assert plasmid.sequence.bases.lower() not in text
