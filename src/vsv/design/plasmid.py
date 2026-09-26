"""A plasmid for one vaccine plan.

The plan file has no coding sequence. The plasmid file carries a transcription unit on one insert
and the cited Indiana antigenome the unit replaces or is inserted into.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import sbol3
from lab.protocols import Part

from vsv import NAMESPACE
from vsv.design.architecture import Replication
from vsv.design.rules import CheckResult, Violation
from vsv.errors import InputError
from vsv.genes import Gene
from vsv.genome.annotation import GENE_END_SIGNAL
from vsv.genome.sbol import (
    ROLE_GENE_END,
    ROLE_GENE_START,
    SO_CDS,
    SO_GENE,
    SO_REGION,
    SO_REGULATORY_REGION,
    SO_TERMINATOR,
    new_document,
    role_key,
)
from vsv.genome.sequence import (
    STOP_CODONS,
    Nucleotides,
    hydrophobic_window,
    read_fasta,
    translate,
)

_ACCESSION = "NC_001560.1"
_GENE_START = "AACAGAGATC"
_START_PATTERN = re.compile(r"AACAG[ACGT]{2}ATC\Z")
_SLOT_NAMES = {"pvsv-dg": "g", "pvsv-g-l": "g-l", "pvsv-leader": "leader"}
_G_SLOT = (3048, 4720)
_GL_CUT = 4720
_LEADER_CUT = 50
_LEADER_FLANK = 20
_PRODUCT_KEYS = ("plan", "edit", "slot", "insert", "replication", "accession")
_BACKBONE_KEYS = (
    "accession",
    "citation",
    "slot-start",
    "slot-end",
    "upstream",
    "downstream",
    "gene-start",
    "gene-end",
)
_REPLACE_G = f"replace {Gene.GLYCOPROTEIN.short}"


@dataclass(frozen=True, slots=True)
class TranscriptionUnit:
    """Antigenome-sense transcription unit: gene-start, coding sequence, and gene-end."""

    gene_start: Nucleotides
    coding_sequence: Nucleotides
    gene_end: Nucleotides

    def __post_init__(self) -> None:
        if _START_PATTERN.fullmatch(self.gene_start.bases) is None:
            raise InputError("A gene-start is AACAG, two bases, and ATC.")
        if self.gene_end.bases != GENE_END_SIGNAL:
            raise InputError("A gene-end is TATG followed by seven A's.")

    @property
    def sequence(self) -> Nucleotides:
        """Gene-start, coding sequence, and gene-end, with nothing between them."""
        return Nucleotides(self.gene_start.bases + self.coding_sequence.bases + self.gene_end.bases)


@dataclass(frozen=True, slots=True)
class Slot:
    """Half-open span of a backbone sequence that a transcription unit replaces.

    An empty span (``start == end``) inserts. ``upstream`` ends at ``start``. ``downstream``
    begins at ``end``.
    """

    name: str
    start: int
    end: int
    upstream: Nucleotides
    downstream: Nucleotides
    gene_start: Nucleotides
    gene_end: Nucleotides


@dataclass(frozen=True, slots=True)
class Backbone(Part):
    """A ``Part`` the Lab compiler can take as ``backbone=``, plus the cited sequence."""

    id: str
    accession: str
    citation: str
    sequence: Nucleotides
    slot: Slot

    def __post_init__(self) -> None:
        # This slots dataclass does not bind super() to the instance class.
        Part.__post_init__(self)
        if self.iri != f"{NAMESPACE}/backbone/{self.id}":
            raise InputError(f"{self.id}: the backbone IRI does not match its id.")
        expected = _SLOT_NAMES.get(self.id)
        if expected is not None and self.slot.name != expected:
            raise InputError(f"{self.id}: the slot is {self.slot.name}.")
        bases = self.sequence.bases
        slot = self.slot
        if not 0 <= slot.start <= slot.end <= len(bases):
            raise InputError(f"{self.id}: the slot is outside the backbone sequence.")
        upstream = slot.upstream.bases
        downstream = slot.downstream.bases
        if slot.start < len(upstream) or bases[slot.start - len(upstream) : slot.start] != upstream:
            raise InputError(f"{self.id}: the upstream flank is not on the backbone sequence.")
        end = slot.end + len(downstream)
        if end > len(bases) or bases[slot.end : end] != downstream:
            raise InputError(f"{self.id}: the downstream flank is not on the backbone sequence.")


@dataclass(frozen=True, slots=True)
class Plasmid:
    """The transcription unit, cited backbone, and product for one plan."""

    plan_id: str
    edit: str
    slot: str
    insert_id: str
    insert: Part
    unit: TranscriptionUnit
    sequence: Nucleotides
    replication: Replication
    backbone: Backbone
    antigenome: Nucleotides
    product: Part

    def __post_init__(self) -> None:
        if self.sequence.bases != self.unit.sequence.bases:
            raise InputError(f"{self.plan_id}: the plasmid sequence is not the transcription unit.")
        if self.slot != self.backbone.slot.name:
            raise InputError(f"{self.plan_id}: the slot does not match the backbone.")
        if self.unit.gene_start != self.backbone.slot.gene_start:
            raise InputError(f"{self.plan_id}: the gene-start does not match the slot.")
        if self.unit.gene_end != self.backbone.slot.gene_end:
            raise InputError(f"{self.plan_id}: the gene-end does not match the slot.")
        if self.antigenome.bases != _antigenome(self.backbone, self.unit).bases:
            raise InputError(f"{self.plan_id}: the antigenome does not place the unit in the slot.")

    def check(self) -> CheckResult:
        """Sequence rules for this transcription unit. True when every rule passes."""
        violations: list[Violation] = []
        coding = self.unit.coding_sequence.bases
        closed = _closed_orf(coding)
        if not closed:
            violations.append(
                Violation(
                    "open-reading-frame",
                    f"{self.plan_id}: the coding sequence is not a closed open reading frame.",
                )
            )
        if not _terminal_gene_end(self.sequence.bases):
            violations.append(
                Violation(
                    "gene-end",
                    f"{self.plan_id}: a gene-end signal occurs inside the transcription unit.",
                )
            )
        if closed and self.replication is Replication.COMPETENT and self.edit == _REPLACE_G:
            protein = translate(coding, to_stop=True)
            if not hydrophobic_window(protein):
                violations.append(
                    Violation(
                        "transmembrane-anchor",
                        f"{self.plan_id}: a replication-competent glycoprotein replacement "
                        "has no transmembrane anchor.",
                    )
                )
        return CheckResult(tuple(violations))

    def write(self, path: str | Path) -> None:
        """Write this plasmid as an SBOL3 file. The suffix selects the format."""
        _document(self).write(path)


def read_plasmid(path: str | Path) -> Plasmid:
    """Read the plasmid in an SBOL3 file."""
    document = sbol3.Document()
    try:
        document.read(str(path))
    except (OSError, ValueError) as error:
        raise InputError(f"{path}: {error}") from error
    products = [
        obj
        for obj in document.objects
        if isinstance(obj, sbol3.Component) and str(obj.identity).endswith("/plasmid")
    ]
    if len(products) != 1:
        raise InputError(f"{path}: expected one plasmid, found {len(products)}.")
    product = products[0]
    facts = _require(_facts(str(product.description or "")), _PRODUCT_KEYS, path)
    roles = [str(role) for role in product.roles if str(role).startswith(f"{NAMESPACE}/backbone/")]
    if len(roles) != 1:
        raise InputError(f"{path}: expected one backbone, found {len(roles)}.")
    located = [feature for feature in product.features if isinstance(feature, sbol3.SubComponent)]
    if len(located) != 1:
        raise InputError(f"{path}: expected one insert, found {len(located)}.")
    insert = document.find(str(located[0].instance_of))
    if not isinstance(insert, sbol3.Component):
        raise InputError(f"{path}: the insert component is missing.")
    # The part IRI ends in a hyphenated id, which SBOL rejects as a display id.
    # The component uses a display id; the product role keeps the part IRI.
    cited_components = [
        obj
        for obj in document.objects
        if isinstance(obj, sbol3.Component) and "slot-start=" in str(obj.description or "")
    ]
    if len(cited_components) != 1:
        raise InputError(f"{path}: expected one backbone, found {len(cited_components)}.")
    backbone_component = cited_components[0]
    record = _require(_facts(str(backbone_component.description or "")), _BACKBONE_KEYS, path)
    if record["accession"] != facts["accession"]:
        raise InputError(f"{path}: the backbone accession does not match the plasmid.")
    sequence = Nucleotides(_bases(document, backbone_component, path))
    slot = Slot(
        facts["slot"],
        _coordinate(record["slot-start"], path),
        _coordinate(record["slot-end"], path),
        Nucleotides(record["upstream"]),
        Nucleotides(record["downstream"]),
        Nucleotides(record["gene-start"]),
        Nucleotides(record["gene-end"]),
    )
    _check_slot_feature(backbone_component, slot, path)
    cited = Backbone(
        roles[0],
        roles[0].rsplit("/", 1)[-1],
        record["accession"],
        record["citation"],
        sequence,
        slot,
    )
    unit = _unit(document, insert, path)
    try:
        replication = Replication(facts["replication"])
    except ValueError as error:
        raise InputError(f"{path}: unknown replication {facts['replication']}.") from error
    return Plasmid(
        facts["plan"],
        facts["edit"],
        facts["slot"],
        facts["insert"],
        Part(str(insert.identity)),
        unit,
        unit.sequence,
        replication,
        cited,
        _antigenome(cited, unit),
        Part(str(product.identity)),
    )


def backbone(identifier: str) -> Backbone:
    """The cited backbone ``identifier``: ``pvsv-dg``, ``pvsv-g-l``, or ``pvsv-leader``."""
    try:
        return _backbones()[identifier]
    except KeyError:
        raise InputError(f"Unknown backbone {identifier}.") from None


def _antigenome(cited: Backbone, unit: TranscriptionUnit) -> Nucleotides:
    bases = cited.sequence.bases
    slot = cited.slot
    return Nucleotides(bases[: slot.start] + unit.sequence.bases + bases[slot.end :])


def _closed_orf(coding: str) -> bool:
    if len(coding) < 6 or len(coding) % 3 != 0 or not coding.startswith("ATG"):
        return False
    if coding[-3:] not in STOP_CODONS:
        return False
    return all(
        coding[index : index + 3] not in STOP_CODONS for index in range(0, len(coding) - 3, 3)
    )


def _terminal_gene_end(unit: str) -> bool:
    return unit.endswith(GENE_END_SIGNAL) and unit.count(GENE_END_SIGNAL) == 1


@cache
def _backbones() -> dict[str, Backbone]:
    bases = _indiana_sequence().bases
    _require_accession(bases)
    signals = (Nucleotides(_GENE_START), Nucleotides(GENE_END_SIGNAL))
    shared = {
        "pvsv-dg": Slot(
            "g", *_G_SLOT, _flank(bases, 3035, 3048), _flank(bases, 4720, 4732), *signals
        ),
        "pvsv-g-l": Slot(
            "g-l", _GL_CUT, _GL_CUT, _flank(bases, 4709, 4720), _flank(bases, 4720, 4732), *signals
        ),
        "pvsv-leader": Slot(
            "leader",
            _LEADER_CUT,
            _LEADER_CUT,
            _flank(bases, _LEADER_CUT - _LEADER_FLANK, _LEADER_CUT),
            _flank(bases, _LEADER_CUT, _LEADER_CUT + len(_GENE_START)),
            *signals,
        ),
    }
    citations = {
        "pvsv-dg": "https://doi.org/10.1128/jvi.78.10.5458-5465.2004",
        "pvsv-g-l": "https://doi.org/10.1128/jvi.71.8.5982-5989.1997",
        "pvsv-leader": "https://doi.org/10.1128/jvi.00019-14",
    }
    return {
        identifier: Backbone(
            f"{NAMESPACE}/backbone/{identifier}",
            identifier,
            _ACCESSION,
            citations[identifier],
            Nucleotides(bases),
            slot,
        )
        for identifier, slot in shared.items()
    }


def _flank(bases: str, start: int, end: int) -> Nucleotides:
    return Nucleotides(bases[start:end])


def _require_accession(bases: str) -> None:
    expected = {
        (3035, 3048): "TATGAAAAAAACT",
        (3048, 3058): _GENE_START,
        (4709, 4720): GENE_END_SIGNAL,
        (4720, 4732): "CTAACAGCAATC",
        (50, 60): "AACAGTAATC",
    }
    for (start, end), sequence in expected.items():
        if bases[start:end] != sequence:
            raise InputError(f"{_ACCESSION}: the cited signal at {start} is not on this sequence.")
    if len(bases) != 11161:
        raise InputError(f"{_ACCESSION}: expected 11161 nt, found {len(bases)}.")


def _indiana_sequence() -> Nucleotides:
    path = _antigenome_path()
    try:
        records = read_fasta(path)
    except InputError as error:
        raise InputError(f"{path}: {error}") from error
    if len(records) != 1 or records[0].id.split()[0] != _ACCESSION:
        raise InputError(f"{path}: expected one {_ACCESSION} record.")
    return Nucleotides(records[0].sequence)


def _antigenome_path() -> Path:
    root = Path(__file__).resolve().parents[3]
    checkout = root / "data" / "raw" / "NC_001560.fasta"
    if (root / "pyproject.toml").is_file() and checkout.is_file():
        return checkout
    packaged = Path(__file__).resolve().parents[1] / "NC_001560.fasta"
    if packaged.is_file():
        return packaged
    raise InputError("The Indiana antigenome is not installed.")


def _document(plasmid: Plasmid) -> sbol3.Document:
    document = new_document()
    insert_sequence = sbol3.Sequence(
        _local_id(plasmid.insert.iri, "seq"),
        elements=plasmid.sequence.bases.lower(),
        encoding=sbol3.IUPAC_DNA_ENCODING,
    )
    coding_start = len(plasmid.unit.gene_start.bases)
    coding_end = coding_start + len(plasmid.unit.coding_sequence.bases)
    insert = sbol3.Component(
        plasmid.insert.iri,
        [sbol3.SBO_DNA],
        roles=[SO_GENE],
        sequences=[insert_sequence],
        name=plasmid.insert_id,
        features=[
            _feature(
                insert_sequence,
                "gene start",
                0,
                coding_start,
                [SO_REGULATORY_REGION, ROLE_GENE_START],
            ),
            _feature(insert_sequence, "coding sequence", coding_start, coding_end, [SO_CDS]),
            _feature(
                insert_sequence,
                "gene end",
                coding_end,
                len(plasmid.sequence.bases),
                [SO_TERMINATOR, ROLE_GENE_END],
            ),
        ],
    )
    backbone_sequence = sbol3.Sequence(
        _local_id(plasmid.backbone.iri, "seq"),
        elements=plasmid.backbone.sequence.bases.lower(),
        encoding=sbol3.IUPAC_DNA_ENCODING,
    )
    slot = plasmid.backbone.slot
    features = []
    if slot.start < slot.end:
        features.append(_feature(backbone_sequence, "slot", slot.start, slot.end, [SO_REGION]))
    cited = sbol3.Component(
        _local_id(plasmid.backbone.iri, "backbone"),
        [sbol3.SBO_DNA],
        sequences=[backbone_sequence],
        name=plasmid.backbone.id,
        description=_backbone_facts(plasmid.backbone),
        features=features,
    )
    product = sbol3.Component(
        plasmid.product.iri,
        [sbol3.SBO_DNA],
        roles=[plasmid.backbone.iri],
        name=plasmid.plan_id,
        description=_product_facts(plasmid),
        features=[sbol3.SubComponent(insert, name=plasmid.insert_id)],
    )
    document.add(insert_sequence)
    document.add(insert)
    document.add(backbone_sequence)
    document.add(cited)
    document.add(product)
    return document


def _feature(
    sequence: sbol3.Sequence, name: str, start: int, end: int, roles: list[str]
) -> sbol3.SequenceFeature:
    location = sbol3.Range(sequence, start + 1, end, orientation=sbol3.SBOL_INLINE)
    return sbol3.SequenceFeature([location], roles=roles, orientation=sbol3.SBOL_INLINE, name=name)


def _product_facts(plasmid: Plasmid) -> str:
    return (
        f"plan={plasmid.plan_id};edit={plasmid.edit};slot={plasmid.slot};"
        f"insert={plasmid.insert_id};replication={plasmid.replication.value};"
        f"accession={plasmid.backbone.accession}"
    )


def _backbone_facts(cited: Backbone) -> str:
    slot = cited.slot
    return (
        f"accession={cited.accession};citation={cited.citation};"
        f"slot-start={slot.start};slot-end={slot.end};"
        f"upstream={slot.upstream.bases};downstream={slot.downstream.bases};"
        f"gene-start={slot.gene_start.bases};gene-end={slot.gene_end.bases}"
    )


def _unit(document: sbol3.Document, insert: sbol3.Component, path: str | Path) -> TranscriptionUnit:
    bases = _bases(document, insert, path)
    return TranscriptionUnit(
        Nucleotides(_feature_bases(insert, bases, ROLE_GENE_START, path)),
        Nucleotides(_feature_bases(insert, bases, SO_CDS, path)),
        Nucleotides(_feature_bases(insert, bases, ROLE_GENE_END, path)),
    )


def _feature_bases(component: sbol3.Component, bases: str, role: str, path: str | Path) -> str:
    found = [
        feature
        for feature in component.features
        if any(role_key(str(item)) == role_key(role) for item in feature.roles)
    ]
    if len(found) != 1:
        raise InputError(f"{path}: expected one {role_key(role)} feature, found {len(found)}.")
    location = next((item for item in found[0].locations if isinstance(item, sbol3.Range)), None)
    if location is None:
        raise InputError(f"{path}: {role_key(role)} has no range.")
    start, end = int(location.start) - 1, int(location.end)
    if not 0 <= start < end <= len(bases):
        raise InputError(f"{path}: {role_key(role)} is outside the insert.")
    return bases[start:end]


def _check_slot_feature(component: sbol3.Component, slot: Slot, path: str | Path) -> None:
    found = [feature for feature in component.features if feature.name == "slot"]
    if slot.start == slot.end:
        if found:
            raise InputError(f"{path}: an insertion slot has a range.")
        return
    if len(found) != 1:
        raise InputError(f"{path}: expected one slot feature, found {len(found)}.")
    location = next((item for item in found[0].locations if isinstance(item, sbol3.Range)), None)
    if location is None or (int(location.start) - 1, int(location.end)) != (slot.start, slot.end):
        raise InputError(f"{path}: the slot range does not match the backbone.")


def _bases(document: sbol3.Document, component: sbol3.Component, path: str | Path) -> str:
    if len(component.sequences) != 1:
        raise InputError(f"{path}: expected one sequence, found {len(component.sequences)}.")
    sequence = document.find(str(component.sequences[0]))
    if not isinstance(sequence, sbol3.Sequence) or not sequence.elements:
        raise InputError(f"{path}: the sequence is missing.")
    return str(sequence.elements).upper()


def _facts(description: str) -> dict[str, str]:
    facts: dict[str, str] = {}
    for part in description.split(";"):
        key, separator, value = part.partition("=")
        if separator:
            facts[key] = value
    return facts


def _require(facts: dict[str, str], keys: tuple[str, ...], path: str | Path) -> dict[str, str]:
    missing = [key for key in keys if key not in facts]
    if missing:
        raise InputError(f"{path}: missing {', '.join(missing)}.")
    return facts


def _local_id(iri: str, suffix: str) -> str:
    leaf = iri.rstrip("/").rsplit("/", 1)[-1].replace("-", "_")
    return f"{leaf}_{suffix}"


def _coordinate(value: str, path: str | Path) -> int:
    if not value.isdecimal():
        raise InputError(f"{path}: slot coordinate {value!r} is not an integer.")
    return int(value)
