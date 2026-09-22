"""Observations and the chassis card for a VSV design.

An observation is one published result: a phenotype, the host it was measured in, the assay, and
a citation. The card is computed from the architecture and those observations, so it cannot
disagree with the design it describes.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum

from vsv.architecture import Architecture, ModuleKind, Replication
from vsv.rules import Host, is_citation
from vsv.sequence import InputError


class Phenotype(StrEnum):
    """A kind of result an experiment can report."""

    EXPRESSION = "expression"
    GROWTH = "growth"
    IMMUNOGENICITY = "immunogenicity"
    SAFETY = "safety"
    PROTECTION = "protection"


class Outcome(StrEnum):
    """Whether the paper reports the phenotype, does not, or reports a harm."""

    REPORTED = "reported"
    NOT_REPORTED = "not-reported"
    ADVERSE = "adverse"


class Assay(StrEnum):
    """The measurement named by an observation."""

    LETHAL_CHALLENGE = "lethal-challenge"
    NEUTRALIZATION = "neutralization"
    IMMUNOGENICITY = "immunogenicity"
    SAFETY = "safety"
    REACTOGENICITY = "reactogenicity"


class Stage(StrEnum):
    """Where in the development lineage an implementation sits."""

    PRECLINICAL = "preclinical"
    CLINICAL = "clinical"


@dataclass(frozen=True, slots=True)
class Unmodeled:
    """A field the current schema does not name, kept with the sentence it came from."""

    key: str
    value: str
    source: str = ""


@dataclass(frozen=True, slots=True)
class Implementation:
    """A preclinical or clinical realization of an architecture."""

    id: str
    stage: Stage
    citation: str
    note: str = ""

    def __post_init__(self) -> None:
        if not self.id:
            raise InputError("An implementation needs an id.")
        if not is_citation(self.citation):
            raise InputError(f"{self.id}: citation must be a DOI or PubMed identifier.")


class EvidenceLevel(StrEnum):
    """The highest biological context among a design's observations."""

    NONE = "none"
    IN_VITRO = "in-vitro"
    RODENT = "rodent"
    NONHUMAN_PRIMATE = "nonhuman-primate"
    HUMAN = "human"


_HOST_LEVEL = {
    Host.CELL_LINE: EvidenceLevel.IN_VITRO,
    Host.MOUSE: EvidenceLevel.RODENT,
    Host.HAMSTER: EvidenceLevel.RODENT,
    Host.NONHUMAN_PRIMATE: EvidenceLevel.NONHUMAN_PRIMATE,
    Host.HUMAN: EvidenceLevel.HUMAN,
}
_LEVEL_RANK = {level: index for index, level in enumerate(EvidenceLevel)}


@dataclass(frozen=True, slots=True)
class Observation:
    """One cited result for an architecture."""

    id: str
    architecture_id: str
    implementation_id: str
    phenotype: Phenotype
    host: Host
    assay: Assay
    outcome: Outcome
    citation: str
    note: str = ""

    def __post_init__(self) -> None:
        if not self.id or not self.implementation_id:
            raise InputError("An observation needs an id and an implementation.")
        if not is_citation(self.citation):
            raise InputError(f"{self.id}: citation must be a DOI or PubMed identifier.")


@dataclass(frozen=True, slots=True)
class ChassisCard:
    """The datasheet computed for one architecture."""

    architecture_id: str
    replication: Replication
    antigen_classes: tuple[str, ...]
    hosts: tuple[str, ...]
    evidence_level: EvidenceLevel


def chassis_card(architecture: Architecture, observations: Sequence[Observation]) -> ChassisCard:
    """Build the card for ``architecture`` from ``observations``."""
    classes: list[str] = []
    for module in architecture.modules:
        if module.antigen_class and module.antigen_class not in classes:
            classes.append(module.antigen_class)
    hosts: list[str] = []
    level = EvidenceLevel.NONE
    for observation in observations:
        if observation.architecture_id != architecture.id:
            raise InputError(
                f"{observation.id} belongs to {observation.architecture_id}, not {architecture.id}."
            )
        if observation.host.value not in hosts:
            hosts.append(observation.host.value)
        candidate = _HOST_LEVEL.get(observation.host, EvidenceLevel.NONE)
        if _LEVEL_RANK[candidate] > _LEVEL_RANK[level]:
            level = candidate
    return ChassisCard(
        architecture.id,
        architecture.replication,
        tuple(classes),
        tuple(hosts),
        level,
    )


def antigen_classes(architecture: Architecture) -> tuple[str, ...]:
    """Antigen classes on ``architecture``, in module order, without duplicates."""
    return chassis_card(architecture, ()).antigen_classes


def has_antigen(architecture: Architecture, antigen_class: str) -> bool:
    """Whether any module carries ``antigen_class``."""
    return any(
        module.kind is ModuleKind.ANTIGEN and module.antigen_class == antigen_class
        for module in architecture.modules
    )
