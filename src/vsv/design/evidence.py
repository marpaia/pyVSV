"""Observations and the design they belong to.

An observation is one published result: a phenotype, the host it was measured in, the assay, and
a citation. A design is an architecture plus the implementations and observations recorded for
it. The chassis card is computed from those, so it cannot disagree with the design it describes.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import TYPE_CHECKING

from vsv.design.architecture import Architecture, Replication
from vsv.errors import InputError

if TYPE_CHECKING:
    from vsv.design.rules import Property, Violation

_DOI = re.compile(r"(?:https://doi\.org/)?10\.\d{4,9}/\S+\Z", re.IGNORECASE)
_PMID = re.compile(r"(?:PMID:|https://pubmed\.ncbi\.nlm\.nih\.gov/)?(\d+)/?\Z", re.IGNORECASE)


class Host(StrEnum):
    """The biological context an observation was made in."""

    CELL_LINE = "cell-line"
    MOUSE = "mouse"
    HAMSTER = "hamster"
    NONHUMAN_PRIMATE = "nonhuman-primate"
    HUMAN = "human"
    OTHER = "other"


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


class EvidenceLevel(StrEnum):
    """The highest biological context among a design's observations."""

    NONE = "none"
    IN_VITRO = "in-vitro"
    RODENT = "rodent"
    NONHUMAN_PRIMATE = "nonhuman-primate"
    HUMAN = "human"

    @property
    def rank(self) -> int:
        """Position in the evidence order, from none up to human."""
        return list(EvidenceLevel).index(self)


_HOST_LEVEL = {
    Host.CELL_LINE: EvidenceLevel.IN_VITRO,
    Host.MOUSE: EvidenceLevel.RODENT,
    Host.HAMSTER: EvidenceLevel.RODENT,
    Host.NONHUMAN_PRIMATE: EvidenceLevel.NONHUMAN_PRIMATE,
    Host.HUMAN: EvidenceLevel.HUMAN,
}


def is_citation(citation: str) -> bool:
    """Whether ``citation`` is a DOI or a PubMed identifier."""
    if _DOI.fullmatch(citation) is not None:
        return True
    return _PMID.fullmatch(citation) is not None


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


@dataclass(frozen=True, slots=True)
class Observation:
    """One cited result for an implementation of a design."""

    id: str
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
    hosts: list[str] = []
    level = EvidenceLevel.NONE
    for observation in observations:
        if observation.host.value not in hosts:
            hosts.append(observation.host.value)
        candidate = _HOST_LEVEL.get(observation.host, EvidenceLevel.NONE)
        if candidate.rank > level.rank:
            level = candidate
    return ChassisCard(
        architecture.id,
        architecture.replication,
        architecture.antigen_classes,
        tuple(hosts),
        level,
    )


@dataclass(frozen=True, slots=True)
class Design:
    """One curated design: its architecture, implementations, and observations."""

    architecture: Architecture
    observations: tuple[Observation, ...] = ()
    implementations: tuple[Implementation, ...] = ()
    realization: str | None = None
    model_source: str | None = None
    unmodeled: tuple[Unmodeled, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "observations", tuple(self.observations))
        object.__setattr__(self, "implementations", tuple(self.implementations))
        object.__setattr__(self, "unmodeled", tuple(self.unmodeled))
        known = {item.id for item in self.implementations}
        for observation in self.observations:
            if observation.implementation_id not in known:
                message = (
                    f"{observation.id} cites unknown implementation "
                    f"{observation.implementation_id}."
                )
                raise InputError(message)

    @property
    def card(self) -> ChassisCard:
        return chassis_card(self.architecture, self.observations)

    def violations(self) -> tuple[Violation, ...]:
        """Design-rule failures for this architecture and its observations."""
        from vsv.design.rules import check

        return tuple(check(self.architecture, self.observations))

    def satisfies(self, requirement: Property) -> bool:
        """Whether this design meets ``requirement``."""
        from vsv.design.rules import satisfies

        return satisfies(self.architecture, self.observations, requirement)
