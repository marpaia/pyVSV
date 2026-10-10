"""Published results recorded for a vaccine.

An observation is one result: a phenotype, the host it was measured in, the assay, and a
citation. An implementation is a preclinical or clinical preparation of that vaccine.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import StrEnum

from vsv.errors import InputError

_DOI = re.compile(r"(?:https://doi\.org/)?10\.\d{4,9}/\S+\Z", re.IGNORECASE)
_PMID = re.compile(r"(?:PMID:|https://pubmed\.ncbi\.nlm\.nih\.gov/)?(\d+)/?\Z", re.IGNORECASE)


class Host(StrEnum):
    """The biological context an observation was made in.

    ``cynomolgus`` and ``rhesus`` are macaque species. A query for ``nonhuman-primate``
    matches either species and also an observation that names no species.
    """

    CELL_LINE = "cell-line"
    MOUSE = "mouse"
    HAMSTER = "hamster"
    NONHUMAN_PRIMATE = "nonhuman-primate"
    HUMAN = "human"
    OTHER = "other"
    GUINEA_PIG = "guinea-pig"
    FERRET = "ferret"
    CYNOMOLGUS = "cynomolgus"
    RHESUS = "rhesus"


class Timing(StrEnum):
    """Whether the vaccine was given before challenge or after exposure."""

    PROPHYLACTIC = "prophylactic"
    POST_EXPOSURE = "post-exposure"
    UNKNOWN = "unknown"


class DoseUnit(StrEnum):
    """The unit a paper gives for a vaccine dose or a challenge dose."""

    PFU = "pfu"
    LD50 = "ld50"
    TCID50 = "tcid50"
    UNSPECIFIED = "unspecified"


class Route(StrEnum):
    """The route of immunization, when the paper states one."""

    UNKNOWN = "unknown"
    INTRAMUSCULAR = "intramuscular"
    INTRAVENOUS = "intravenous"
    INTRANASAL = "intranasal"
    ORAL = "oral"
    INTRAPERITONEAL = "intraperitoneal"
    SUBCUTANEOUS = "subcutaneous"
    INTRADERMAL = "intradermal"
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
    FIELD_EFFICACY = "field-efficacy"


class Stage(StrEnum):
    """Where in the development lineage an implementation sits."""

    PRECLINICAL = "preclinical"
    CLINICAL = "clinical"


class EvidenceLevel(StrEnum):
    """The highest biological context among a vaccine's observations.

    ``rodent`` is the small-animal tier: mouse, hamster, guinea pig, and ferret.
    """

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
    Host.GUINEA_PIG: EvidenceLevel.RODENT,
    Host.FERRET: EvidenceLevel.RODENT,
    Host.NONHUMAN_PRIMATE: EvidenceLevel.NONHUMAN_PRIMATE,
    Host.CYNOMOLGUS: EvidenceLevel.NONHUMAN_PRIMATE,
    Host.RHESUS: EvidenceLevel.NONHUMAN_PRIMATE,
    Host.HUMAN: EvidenceLevel.HUMAN,
}
_NONHUMAN_PRIMATE_HOSTS = frozenset({Host.NONHUMAN_PRIMATE, Host.CYNOMOLGUS, Host.RHESUS})
_HOSTS = frozenset(host.value for host in Host)


def host_matches(observed: Host, query: Host) -> bool:
    """Whether an observation host satisfies a host query.

    A nonhuman-primate query matches cynomolgus, rhesus, and an observation that
    names no macaque species. A query for one species does not match the other.
    """
    if query is Host.NONHUMAN_PRIMATE:
        return observed in _NONHUMAN_PRIMATE_HOSTS
    return observed is query


def is_citation(citation: str) -> bool:
    """Whether ``citation`` is a DOI or a PubMed identifier."""
    if _DOI.fullmatch(citation) is not None:
        return True
    return _PMID.fullmatch(citation) is not None


def evidence_level(hosts: tuple[Host, ...] | list[Host]) -> EvidenceLevel:
    """The highest evidence tier among ``hosts``."""
    level = EvidenceLevel.NONE
    for host in hosts:
        candidate = _HOST_LEVEL.get(host, EvidenceLevel.NONE)
        if candidate.rank > level.rank:
            level = candidate
    return level


@dataclass(frozen=True, slots=True)
class Violation:
    """One design rule a vaccine does not meet."""

    rule: str
    message: str


@dataclass(frozen=True, slots=True)
class CheckResult:
    """Rule check of one vaccine. True when every rule passes."""

    violations: tuple[Violation, ...]

    def __bool__(self) -> bool:
        return not self.violations


@dataclass(frozen=True, slots=True)
class Implementation:
    """A preclinical or clinical realization of a vaccine."""

    id: str
    stage: Stage
    citation: str
    note: str = ""
    plasmid: str = ""
    cloning_sites: str = ""

    def __post_init__(self) -> None:
        if not self.id:
            raise InputError("An implementation needs an id.")
        if not is_citation(self.citation):
            raise InputError(f"{self.id}: citation must be a DOI or PubMed identifier.")


@dataclass(frozen=True, slots=True)
class Challenge:
    """A challenge named by an observation."""

    virus: str
    strain: str = ""
    route: Route = Route.UNKNOWN
    amount: str = ""
    unit: DoseUnit = DoseUnit.UNSPECIFIED

    def __post_init__(self) -> None:
        stated = (
            bool(self.virus)
            or bool(self.strain)
            or bool(self.amount)
            or self.route is not Route.UNKNOWN
            or self.unit is not DoseUnit.UNSPECIFIED
        )
        if not stated:
            raise InputError("A challenge needs a virus, a strain, a route, or a dose.")


@dataclass(frozen=True, slots=True)
class Observation:
    """One cited result for an implementation of a vaccine.

    ``day`` is days relative to challenge. A negative day is before challenge.
    """

    id: str
    implementation_id: str
    phenotype: Phenotype
    host: Host
    assay: Assay
    outcome: Outcome
    citation: str
    note: str = ""
    route: Route = Route.UNKNOWN
    timing: Timing = Timing.UNKNOWN
    day: int | None = None
    vaccine_amount: str = ""
    vaccine_unit: DoseUnit = DoseUnit.UNSPECIFIED
    challenge: Challenge | None = None

    def __post_init__(self) -> None:
        if not self.id or not self.implementation_id:
            raise InputError("An observation needs an id and an implementation.")
        if not is_citation(self.citation):
            raise InputError(f"{self.id}: citation must be a DOI or PubMed identifier.")
        if self.host.value not in _HOSTS:
            raise InputError(f"{self.id}: unknown host {self.host!r}.")
        if self.day is not None and self.timing is Timing.UNKNOWN:
            raise InputError(f"{self.id}: a day needs a timing.")
        if self.vaccine_unit is not DoseUnit.UNSPECIFIED and not self.vaccine_amount:
            raise InputError(f"{self.id}: a vaccine unit needs an amount.")
        challenge = self.challenge
        if (
            challenge is not None
            and challenge.unit is not DoseUnit.UNSPECIFIED
            and not challenge.amount
        ):
            message = f"{self.id}: a challenge unit of {challenge.unit.value} needs an amount."
            raise InputError(message)
