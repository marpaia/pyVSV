"""Classification, extraction, and cluster records."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class PaperKind(StrEnum):
    """What a corpus paper is about."""

    VACCINE_VECTOR = "vaccine-vector"
    ONCOLYTIC = "oncolytic"
    REPORTER = "reporter"
    BASIC_VIROLOGY = "basic-virology"
    STRUCTURE = "structure"
    EPIDEMIOLOGY = "epidemiology"
    PSEUDOTYPE = "pseudotype-reagent"
    REVIEW = "review"
    OTHER = "other"


VECTOR_KINDS = frozenset({PaperKind.VACCINE_VECTOR, PaperKind.ONCOLYTIC, PaperKind.REPORTER})


@dataclass(frozen=True, slots=True)
class Classification:
    """The class of one paper. ``review`` is also set when the publication type says so."""

    kind: PaperKind
    review: bool
    status: str


@dataclass(frozen=True, slots=True)
class QuotedValue:
    """A value taken from one sentence, or ``unknown`` when the text does not state it."""

    value: str
    quote: str

    @property
    def known(self) -> bool:
        """Whether the text stated a value."""
        return self.value != "unknown"


UNKNOWN = QuotedValue("unknown", "")


@dataclass(frozen=True, slots=True)
class ExtractedClaim:
    """One result sentence: phenotype, host, assay, and outcome."""

    phenotype: str
    host: str
    assay: str
    outcome: str
    quote: str


@dataclass(frozen=True, slots=True)
class PaperExtraction:
    """Classification and architecture fields for one paper.

    A field is ``unknown`` when no sentence states it. ``conflicts`` names fields whose
    sentences disagree. Those fields are ``unknown`` too.
    """

    pmid: str
    kind: str
    review: bool
    status: str
    backbone: QuotedValue
    replication: QuotedValue
    edit_activity: QuotedValue
    edit_target: QuotedValue
    antigen_name: QuotedValue
    antigen_slug: QuotedValue
    locus: QuotedValue
    stage: QuotedValue
    claims: tuple[ExtractedClaim, ...]
    conflicts: tuple[str, ...]
    text_source: str


@dataclass(frozen=True, slots=True)
class PaperCluster:
    """Papers that share one architecture key."""

    id: str
    name: str
    backbone: str
    replication: str
    edit_activity: str
    edit_target: str
    antigen_name: str
    antigen_slug: str
    locus: str
    pmids: tuple[str, ...]
    hosts: tuple[str, ...]
    phenotypes: tuple[str, ...]
    evidence: str
    gaps: tuple[str, ...]
    match: str | None
