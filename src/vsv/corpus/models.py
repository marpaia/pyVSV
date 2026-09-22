"""Records stored for one VSV paper and the design facts extracted from it."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Paper:
    """One PubMed record. Abstract text has long nucleotide runs removed."""

    pmid: str
    doi: str
    pmcid: str
    year: str
    pub_types: tuple[str, ...]
    mesh: tuple[str, ...]
    title: str
    abstract: str
    open_access: bool
    source: str
    retrieved: str
    in_europepmc: bool


@dataclass(frozen=True, slots=True)
class Ledger:
    """The PubMed result set that ``papers.jsonl`` is expected to hold."""

    query: str
    residual_query: str
    retrieved: str
    pubmed_count: int
    residual_count: int
    pmids: tuple[str, ...]
    europepmc_misses: tuple[str, ...]
    digest: str


@dataclass(frozen=True, slots=True)
class Quoted:
    """A value taken from one sentence, or ``unknown`` when the text does not state it."""

    value: str
    quote: str

    @property
    def known(self) -> bool:
        """Whether the text stated a value."""
        return self.value != "unknown"


UNKNOWN = Quoted("unknown", "")


@dataclass(frozen=True, slots=True)
class Claim:
    """One result sentence: phenotype, host, assay, and outcome."""

    phenotype: str
    host: str
    assay: str
    outcome: str
    quote: str


@dataclass(frozen=True, slots=True)
class Extraction:
    """Classification and architecture fields for one paper.

    A field is ``unknown`` when no sentence states it. ``conflicts`` names fields whose
    sentences disagree. Those fields are ``unknown`` too.
    """

    pmid: str
    kind: str
    review: bool
    status: str
    backbone: Quoted
    replication: Quoted
    edit_activity: Quoted
    edit_target: Quoted
    antigen_name: Quoted
    antigen_slug: Quoted
    locus: Quoted
    stage: Quoted
    claims: tuple[Claim, ...]
    conflicts: tuple[str, ...]
    text_source: str
