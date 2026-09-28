"""Records for one cited paper and the file downloaded for it."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class CitedPaper:
    """A PMID cited by one or more index reviews."""

    pmid: str
    pmcid: str
    doi: str
    title: str
    year: str
    open_access: bool
    in_corpus: bool
    cited_by: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class CensusIndex:
    """The review list and the deduplicated papers those reviews cite."""

    retrieved: str
    papers: tuple[CitedPaper, ...]


@dataclass(frozen=True, slots=True)
class ManifestRow:
    """One download attempt. ``sha256`` is empty when no file was kept."""

    pmid: str
    pmcid: str
    doi: str
    title: str
    year: str
    source: str
    retrieved: str
    sha256: str
    nbytes: int
    body_characters: int
    license: str
    status: str
    detail: str
    cited_by: tuple[str, ...]
