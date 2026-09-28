"""Shared bibliographic, collection, and document records.

One paper has one bibliographic identity. Collection membership, an acquired
file, and the text parsed from that file are separate records.
"""

from __future__ import annotations

from dataclasses import dataclass

PUBMED = "pubmed"
REVIEW_REFERENCES = "review-references"


@dataclass(frozen=True, slots=True)
class Paper:
    """Bibliographic identity of one article, independent of any collection."""

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
class Membership:
    """One paper admitted to a collection on the date of that snapshot."""

    pmid: str
    as_of: str


@dataclass(frozen=True, slots=True)
class Collection:
    """A named selection of papers, the rule that admits them, and its provenance."""

    name: str
    rule: str
    provenance: str
    members: tuple[Membership, ...]

    @property
    def pmids(self) -> tuple[str, ...]:
        """PMIDs in snapshot order."""
        return tuple(member.pmid for member in self.members)

    @property
    def retrieved(self) -> str:
        """Date of the membership snapshot. Empty when the collection has no members."""
        return self.members[0].as_of if self.members else ""


@dataclass(frozen=True, slots=True)
class Citation:
    """One paper citing another. Membership in a collection is a separate fact."""

    source: str
    target: str


@dataclass(frozen=True, slots=True)
class AcquisitionAttempt:
    """One attempt to obtain a file."""

    retrieved: str
    source: str
    method: str
    outcome: str
    detail: str


@dataclass(frozen=True, slots=True)
class Document:
    """An acquired file, or the record that no file was kept."""

    pmid: str
    format: str
    location: str
    sha256: str
    nbytes: int
    source: str
    retrieved: str
    availability: str
    license: str
    body_characters: int
    method: str
    detail: str
    history: tuple[AcquisitionAttempt, ...]


@dataclass(frozen=True, slots=True)
class TextBlock:
    """One quotable span of an article."""

    id: str
    section: str
    kind: str
    text: str


@dataclass(frozen=True, slots=True)
class DocumentText:
    """Parsed blocks tied to one document version.

    ``document`` is the SHA-256 of the file that was parsed, or empty when the
    stored text does not name a version.
    """

    pmid: str
    document: str
    blocks: tuple[TextBlock, ...]
