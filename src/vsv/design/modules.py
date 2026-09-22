"""The parts of a VSV architecture.

A module is a leader, a trailer, a viral gene, an antigen, or a reporter. An antigen class
exists only on an antigen. A transcription index exists only on a transcribed module.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from enum import StrEnum

from vsv.errors import InputError
from vsv.genome.sbol import display_id

_SLUG = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*\Z")


class ModuleKind(StrEnum):
    """What a module contributes to the vector."""

    LEADER = "leader"
    TRAILER = "trailer"
    VIRAL_GENE = "viral-gene"
    ANTIGEN = "antigen"
    REPORTER = "reporter"


class LocusStatus(StrEnum):
    """Whether a module's transcription index is known.

    This is derived from the module. A transcribed module with no locus is unknown. The leader
    and the trailer are not transcribed.
    """

    KNOWN = "known"
    UNKNOWN = "unknown"
    NOT_TRANSCRIBED = "not-transcribed"


def require_display_id(identifier: str, what: str) -> None:
    """Reject an identifier pySBOL3 would rewrite."""
    if not identifier or display_id(identifier) != identifier:
        raise InputError(f"Invalid {what} {identifier!r}.")


def _identity(module_id: str, name: str) -> None:
    require_display_id(module_id, "module id")
    if not name:
        raise InputError(f"Module {module_id} needs a name.")


def _nonnegative(module_id: str, locus: int | None) -> None:
    if locus is not None and locus < 0:
        raise InputError(f"Module {module_id} has a negative locus.")


@dataclass(frozen=True, slots=True)
class Leader:
    """The genomic leader. It is not transcribed."""

    name: str = "leader"
    id: str = "leader"
    feature_identity: str | None = None

    def __post_init__(self) -> None:
        _identity(self.id, self.name)

    @property
    def kind(self) -> ModuleKind:
        return ModuleKind.LEADER

    @property
    def locus(self) -> int | None:
        return None

    @property
    def locus_status(self) -> LocusStatus:
        return LocusStatus.NOT_TRANSCRIBED

    @property
    def transcribed(self) -> bool:
        return False

    @property
    def antigen_class(self) -> str | None:
        return None


@dataclass(frozen=True, slots=True)
class Trailer:
    """The genomic trailer. It is not transcribed."""

    name: str = "trailer"
    id: str = "trailer"
    feature_identity: str | None = None

    def __post_init__(self) -> None:
        _identity(self.id, self.name)

    @property
    def kind(self) -> ModuleKind:
        return ModuleKind.TRAILER

    @property
    def locus(self) -> int | None:
        return None

    @property
    def locus_status(self) -> LocusStatus:
        return LocusStatus.NOT_TRANSCRIBED

    @property
    def transcribed(self) -> bool:
        return False

    @property
    def antigen_class(self) -> str | None:
        return None


@dataclass(frozen=True, slots=True)
class ViralGene:
    """One of the viral genes, with its transcription index when the source states it."""

    id: str
    name: str
    locus: int | None = None
    feature_identity: str | None = None

    def __post_init__(self) -> None:
        _identity(self.id, self.name)
        _nonnegative(self.id, self.locus)

    @property
    def kind(self) -> ModuleKind:
        return ModuleKind.VIRAL_GENE

    @property
    def locus_status(self) -> LocusStatus:
        return LocusStatus.UNKNOWN if self.locus is None else LocusStatus.KNOWN

    @property
    def transcribed(self) -> bool:
        return True

    @property
    def antigen_class(self) -> str | None:
        return None


@dataclass(frozen=True, slots=True)
class Antigen:
    """An antigen module. ``antigen_class`` is a slug such as ``filovirus-glycoprotein``."""

    id: str
    name: str
    antigen_class: str
    locus: int | None = None
    feature_identity: str | None = None

    def __post_init__(self) -> None:
        _identity(self.id, self.name)
        _nonnegative(self.id, self.locus)
        if _SLUG.fullmatch(self.antigen_class) is None:
            raise InputError(f"Module {self.id} needs an antigen class slug.")

    @property
    def kind(self) -> ModuleKind:
        return ModuleKind.ANTIGEN

    @property
    def locus_status(self) -> LocusStatus:
        return LocusStatus.UNKNOWN if self.locus is None else LocusStatus.KNOWN

    @property
    def transcribed(self) -> bool:
        return True


@dataclass(frozen=True, slots=True)
class Reporter:
    """A reporter module, with its transcription index when the source states it."""

    id: str
    name: str
    locus: int | None = None
    feature_identity: str | None = None

    def __post_init__(self) -> None:
        _identity(self.id, self.name)
        _nonnegative(self.id, self.locus)

    @property
    def kind(self) -> ModuleKind:
        return ModuleKind.REPORTER

    @property
    def locus_status(self) -> LocusStatus:
        return LocusStatus.UNKNOWN if self.locus is None else LocusStatus.KNOWN

    @property
    def transcribed(self) -> bool:
        return True

    @property
    def antigen_class(self) -> str | None:
        return None


Module = Leader | Trailer | ViralGene | Antigen | Reporter
Transcribed = ViralGene | Antigen | Reporter


def at_locus(module: Module, locus: int | None) -> Module:
    """A copy of ``module`` at ``locus``, with no link back to a sequence feature.

    The leader and the trailer keep no locus. ``locus`` is ignored for them.
    """
    match module:
        case Leader() | Trailer():
            return replace(module, feature_identity=None)
        case ViralGene() | Antigen() | Reporter():
            return replace(module, locus=locus, feature_identity=None)
