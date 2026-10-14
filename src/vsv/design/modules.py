"""Genes on a VSV vaccine plan, and the edit that places an antigen."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Literal

from vsv.design.vector import GlycoproteinForm
from vsv.errors import InputError
from vsv.genes import Gene
from vsv.genome.sequence import Nucleotides

_DISPLAY_ID = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")
_SLUG = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*\Z")


def _require_id(identifier: str, what: str) -> None:
    if _DISPLAY_ID.fullmatch(identifier) is None:
        raise InputError(f"Invalid {what} {identifier!r}.")


@dataclass(frozen=True, slots=True)
class Viral:
    """One of the five VSV genes."""

    gene: Gene

    @property
    def id(self) -> str:
        return self.gene.value

    @property
    def name(self) -> str:
        return self.gene.value


@dataclass(frozen=True, slots=True)
class Antigen:
    """A non-viral gene. ``coding`` is kept in memory and written only when it is set."""

    id: str
    name: str
    antigen_class: str
    strain: str = ""
    accession: str = ""
    form: GlycoproteinForm = GlycoproteinForm.UNSPECIFIED
    coding: str = ""

    def __post_init__(self) -> None:
        _require_id(self.id, "module id")
        if not self.name:
            raise InputError(f"Module {self.id} needs a name.")
        if _SLUG.fullmatch(self.antigen_class) is None:
            raise InputError(f"Module {self.id} needs an antigen class slug.")
        if self.coding:
            Nucleotides(self.coding)


Module = Viral | Antigen


@dataclass(frozen=True, slots=True)
class Replace:
    """The antigen ``module`` stands where the VSV glycoprotein stood."""

    module: str


@dataclass(frozen=True, slots=True)
class Insert:
    """``module`` was added after the leader or after the glycoprotein."""

    after: Literal["leader", "g"]
    module: str

    def __post_init__(self) -> None:
        if self.after not in {"leader", "g"}:
            raise InputError(f"Unknown insertion site {self.after!r}.")


Edit = Replace | Insert
