"""Structural design rules for VSV architectures.

Rules inspect an architecture and the observations attached to it. They do not rewrite the
design. :func:`satisfies` answers whether a design meets one of the named properties.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Protocol

from vsv.architecture import Architecture, ModuleKind, ParentStatus, Replication
from vsv.genes import Gene

_DOI = re.compile(r"(?:https://doi\.org/)?10\.\d{4,9}/\S+\Z", re.IGNORECASE)
_PMID = re.compile(r"(?:PMID:|https://pubmed\.ncbi\.nlm\.nih\.gov/)?(\d+)/?\Z", re.IGNORECASE)
_ORDER_RULES = frozenset({"loci", "viral-gene-order", "termini"})
_TERMINAL = frozenset({ModuleKind.LEADER, ModuleKind.TRAILER})


class Host(StrEnum):
    """The biological context an observation was made in."""

    CELL_LINE = "cell-line"
    MOUSE = "mouse"
    HAMSTER = "hamster"
    NONHUMAN_PRIMATE = "nonhuman-primate"
    HUMAN = "human"
    OTHER = "other"


class Property(StrEnum):
    """A property a design can satisfy."""

    CANONICAL_ORDER = "canonical-order"
    CITED = "cited"
    PRECLINICAL_OR_ABOVE = "preclinical-or-above"


_PRECLINICAL = frozenset({Host.MOUSE, Host.HAMSTER, Host.NONHUMAN_PRIMATE, Host.HUMAN})
_HOSTS = frozenset(host.value for host in Host)
_GENE_ORDER = {gene.value: index for index, gene in enumerate(Gene)}
_GENE_ORDER_TEXT = ", ".join(Gene)


class ObservationRecord(Protocol):
    """The citation and host :func:`check` reads from an observation."""

    @property
    def citation(self) -> str: ...

    @property
    def host(self) -> str: ...


@dataclass(frozen=True, slots=True)
class Violation:
    """One design rule an architecture or its observations do not meet."""

    rule: str
    message: str


def is_citation(citation: str) -> bool:
    """Whether ``citation`` is a DOI or a PubMed identifier."""
    if _DOI.fullmatch(citation) is not None:
        return True
    match = _PMID.fullmatch(citation)
    return match is not None


def _is_citation(citation: str) -> bool:
    return is_citation(citation)


KNOWLEDGE_GAPS = frozenset({"locus-unknown", "replication-unknown", "parent-unknown"})


def _loci(architecture: Architecture) -> list[Violation]:
    body = [module for module in architecture.modules if module.kind not in _TERMINAL]
    unknown = [module.id for module in body if module.locus is None]
    known = [module.locus for module in body if module.locus is not None]
    violations = []
    if unknown:
        names = ", ".join(unknown)
        violations.append(
            Violation("locus-unknown", f"{architecture.id}: locus unknown for {names}.")
        )
    if not unknown and known != list(range(len(body))):
        violations.append(
            Violation("loci", f"{architecture.id}: transcribed loci are not contiguous from 0.")
        )
    elif unknown and len(known) != len(set(known)):
        violations.append(Violation("loci", f"{architecture.id}: known loci are repeated."))
    terminal_locus = any(
        module.locus is not None for module in architecture.modules if module.kind in _TERMINAL
    )
    if terminal_locus:
        violations.append(
            Violation("loci", f"{architecture.id}: leader and trailer have no locus.")
        )
    return violations


def _viral_gene_order(architecture: Architecture) -> list[Violation]:
    viral = [
        module
        for module in architecture.modules
        if module.kind is ModuleKind.VIRAL_GENE
        and module.id in _GENE_ORDER
        and module.locus is not None
    ]
    positions = [(_GENE_ORDER[module.id], module.locus) for module in viral]
    if positions != sorted(positions):
        return [
            Violation(
                "viral-gene-order",
                f"{architecture.id}: viral genes are not in the order {_GENE_ORDER_TEXT}.",
            )
        ]
    return []


def _termini(architecture: Architecture) -> list[Violation]:
    leaders = [module for module in architecture.modules if module.kind is ModuleKind.LEADER]
    trailers = [module for module in architecture.modules if module.kind is ModuleKind.TRAILER]
    violations = []
    if len(leaders) > 1 or (leaders and architecture.modules[0].kind is not ModuleKind.LEADER):
        violations.append(Violation("termini", f"{architecture.id}: the leader is not first."))
    if len(trailers) > 1 or (trailers and architecture.modules[-1].kind is not ModuleKind.TRAILER):
        violations.append(Violation("termini", f"{architecture.id}: the trailer is not last."))
    return violations


def _cited(architecture: Architecture) -> list[Violation]:
    if architecture.parent_status is ParentStatus.ROOT:
        return []
    if any(_is_citation(citation) for citation in architecture.citations):
        return []
    return [
        Violation(
            "cited",
            f"{architecture.id}: a derived architecture needs a DOI or PMID citation.",
        )
    ]


def _observations(
    architecture: Architecture, observations: Sequence[ObservationRecord]
) -> list[Violation]:
    violations = []
    for index, observation in enumerate(observations):
        if not _is_citation(observation.citation):
            violations.append(
                Violation(
                    "observation",
                    f"{architecture.id}: observation {index} needs a DOI or PMID citation.",
                )
            )
        if observation.host not in _HOSTS:
            violations.append(
                Violation(
                    "observation",
                    f"{architecture.id}: observation {index} names an unknown host.",
                )
            )
    return violations


def _knowledge(architecture: Architecture) -> list[Violation]:
    gaps = []
    if architecture.replication is Replication.UNKNOWN:
        gaps.append(Violation("replication-unknown", f"{architecture.id}: replication is unknown."))
    if architecture.parent_status is ParentStatus.UNKNOWN:
        gaps.append(Violation("parent-unknown", f"{architecture.id}: parent design is unknown."))
    return gaps


def _replication(architecture: Architecture) -> list[Violation]:
    if architecture.replication is Replication.UNKNOWN:
        return []
    glycoprotein = Gene.GLYCOPROTEIN
    viral_g = any(
        module.kind is ModuleKind.VIRAL_GENE and module.id == glycoprotein
        for module in architecture.modules
    )
    has_g = any(module.id == glycoprotein for module in architecture.modules)
    derivation = architecture.derivation
    replaced_g = (
        derivation is not None
        and derivation.activity == "replace-module"
        and derivation.target == glycoprotein
        and derivation.module_id is not None
        and any(module.id == derivation.module_id for module in architecture.modules)
    )
    if architecture.replication is Replication.SINGLE_CYCLE and viral_g:
        return [
            Violation(
                "replication",
                f"{architecture.id}: a single-cycle architecture still carries {glycoprotein}.",
            )
        ]
    if architecture.replication is Replication.COMPETENT and not (has_g or replaced_g):
        return [
            Violation(
                "replication",
                f"{architecture.id}: a replication-competent architecture is missing "
                f"{glycoprotein}.",
            )
        ]
    return []


def check(
    architecture: Architecture, observations: Sequence[ObservationRecord] = ()
) -> list[Violation]:
    """Every design rule ``architecture`` and ``observations`` fail, in rule order."""
    return [
        *_loci(architecture),
        *_viral_gene_order(architecture),
        *_termini(architecture),
        *_cited(architecture),
        *_observations(architecture, observations),
        *_replication(architecture),
        *_knowledge(architecture),
    ]


def satisfies(
    architecture: Architecture,
    observations: Sequence[ObservationRecord],
    requirement: Property,
) -> bool:
    """Whether ``architecture`` meets ``requirement`` given ``observations``."""
    violations = check(architecture, observations)
    if requirement is Property.CANONICAL_ORDER:
        return not any(violation.rule in _ORDER_RULES for violation in violations)
    if requirement is Property.CITED:
        return not any(violation.rule == "cited" for violation in violations)
    return any(observation.host in _PRECLINICAL for observation in observations)
