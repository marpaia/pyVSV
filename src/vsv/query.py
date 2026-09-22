"""Select published VSV designs by antigen class, phenotype, host, and evidence."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from vsv.architecture import Architecture
from vsv.evidence import (
    ChassisCard,
    EvidenceLevel,
    Implementation,
    Observation,
    Phenotype,
    Unmodeled,
    chassis_card,
)
from vsv.rules import Host

_LEVEL_RANK = {level: index for index, level in enumerate(EvidenceLevel)}


@dataclass(frozen=True, slots=True)
class DesignRecord:
    """One curated design: its architecture, observations, and optional sequence accession."""

    architecture: Architecture
    observations: tuple[Observation, ...]
    implementations: tuple[Implementation, ...] = ()
    realization: str | None = None
    model_source: str | None = None
    unmodeled: tuple[Unmodeled, ...] = ()

    @property
    def card(self) -> ChassisCard:
        return chassis_card(self.architecture, self.observations)


def select(
    records: Sequence[DesignRecord],
    *,
    antigen_class: str | None = None,
    phenotype: Phenotype | None = None,
    host: Host | None = None,
    min_evidence: EvidenceLevel | None = None,
) -> list[DesignRecord]:
    """Designs that meet every filter that is set.

    A returned record keeps the observations that satisfied the phenotype and host filters.
    With neither of those filters set, it keeps every observation. Evidence is judged from the
    full card, before that filtering.
    """
    found: list[DesignRecord] = []
    threshold = None if min_evidence is None else _LEVEL_RANK[min_evidence]
    for record in records:
        if threshold is not None and _LEVEL_RANK[record.card.evidence_level] < threshold:
            continue
        if antigen_class is not None and antigen_class not in record.card.antigen_classes:
            continue
        observations = tuple(
            observation
            for observation in record.observations
            if (phenotype is None or observation.phenotype is phenotype)
            and (host is None or observation.host is host)
        )
        if (phenotype is not None or host is not None) and not observations:
            continue
        if phenotype is None and host is None:
            observations = record.observations
        found.append(
            DesignRecord(
                record.architecture,
                observations,
                record.implementations,
                record.realization,
                record.model_source,
                record.unmodeled,
            )
        )
    return found
