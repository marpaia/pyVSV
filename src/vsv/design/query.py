"""Select published VSV designs by antigen class, phenotype, host, and evidence."""

from __future__ import annotations

from collections.abc import Sequence

from vsv.design.evidence import Design, EvidenceLevel, Host, Phenotype


def select(
    records: Sequence[Design],
    *,
    antigen_class: str | None = None,
    phenotype: Phenotype | None = None,
    host: Host | None = None,
    min_evidence: EvidenceLevel | None = None,
) -> list[Design]:
    """Designs that meet every filter that is set.

    A returned design keeps the observations that satisfied the phenotype and host filters.
    With neither of those filters set, it keeps every observation. Evidence is judged from the
    full card, before that filtering.
    """
    found: list[Design] = []
    for record in records:
        if min_evidence is not None and record.card.evidence_level.rank < min_evidence.rank:
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
            Design(
                record.architecture,
                observations,
                record.implementations,
                record.realization,
                record.model_source,
                record.unmodeled,
            )
        )
    return found
