"""Select published VSV designs by antigen class, phenotype, host, and evidence.

Antigen-class subclasses are the ``rdfs:subClassOf`` edges in ``ontology/synthetic-virology.ttl``.
"""

from __future__ import annotations

from collections.abc import Sequence
from functools import cache
from pathlib import Path

import rdflib
from rdflib import RDFS

from vsv.design.evidence import Design, EvidenceLevel, Host, Phenotype
from vsv.errors import InputError

_ANTIGEN = "https://vsv.bio/antigen-class/"


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
                record.aliases,
            )
        )
    return found


def antigen_matches(slug: str, ancestor: str) -> bool:
    """Whether ``slug`` is ``ancestor`` or a subclass of it in the ontology."""
    seen: set[str] = set()
    current: str | None = slug
    parents = _antigen_parents()
    while current is not None and current not in seen:
        if current == ancestor:
            return True
        seen.add(current)
        current = parents.get(current)
    return False


def _ontology_path() -> Path:
    root = Path(__file__).resolve().parents[3]
    checkout = root / "ontology" / "synthetic-virology.ttl"
    if (root / "pyproject.toml").is_file() and checkout.is_file():
        return checkout
    packaged = Path(__file__).resolve().parents[1] / "synthetic-virology.ttl"
    if packaged.is_file():
        return packaged
    raise InputError("The synthetic virology ontology is not installed.")


@cache
def _antigen_parents() -> dict[str, str]:
    graph = rdflib.Graph()
    graph.parse(_ontology_path())
    parents: dict[str, str] = {}
    for child, parent in graph.subject_objects(RDFS.subClassOf):
        child_iri, parent_iri = str(child), str(parent)
        if child_iri.startswith(_ANTIGEN) and parent_iri.startswith(_ANTIGEN):
            parents[child_iri.removeprefix(_ANTIGEN)] = parent_iri.removeprefix(_ANTIGEN)
    return parents
