"""Antigen-class relationships read from the local ontology.

Subclass means ``rdfs:subClassOf`` between two ``https://vsv.bio/antigen-class/`` terms.
The file is parsed once. Imports are not resolved and no reasoner runs.
"""

from __future__ import annotations

from functools import cache
from pathlib import Path

import rdflib
from rdflib import RDFS, URIRef

from vsv.errors import InputError

_CLASS = "https://vsv.bio/antigen-class/"
_FORM = "https://vsv.bio/glycoprotein-form/"
_EDITING_RUN = URIRef("https://vsv.bio/editing-run")


def antigen_matches(slug: str, ancestor: str) -> bool:
    """Whether ``slug`` is ``ancestor`` or a subclass of it."""
    return ancestor in _lineage(slug)


def is_leaf(slug: str) -> bool:
    """Whether no antigen class declares ``slug`` as its parent."""
    return slug not in set(_parents().values())


def class_label(slug: str) -> str:
    """The ``rdfs:label`` of an antigen class, or the slug when it has none."""
    label = _graph().value(URIRef(_CLASS + slug), RDFS.label)
    return str(label) if label is not None else slug


def editing_length(slug: str, form: str) -> int | None:
    """The adenosine-run length declared for ``form`` on ``slug`` or an ancestor.

    The class itself wins over its parent. ``None`` means this form is not inspected.
    """
    graph = _graph()
    for current in _lineage(slug):
        for node in graph.objects(URIRef(_CLASS + current), _EDITING_RUN):
            value = graph.value(node, URIRef(_FORM + form))
            if value is not None:
                return int(str(value))
    return None


@cache
def _graph() -> rdflib.Graph:
    graph = rdflib.Graph()
    graph.parse(_ontology_path())
    return graph


def _slug(term: object) -> str | None:
    text = str(term)
    if not text.startswith(_CLASS):
        return None
    return text.removeprefix(_CLASS)


def _parents() -> dict[str, str]:
    """Child slug to parent slug, for edges that stay inside antigen-class."""
    edges: dict[str, str] = {}
    for child, parent in _graph().subject_objects(RDFS.subClassOf):
        child_slug, parent_slug = _slug(child), _slug(parent)
        if child_slug is None or parent_slug is None:
            continue
        edges[child_slug] = parent_slug
    return edges


def _lineage(slug: str) -> list[str]:
    """``slug``, then its parent, and so on. A cycle stops the walk."""
    lineage: list[str] = []
    seen: set[str] = set()
    current: str | None = slug
    parents = _parents()
    while current is not None and current not in seen:
        lineage.append(current)
        seen.add(current)
        current = parents.get(current)
    return lineage


def _ontology_path() -> Path:
    root = Path(__file__).resolve().parents[3]
    checkout = root / "ontology" / "synthetic-virology.ttl"
    if (root / "pyproject.toml").is_file() and checkout.is_file():
        return checkout
    packaged = Path(__file__).resolve().parents[1] / "synthetic-virology.ttl"
    if packaged.is_file():
        return packaged
    raise InputError("The ontology is not installed.")
