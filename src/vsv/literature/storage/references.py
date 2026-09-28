"""Reference-snapshot records and their serialization."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.discovery.reviews import REVIEWS
from vsv.literature.storage.paths import INDEX_NAME


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
class ReferenceSnapshot:
    """The review list and the deduplicated papers those reviews cite."""

    retrieved: str
    papers: tuple[CitedPaper, ...]


def index_path(root: Path) -> Path:
    """Path of the bibliography."""
    return root / INDEX_NAME


def write_index(root: Path, index: ReferenceSnapshot) -> None:
    """Write the bibliography. Review order follows ``REVIEWS``."""
    root.mkdir(parents=True, exist_ok=True)
    payload = {
        "retrieved": index.retrieved,
        "reviews": [
            {"pmid": review.pmid, "pmcid": review.pmcid, "label": review.label}
            for review in REVIEWS
        ],
        "papers": [_paper_json(paper) for paper in index.papers],
    }
    index_path(root).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def read_index(root: Path) -> ReferenceSnapshot:
    """Read the bibliography written by :func:`write_index`."""
    path = index_path(root)
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise InputError(f"{path}: {error}") from error
    if not isinstance(raw, dict):
        raise InputError(f"{path}: expected an object.")
    retrieved = raw.get("retrieved")
    papers = raw.get("papers")
    if not isinstance(retrieved, str) or not isinstance(papers, list):
        raise InputError(f"{path}: expected retrieved and papers.")
    return ReferenceSnapshot(retrieved, tuple(_paper(item, path) for item in papers))


def _paper_json(paper: CitedPaper) -> dict[str, object]:
    return {
        "pmid": paper.pmid,
        "pmcid": paper.pmcid,
        "doi": paper.doi,
        "title": paper.title,
        "year": paper.year,
        "open_access": paper.open_access,
        "in_corpus": paper.in_corpus,
        "cited_by": list(paper.cited_by),
    }


def _paper(value: object, path: Path) -> CitedPaper:
    raw = _object(value, path)
    return CitedPaper(
        _text(raw, "pmid", path),
        _text(raw, "pmcid", path),
        _text(raw, "doi", path),
        _text(raw, "title", path),
        _text(raw, "year", path),
        _flag(raw, "open_access", path),
        _flag(raw, "in_corpus", path),
        tuple(_strings(raw.get("cited_by"), path)),
    )


def _object(value: object, path: Path) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise InputError(f"{path}: expected an object.")
    return value


def _text(raw: Mapping[str, object], key: str, path: Path) -> str:
    value = raw.get(key, "")
    if not isinstance(value, str):
        raise InputError(f"{path}: {key} must be a string.")
    return value


def _flag(raw: Mapping[str, object], key: str, path: Path) -> bool:
    value = raw.get(key)
    if not isinstance(value, bool):
        raise InputError(f"{path}: {key} must be a boolean.")
    return value


def _strings(value: object, path: Path) -> list[str]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise InputError(f"{path}: cited_by must be a list of strings.")
    return list(value)
