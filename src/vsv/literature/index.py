"""Read and write the construction-paper index and the download manifest."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.catalog import CensusIndex, CitedPaper, ManifestRow
from vsv.literature.query import CORPUS_ROOT
from vsv.literature.reviews import REVIEWS

INDEX_NAME = "index.json"
MANIFEST_NAME = "manifest.jsonl"
CENSUS_ROOT = CORPUS_ROOT.parent / "census"


def index_path(root: Path) -> Path:
    """Path of the bibliography."""
    return root / INDEX_NAME


def manifest_path(root: Path) -> Path:
    """Path of the download manifest."""
    return root / MANIFEST_NAME


def write_index(root: Path, index: CensusIndex) -> None:
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


def read_index(root: Path) -> CensusIndex:
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
    return CensusIndex(retrieved, tuple(_paper(item, path) for item in papers))


def read_manifest(root: Path) -> tuple[ManifestRow, ...]:
    """Read the download manifest."""
    path = manifest_path(root)
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise InputError(f"{path}: {error}") from error
    rows: list[ManifestRow] = []
    for line in lines:
        if not line.strip():
            continue
        try:
            raw = json.loads(line)
        except json.JSONDecodeError as error:
            raise InputError(f"{path}: {error}") from error
        rows.append(_row(raw, path))
    return tuple(rows)


def write_manifest(root: Path, rows: Sequence[ManifestRow]) -> None:
    """Write one manifest row per line, in PMID order."""
    root.mkdir(parents=True, exist_ok=True)
    ordered = sorted(rows, key=lambda row: int(row.pmid))
    lines = [json.dumps(_row_json(row), sort_keys=True) for row in ordered]
    text = "\n".join(lines)
    if text:
        text += "\n"
    manifest_path(root).write_text(text, encoding="utf-8")


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


def _row_json(row: ManifestRow) -> dict[str, object]:
    return {
        "pmid": row.pmid,
        "pmcid": row.pmcid,
        "doi": row.doi,
        "title": row.title,
        "year": row.year,
        "source": row.source,
        "retrieved": row.retrieved,
        "sha256": row.sha256,
        "bytes": row.nbytes,
        "body_characters": row.body_characters,
        "license": row.license,
        "status": row.status,
        "detail": row.detail,
        "cited_by": list(row.cited_by),
    }


def _row(value: object, path: Path) -> ManifestRow:
    raw = _object(value, path)
    nbytes = raw.get("bytes")
    body = raw.get("body_characters")
    if not isinstance(nbytes, int) or not isinstance(body, int):
        raise InputError(f"{path}: bytes and body_characters must be integers.")
    return ManifestRow(
        _text(raw, "pmid", path),
        _text(raw, "pmcid", path),
        _text(raw, "doi", path),
        _text(raw, "title", path),
        _text(raw, "year", path),
        _text(raw, "source", path),
        _text(raw, "retrieved", path),
        _text(raw, "sha256", path),
        nbytes,
        body,
        _text(raw, "license", path),
        _text(raw, "status", path),
        _text(raw, "detail", path),
        tuple(_strings(raw.get("cited_by"), path)),
    )


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
