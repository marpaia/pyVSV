"""Manifest records, their serialization, and the documents found on disk."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.records import AcquisitionAttempt, Document, Paper
from vsv.literature.storage.paths import MANIFEST_NAME, acquired_dir, census_dir, fulltext_dir


@dataclass(frozen=True, slots=True)
class ManifestRow:
    """One download attempt. ``sha256`` is empty when no file was kept."""

    pmid: str
    pmcid: str
    doi: str
    title: str
    year: str
    source: str
    retrieved: str
    sha256: str
    nbytes: int
    body_characters: int
    license: str
    status: str
    detail: str
    cited_by: tuple[str, ...]


def manifest_path(root: Path) -> Path:
    """Path of the download manifest."""
    return root / MANIFEST_NAME


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


def document_from_row(row: ManifestRow) -> Document:
    """Split a manifest status string into format, method, and availability."""
    kind, availability, method = _decode(row.status, row.sha256, row.source)
    location = _location(row, kind, availability)
    attempt = AcquisitionAttempt(row.retrieved, row.source, method, availability, row.detail)
    return Document(
        row.pmid,
        kind,
        location,
        row.sha256,
        row.nbytes,
        row.source,
        row.retrieved,
        availability,
        row.license,
        row.body_characters,
        method,
        row.detail,
        (attempt,),
    )


def load_documents(root: Path, papers: Mapping[str, Paper]) -> tuple[Document, ...]:
    """Manifest rows plus article files saved beside the manifest."""
    found: list[Document] = []
    covered: set[str] = set()
    manifest = manifest_path(census_dir(root))
    if manifest.is_file():
        for row in read_manifest(census_dir(root)):
            document = document_from_row(row)
            found.append(document)
            if document.location:
                covered.add(document.location)
    acquired = acquired_dir(root)
    if acquired.is_dir():
        for kind, suffix in (("html", ".html"), ("pdf", ".pdf")):
            for path in sorted(acquired.glob(f"*{suffix}")):
                location = f"census/acquired/{path.name}"
                if location in covered or not path.stem.isdigit():
                    continue
                found.append(_file_document(path, path.stem, kind, location, "browser"))
                covered.add(location)
    fulltext = fulltext_dir(root)
    if fulltext.is_dir():
        by_pmcid = {paper.pmcid: paper.pmid for paper in papers.values() if paper.pmcid}
        for path in sorted(fulltext.glob("*.xml")):
            location = f"corpus/fulltext/{path.name}"
            if location in covered:
                continue
            pmid = by_pmcid.get(path.stem, "")
            found.append(_file_document(path, pmid, "xml", location, "download"))
            covered.add(location)
    found.sort(key=lambda item: (int(item.pmid or 0), item.location))
    return tuple(found)


def _file_document(path: Path, pmid: str, kind: str, location: str, method: str) -> Document:
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    attempt = AcquisitionAttempt("", "", method, "available", "")
    return Document(
        pmid,
        kind,
        location,
        digest,
        path.stat().st_size,
        "",
        "",
        "available",
        "",
        0,
        method,
        "",
        (attempt,),
    )


def _decode(status: str, sha256: str, source: str) -> tuple[str, str, str]:
    if status == "xml-body":
        method = "cache" if source == "cache" else "download"
        return "xml", "available", method
    if status == "browser-html":
        return "html", "available", "browser"
    if status == "browser-pdf":
        return "pdf", "available", "browser"
    if status == "front-matter-only":
        return "xml", "partial", "download"
    if status == "needs-login":
        return "", "needs-login", "browser"
    if status == "unavailable":
        kind = "html" if sha256 else ""
        return kind, "unavailable", "browser"
    if status == "error":
        return "", "error", ""
    if status == "no-pmc":
        return "", "unavailable", ""
    return "", "unavailable", ""


def _location(row: ManifestRow, kind: str, availability: str) -> str:
    if kind == "xml" and availability == "available" and row.pmcid:
        return f"corpus/fulltext/{row.pmcid}.xml"
    if kind == "html" and row.pmid:
        return f"census/acquired/{row.pmid}.html"
    if kind == "pdf" and row.pmid:
        return f"census/acquired/{row.pmid}.pdf"
    return ""


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


def _object(value: object, path: Path) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise InputError(f"{path}: expected an object.")
    return value


def _text(raw: Mapping[str, object], key: str, path: Path) -> str:
    value = raw.get(key, "")
    if not isinstance(value, str):
        raise InputError(f"{path}: {key} must be a string.")
    return value


def _strings(value: object, path: Path) -> list[str]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise InputError(f"{path}: cited_by must be a list of strings.")
    return list(value)
