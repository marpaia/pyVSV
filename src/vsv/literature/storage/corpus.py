"""Paper records, the PubMed snapshot, extraction records, and audit labels."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.analysis.records import ExtractedClaim, PaperExtraction, QuotedValue
from vsv.literature.records import Paper
from vsv.literature.storage.paths import (
    AUDIT_NAME,
    CORPUS_ROOT,
    EXTRACTIONS_NAME,
    LEDGER_NAME,
    PAPERS_NAME,
)


@dataclass(frozen=True, slots=True)
class PubMedSnapshot:
    """The PubMed result set that ``papers.jsonl`` is expected to hold."""

    query: str
    residual_query: str
    retrieved: str
    pubmed_count: int
    residual_count: int
    pmids: tuple[str, ...]
    europepmc_misses: tuple[str, ...]
    digest: str


def digest_pmids(pmids: Sequence[str]) -> str:
    """SHA-256 of the sorted PMID list, one identifier per line."""
    body = "\n".join(pmids).encode()
    return hashlib.sha256(body).hexdigest()


def papers_path(root: Path | None = None) -> Path:
    """Path of the paper JSONL file."""
    return _root(root) / PAPERS_NAME


def extractions_path(root: Path | None = None) -> Path:
    """Path of the extraction JSONL file."""
    return _root(root) / EXTRACTIONS_NAME


def ledger_path(root: Path | None = None) -> Path:
    """Path of the coverage ledger."""
    return _root(root) / LEDGER_NAME


def audit_path(root: Path | None = None) -> Path:
    """Path of the hand-labeled audit sample."""
    return _root(root) / AUDIT_NAME


def write_papers(papers: Sequence[Paper], root: Path | None = None) -> None:
    """Replace the paper JSONL with ``papers``, in PMID order."""
    ordered = sorted(papers, key=lambda paper: int(paper.pmid))
    _write_jsonl(papers_path(root), [_paper_json(paper) for paper in ordered])


def load_papers(root: Path | None = None) -> list[Paper]:
    """Load the paper JSONL."""
    path = papers_path(root)
    if not path.is_file():
        raise InputError(f"{path}: corpus papers not found.")
    return [_paper(row, path) for row in _read_jsonl(path)]


def write_extractions(extractions: Sequence[PaperExtraction], root: Path | None = None) -> None:
    """Replace the extraction JSONL with ``extractions``, in PMID order."""
    ordered = sorted(extractions, key=lambda item: int(item.pmid))
    _write_jsonl(extractions_path(root), [_extraction_json(item) for item in ordered])


def load_extractions(root: Path | None = None) -> list[PaperExtraction]:
    """Load the extraction JSONL."""
    path = extractions_path(root)
    if not path.is_file():
        raise InputError(f"{path}: extractions not found.")
    return [_extraction(row, path) for row in _read_jsonl(path)]


def write_ledger(ledger: PubMedSnapshot, root: Path | None = None) -> None:
    """Replace the coverage ledger."""
    path = ledger_path(root)
    payload = {
        "query": ledger.query,
        "residual_query": ledger.residual_query,
        "retrieved": ledger.retrieved,
        "pubmed_count": ledger.pubmed_count,
        "residual_count": ledger.residual_count,
        "digest": ledger.digest,
        "europepmc_misses": list(ledger.europepmc_misses),
        "pmids": list(ledger.pmids),
    }
    _replace(path, json.dumps(payload, indent=2) + "\n")


def load_ledger(root: Path | None = None) -> PubMedSnapshot:
    """Load the coverage ledger."""
    path = ledger_path(root)
    if not path.is_file():
        raise InputError(f"{path}: corpus ledger not found.")
    raw = _object(_load(path), path)
    pmids = _strings(raw.get("pmids"), path, "pmids")
    misses = _strings(raw.get("europepmc_misses"), path, "europepmc_misses")
    return PubMedSnapshot(
        _text(raw, "query", path),
        _text(raw, "residual_query", path),
        _text(raw, "retrieved", path),
        _count(raw, "pubmed_count", path),
        _count(raw, "residual_count", path),
        tuple(pmids),
        tuple(misses),
        _text(raw, "digest", path),
    )


def load_audit(root: Path | None = None) -> list[dict[str, object]] | None:
    """The audit labels, or ``None`` when the sample has not been written."""
    path = audit_path(root)
    if not path.is_file():
        return None
    raw = _load(path)
    if not isinstance(raw, dict):
        raise InputError(f"{path}: audit file must be an object.")
    labels = raw.get("labels")
    if not isinstance(labels, list):
        raise InputError(f"{path}: audit labels must be a list.")
    parsed: list[dict[str, object]] = []
    for item in labels:
        if not isinstance(item, dict):
            raise InputError(f"{path}: each audit label must be an object.")
        parsed.append(item)
    return parsed


def _root(root: Path | None) -> Path:
    return CORPUS_ROOT if root is None else root


def _paper_json(paper: Paper) -> dict[str, object]:
    return {
        "pmid": paper.pmid,
        "doi": paper.doi,
        "pmcid": paper.pmcid,
        "year": paper.year,
        "pub_types": list(paper.pub_types),
        "mesh": list(paper.mesh),
        "title": paper.title,
        "abstract": paper.abstract,
        "open_access": paper.open_access,
        "source": paper.source,
        "retrieved": paper.retrieved,
        "in_europepmc": paper.in_europepmc,
    }


def _quoted_json(quoted: QuotedValue) -> dict[str, str]:
    return {"value": quoted.value, "quote": quoted.quote}


def _extraction_json(item: PaperExtraction) -> dict[str, object]:
    return {
        "pmid": item.pmid,
        "kind": item.kind,
        "review": item.review,
        "status": item.status,
        "backbone": _quoted_json(item.backbone),
        "replication": _quoted_json(item.replication),
        "edit_activity": _quoted_json(item.edit_activity),
        "edit_target": _quoted_json(item.edit_target),
        "antigen_name": _quoted_json(item.antigen_name),
        "antigen_slug": _quoted_json(item.antigen_slug),
        "locus": _quoted_json(item.locus),
        "stage": _quoted_json(item.stage),
        "claims": [
            {
                "phenotype": claim.phenotype,
                "host": claim.host,
                "assay": claim.assay,
                "outcome": claim.outcome,
                "quote": claim.quote,
            }
            for claim in item.claims
        ],
        "conflicts": list(item.conflicts),
        "text_source": item.text_source,
    }


def _paper(raw: Mapping[str, object], path: Path) -> Paper:
    return Paper(
        _text(raw, "pmid", path),
        _text(raw, "doi", path),
        _text(raw, "pmcid", path),
        _text(raw, "year", path),
        tuple(_strings(raw.get("pub_types"), path, "pub_types")),
        tuple(_strings(raw.get("mesh"), path, "mesh")),
        _text(raw, "title", path),
        _text(raw, "abstract", path),
        _flag(raw, "open_access", path),
        _text(raw, "source", path),
        _text(raw, "retrieved", path),
        _flag(raw, "in_europepmc", path),
    )


def _quoted(value: object, path: Path, field: str) -> QuotedValue:
    raw = _object(value, path, field)
    return QuotedValue(_text(raw, "value", path), _text(raw, "quote", path))


def _extraction(raw: Mapping[str, object], path: Path) -> PaperExtraction:
    claims_raw = raw.get("claims")
    if not isinstance(claims_raw, list):
        raise InputError(f"{path}: claims must be a list.")
    claims = []
    for item in claims_raw:
        claim = _object(item, path, "claim")
        claims.append(
            ExtractedClaim(
                _text(claim, "phenotype", path),
                _text(claim, "host", path),
                _text(claim, "assay", path),
                _text(claim, "outcome", path),
                _text(claim, "quote", path),
            )
        )
    return PaperExtraction(
        _text(raw, "pmid", path),
        _text(raw, "kind", path),
        _flag(raw, "review", path),
        _text(raw, "status", path),
        _quoted(raw.get("backbone"), path, "backbone"),
        _quoted(raw.get("replication"), path, "replication"),
        _quoted(raw.get("edit_activity"), path, "edit_activity"),
        _quoted(raw.get("edit_target"), path, "edit_target"),
        _quoted(raw.get("antigen_name"), path, "antigen_name"),
        _quoted(raw.get("antigen_slug"), path, "antigen_slug"),
        _quoted(raw.get("locus"), path, "locus"),
        _quoted(raw.get("stage"), path, "stage"),
        tuple(claims),
        tuple(_strings(raw.get("conflicts"), path, "conflicts")),
        _text(raw, "text_source", path),
    )


def _read_jsonl(path: Path) -> list[Mapping[str, object]]:
    rows = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        try:
            parsed = json.loads(line)
        except json.JSONDecodeError as error:
            raise InputError(f"{path}:{number}: {error.msg}.") from error
        rows.append(_object(parsed, path, f"line {number}"))
    return rows


def _write_jsonl(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    text = "".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows)
    _replace(path, text)


def _replace(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def _load(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise InputError(f"{path}: {error.msg}.") from error


def _object(value: object, path: Path, field: str = "record") -> Mapping[str, object]:
    if isinstance(value, dict):
        return value
    raise InputError(f"{path}: {field} must be an object.")


def _text(raw: Mapping[str, object], key: str, path: Path) -> str:
    value = raw.get(key)
    if isinstance(value, str):
        return value
    raise InputError(f"{path}: {key} must be a string.")


def _flag(raw: Mapping[str, object], key: str, path: Path) -> bool:
    value = raw.get(key)
    if isinstance(value, bool):
        return value
    raise InputError(f"{path}: {key} must be a boolean.")


def _count(raw: Mapping[str, object], key: str, path: Path) -> int:
    value = raw.get(key)
    if isinstance(value, int) and value >= 0:
        return value
    raise InputError(f"{path}: {key} must be a non-negative integer.")


def _strings(value: object, path: Path, field: str) -> list[str]:
    if isinstance(value, list) and all(isinstance(item, str) for item in value):
        return list(value)
    raise InputError(f"{path}: {field} must be a list of strings.")
