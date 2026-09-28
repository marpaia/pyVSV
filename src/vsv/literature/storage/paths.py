"""Workspace location and the filenames of stored literature.

This module lives at ``src/vsv/literature/storage/paths.py``. ``parents[4]`` is the
repository root, so the default workspace stays ``<repository>/data``.
"""

from __future__ import annotations

from pathlib import Path

_REPOSITORY = Path(__file__).resolve().parents[4]

PAPERS_NAME = "papers.jsonl"
EXTRACTIONS_NAME = "extractions.jsonl"
LEDGER_NAME = "ledger.json"
AUDIT_NAME = "audit.json"
FULLTEXT_NAME = "fulltext"
INDEX_NAME = "index.json"
MANIFEST_NAME = "manifest.jsonl"


def default_root() -> Path:
    """The data directory that holds the corpus and the census."""
    return _REPOSITORY / "data"


def corpus_dir(root: Path) -> Path:
    """PubMed collection files under ``root``."""
    return root / "corpus"


def census_dir(root: Path) -> Path:
    """Review-reference files under ``root``."""
    return root / "census"


def fulltext_dir(root: Path) -> Path:
    """Open-access XML files for the workspace."""
    return corpus_dir(root) / FULLTEXT_NAME


def acquired_dir(root: Path) -> Path:
    """Publisher files saved for review references."""
    return census_dir(root) / "acquired"


def sections_dir(root: Path) -> Path:
    """Parsed document text for the workspace."""
    return census_dir(root) / "sections"


CORPUS_ROOT = corpus_dir(default_root())
CENSUS_ROOT = census_dir(default_root())
