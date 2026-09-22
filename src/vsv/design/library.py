"""Load the curated VSV design library.

Each design is a directory with ``design.ttl`` and a provenance note. The document describes
the architecture, its implementations, and qualitative observations. It does not carry
nucleotide or amino-acid sequences.
"""

from __future__ import annotations

from pathlib import Path

import sbol3

from vsv.design.document import design_of
from vsv.design.evidence import Design
from vsv.errors import InputError

LIBRARY_ROOT = Path(__file__).resolve().parents[3] / "data" / "designs"


def _load_ttl(path: Path) -> Design:
    document = sbol3.Document()
    try:
        document.read(str(path))
    except (OSError, ValueError) as error:
        raise InputError(f"{path}: {error}") from error
    return design_of(document)


def load_design(directory: Path) -> Design:
    """Load one design directory."""
    turtle = directory / "design.ttl"
    if not turtle.is_file():
        raise InputError(f"{directory}: missing design.ttl.")
    return _load_ttl(turtle)


def load_library(root: Path | None = None) -> list[Design]:
    """Load every design directory under ``root``, in directory-name order."""
    directory = LIBRARY_ROOT if root is None else root
    if not directory.is_dir():
        raise InputError(f"{directory}: design library not found.")
    directories = sorted(path for path in directory.iterdir() if path.is_dir())
    records = [load_design(path) for path in directories]
    identifiers = [record.architecture.id for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise InputError(f"{directory}: duplicate design ids.")
    return records
