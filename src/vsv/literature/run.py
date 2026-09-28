"""Build the review-reference index. Callers pass the HTTP functions."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from vsv.literature.bibliography import FetchJson, FetchText, build_index
from vsv.literature.catalog import CensusIndex
from vsv.literature.documents import download_xml
from vsv.literature.index import write_index
from vsv.literature.records import Paper

__all__ = ["download_xml", "make_index"]


def make_index(
    dest: Path,
    papers: Mapping[str, Paper],
    fetch_json: FetchJson,
    fetch_text: FetchText,
    email: str,
    retrieved: str,
) -> CensusIndex:
    """Write ``dest/index.json`` and return the index."""
    index = build_index(fetch_json, fetch_text, papers, email, retrieved)
    write_index(dest, index)
    return index
