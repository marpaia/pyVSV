"""Bibliographic records, collections, documents, and readings.

Design promotion stays with the design library. This package records where a paper
came from and which files were saved for it.
"""

from vsv.literature.records import (
    PUBMED,
    REVIEW_REFERENCES,
    Citation,
    Collection,
    Document,
    Membership,
    Paper,
    Reading,
)
from vsv.literature.store import Catalog, default_root, load_catalog

__all__ = [
    "PUBMED",
    "REVIEW_REFERENCES",
    "Catalog",
    "Citation",
    "Collection",
    "Document",
    "Membership",
    "Paper",
    "Reading",
    "default_root",
    "load_catalog",
]
