"""Bibliographic records, collections, documents, and parsed text.

Design promotion stays with the design library. This package records where a paper
came from and which files were saved for it.
"""

from vsv.literature.catalog import Catalog, load_catalog
from vsv.literature.records import (
    PUBMED,
    REVIEW_REFERENCES,
    Citation,
    Collection,
    Document,
    DocumentText,
    Membership,
    Paper,
)
from vsv.literature.storage.paths import default_root

__all__ = [
    "PUBMED",
    "REVIEW_REFERENCES",
    "Catalog",
    "Citation",
    "Collection",
    "Document",
    "DocumentText",
    "Membership",
    "Paper",
    "default_root",
    "load_catalog",
]
