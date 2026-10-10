"""SBOL identifiers and display ids shared by a vaccine plan and a genome document."""

from __future__ import annotations

import re

import sbol3

NAMESPACE = "https://vsv.bio"
ROLE = f"{NAMESPACE}/role"
ROLE_LEADER = f"{ROLE}/leader"
ROLE_TRAILER = f"{ROLE}/trailer"
ROLE_GENE_START = f"{ROLE}/gene-start"
ROLE_GENE_END = f"{ROLE}/gene-end"

SO_REGION = "https://identifiers.org/SO:0000001"
SO_GENE = sbol3.SO_GENE
SO_CDS = sbol3.SO_CDS
SO_GENE_START = "https://identifiers.org/SO:0005836"
SO_GENE_END = sbol3.SO_TERMINATOR

DNA_ENCODING = sbol3.IUPAC_DNA_ENCODING
PROTEIN_ENCODING = sbol3.IUPAC_PROTEIN_ENCODING

_ACCESSION_VERSION = re.compile(r"\.\d+\Z")
_INVALID_DISPLAY_ID = re.compile(r"[^A-Za-z0-9_]")

sbol3.set_namespace(NAMESPACE)


def display_id(identifier: str) -> str:
    """The pySBOL3 display id for ``identifier``.

    A trailing ``.digits`` version is removed. Any other character outside ASCII letters,
    digits, and ``_`` becomes ``_``.
    """
    cleaned = _INVALID_DISPLAY_ID.sub("_", _ACCESSION_VERSION.sub("", identifier))
    if cleaned[:1].isdigit():
        cleaned = f"_{cleaned}"
    return cleaned


def pin(identified: sbol3.Identified, identifier: str) -> None:
    """Set the display id pySBOL3 keeps when this object is parented."""
    identified._display_id = display_id(identifier)
