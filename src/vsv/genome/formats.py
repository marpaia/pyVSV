"""Reading genomes from FASTA, GenBank, and SBOL3 files into one SBOL3 document."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from xml.sax import SAXParseException

import sbol3
from sbol2 import SBOLError
from sbol_utilities.conversion import convert_from_genbank

import vsv.genome.sbol as sbol
from vsv import NAMESPACE
from vsv.errors import InputError
from vsv.genome.annotation import name_curated_genes
from vsv.genome.sbol import Genome
from vsv.genome.sequence import read_fasta

FASTA_SUFFIXES = frozenset({".fasta", ".fa", ".fna", ".fas"})
GENBANK_SUFFIXES = frozenset({".gb", ".gbk", ".gbff", ".genbank"})
SBOL_SUFFIXES = {
    ".ttl": sbol3.TURTLE,
    ".nt": sbol3.NTRIPLES,
    ".rdf": sbol3.RDF_XML,
    ".xml": sbol3.RDF_XML,
    ".jsonld": sbol3.JSONLD,
}

# Names accepted for writing SBOL3, mapped to the pySBOL3 serialization.
OUTPUT_FORMATS = {
    "turtle": sbol3.TURTLE,
    "ntriples": sbol3.NTRIPLES,
    "rdfxml": sbol3.RDF_XML,
    "jsonld": sbol3.JSONLD,
}


def read_genomes(
    paths: Sequence[Path], *, name_genes: bool = True
) -> tuple[sbol3.Document, list[Genome]]:
    """Read every genome from ``paths`` into a single SBOL3 document.

    The format of each file is chosen by its suffix: FASTA, GenBank, or SBOL3 (Turtle,
    N-Triples, RDF/XML, or JSON-LD). Genomes are returned in input order. When ``name_genes``
    is set, a GenBank or SBOL record with five genes is named as a VSV genome.
    """
    document = sbol.new_document()
    found: list[Genome] = []
    for path in paths:
        suffix = path.suffix.lower()
        try:
            if suffix in FASTA_SUFFIXES:
                loaded = _read_fasta(document, path)
            elif suffix in GENBANK_SUFFIXES:
                loaded = _read_genbank(document, path, name_genes=name_genes)
            elif suffix in SBOL_SUFFIXES:
                loaded = _read_sbol(document, path, name_genes=name_genes)
            else:
                raise InputError(
                    f"Unrecognized file type {suffix or '(none)'}; expected FASTA "
                    f"({', '.join(sorted(FASTA_SUFFIXES))}), GenBank "
                    f"({', '.join(sorted(GENBANK_SUFFIXES))}), or SBOL3 "
                    f"({', '.join(sorted(SBOL_SUFFIXES))})."
                )
        except OSError as error:
            raise InputError(f"{path}: {error.strerror or error}") from error
        except InputError as error:
            if str(error).startswith(f"{path}:"):
                raise
            raise InputError(f"{path}: {error}") from error
        found.extend(loaded)
    if not found:
        raise InputError("No sequence records were found.")
    return document, found


def _read_fasta(document: sbol3.Document, path: Path) -> list[Genome]:
    return [
        sbol.add_genome(document, record.id, record.description, record.sequence)
        for record in read_fasta(path)
    ]


def _read_genbank(document: sbol3.Document, path: Path, *, name_genes: bool) -> list[Genome]:
    try:
        converted = convert_from_genbank(str(path), NAMESPACE)
    except ValueError as error:
        raise InputError(f"not a readable GenBank file: {error}") from error
    except SBOLError as error:
        # The message describes SBOL2 validation, which does not help with a GenBank file.
        raise InputError("not a readable GenBank file.") from error
    return _merge(document, converted, name_genes=name_genes)


def _read_sbol(document: sbol3.Document, path: Path, *, name_genes: bool) -> list[Genome]:
    incoming = sbol3.Document()
    try:
        incoming.read(str(path), SBOL_SUFFIXES[path.suffix.lower()])
    except (ValueError, SyntaxError, SAXParseException) as error:
        raise InputError(f"not a readable SBOL3 file: {error}") from error
    return _merge(document, incoming, name_genes=name_genes)


def _merge(document: sbol3.Document, incoming: sbol3.Document, *, name_genes: bool) -> list[Genome]:
    """Move every object of ``incoming`` into ``document`` and return its genomes."""
    names = [str(o.identity) for o in incoming.objects]
    duplicates = [name for name in names if document.find(name) is not None]
    if duplicates:
        raise InputError(f"Duplicate record identifier: {duplicates[0]}")
    for obj in list(incoming.objects):
        document.add(obj)
    merged = [g for g in sbol.genomes(document) if str(g.component.identity) in names]
    if name_genes:
        for genome in merged:
            name_curated_genes(genome)
    return merged


def write_document(document: sbol3.Document, output_format: str) -> str:
    """Serialize a document as SBOL3 text in one of :data:`OUTPUT_FORMATS`."""
    return str(document.write_string(OUTPUT_FORMATS[output_format]))
