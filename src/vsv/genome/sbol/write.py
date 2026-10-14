"""Write one antigenome as an SBOL3 document."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import sbol3

from vsv.errors import InputError
from vsv.sbol import (
    DNA_ENCODING,
    NAMESPACE,
    PROTEIN_ENCODING,
    ROLE_GENE_END,
    ROLE_GENE_START,
    ROLE_LEADER,
    ROLE_TRAILER,
    SO_CDS,
    SO_GENE,
    SO_GENE_END,
    SO_GENE_START,
    SO_REGION,
    display_id,
    pin,
)

if TYPE_CHECKING:
    from vsv.genome.record import Genome, Span

_ENCODINGS = {
    DNA_ENCODING: sbol3.SBO_DNA,
    PROTEIN_ENCODING: sbol3.SBO_PROTEIN,
}


def write(genome: Genome, path: str | Path) -> None:
    """Write ``genome`` as an SBOL3 genome document."""
    Path(path).write_text(_text(genome))


def _text(genome: Genome) -> str:
    """Turtle for one genome, including leader, genes, signals, coding sequences, and trailer."""
    return _turtle(_document(genome))


def _document(genome: Genome) -> sbol3.Document:
    sbol3.set_namespace(NAMESPACE)
    document = sbol3.Document()
    pinned = display_id(genome.id)
    sequence = sbol3.Sequence(f"{pinned}_sequence", elements=genome.bases, encoding=DNA_ENCODING)
    document.add(sequence)
    component = sbol3.Component(
        pinned,
        [sbol3.SBO_DNA],
        sequences=[sequence],
        name=genome.id,
        description=genome.description or None,
    )
    document.add(component)
    if genome.leader is not None:
        _feature(component, sequence, "leader", "leader", (SO_REGION, ROLE_LEADER), genome.leader)
    for gene in genome.genes:
        _feature(component, sequence, gene.name, gene.name, (SO_GENE,), gene.unit)
        if gene.gene_start is not None:
            _feature(
                component,
                sequence,
                f"{gene.name}_gene_start",
                f"{gene.name} gene start",
                (SO_GENE_START, ROLE_GENE_START),
                gene.gene_start,
            )
        if gene.gene_end is not None:
            _feature(
                component,
                sequence,
                f"{gene.name}_gene_end",
                f"{gene.name} gene end",
                (SO_GENE_END, ROLE_GENE_END),
                gene.gene_end,
            )
        if gene.cds is not None:
            _feature(
                component,
                sequence,
                f"{gene.name}_cds",
                gene.product,
                (SO_CDS,),
                gene.cds,
            )
    if genome.trailer is not None:
        _feature(
            component, sequence, "trailer", "trailer", (SO_REGION, ROLE_TRAILER), genome.trailer
        )
    return document


def write_sequences(records: Sequence[tuple[str, str, str]]) -> str:
    """Turtle for extracted products. Each record is an id, an IUPAC string, and an encoding."""
    sbol3.set_namespace(NAMESPACE)
    document = sbol3.Document()
    seen: set[str] = set()
    for identifier, elements, encoding in records:
        kind = _ENCODINGS.get(encoding)
        if kind is None:
            raise InputError(f"Unsupported sequence encoding {encoding}.")
        if not elements:
            raise InputError(f"{identifier}: sequence is empty.")
        pinned = display_id(identifier)
        if pinned in seen:
            raise InputError(f"Duplicate record identifier: {identifier}")
        seen.add(pinned)
        sequence = sbol3.Sequence(f"{pinned}_sequence", elements=elements, encoding=encoding)
        document.add(sequence)
        document.add(
            sbol3.Component(
                pinned,
                [kind],
                sequences=[sequence],
                name=identifier,
            )
        )
    return _turtle(document)


def _turtle(document: sbol3.Document) -> str:
    text = document.write_string(sbol3.TURTLE)
    if isinstance(text, bytes):
        return text.decode()
    return str(text)


def _feature(
    component: sbol3.Component,
    sequence: sbol3.Sequence,
    identifier: str,
    name: str,
    roles: tuple[str, ...],
    span: Span,
) -> None:
    orientation = sbol3.SBOL_REVERSE_COMPLEMENT if span.strand == "-" else None
    location = sbol3.Range(sequence, span.start + 1, span.end, orientation=orientation)
    feature = sbol3.SequenceFeature([location], roles=list(roles), name=name)
    pin(feature, identifier)
    component.features.append(feature)
