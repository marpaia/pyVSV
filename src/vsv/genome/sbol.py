"""Read and write a genome as an SBOL3 document.

The file is one antigenome. Callers pass and receive a :class:`~vsv.genome.record.Genome`.
A feature is a single range. Ranges are one-based and inclusive; :class:`~vsv.genome.record.Span`
stays zero-based and half-open. A sequence with no gene feature is annotated from its
transcription signals. A sequence with gene features keeps those spans.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import sbol3

from vsv.errors import InputError
from vsv.genes import Gene
from vsv.genome.annotate import annotate
from vsv.genome.record import GeneSpan, Genome, Span
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

_ENCODINGS = {
    DNA_ENCODING: sbol3.SBO_DNA,
    PROTEIN_ENCODING: sbol3.SBO_PROTEIN,
}


@dataclass(frozen=True, slots=True)
class _Piece:
    kind: str
    name: str
    span: Span


def read_genomes(path: str | Path) -> list[Genome]:
    """Every DNA component in an SBOL genome document that owns one sequence."""
    path = Path(path)
    document = _open(path)
    genomes: list[Genome] = []
    for obj in document.objects:
        if isinstance(obj, sbol3.Component) and _is_dna(obj):
            sequence = _sequence(document, obj, path)
            if sequence is None:
                continue
            genomes.append(_genome(obj, sequence, path))
    if not genomes:
        raise InputError(f"{path}: no genome sequence.")
    return genomes


def write_genome(genome: Genome) -> str:
    """Turtle for one genome, including leader, genes, signals, coding sequences, and trailer."""
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
    return _turtle(document)


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


def _open(path: Path) -> sbol3.Document:
    sbol3.set_namespace(NAMESPACE)
    document = sbol3.Document()
    try:
        document.read(str(path))
    except OSError as error:
        raise InputError(f"{path}: {error.strerror or error}") from error
    except Exception as error:
        raise InputError(f"{path}: not a readable SBOL file.") from error
    return document


def _turtle(document: sbol3.Document) -> str:
    text = document.write_string(sbol3.TURTLE)
    if isinstance(text, bytes):
        return text.decode()
    return str(text)


def _is_dna(component: sbol3.Component) -> bool:
    return sbol3.SBO_DNA in {str(kind) for kind in component.types}


def _sequence(
    document: sbol3.Document, component: sbol3.Component, path: Path
) -> sbol3.Sequence | None:
    uris = [str(uri) for uri in component.sequences]
    if not uris:
        return None
    if len(uris) != 1:
        raise InputError(f"{path}: {component.display_id} owns more than one sequence.")
    sequence = document.find(uris[0])
    if not isinstance(sequence, sbol3.Sequence) or not sequence.elements:
        raise InputError(f"{path}: {component.display_id} is missing its sequence.")
    if str(sequence.encoding) != DNA_ENCODING:
        raise InputError(f"{path}: {component.display_id} is not an IUPAC DNA sequence.")
    return sequence


def _genome(component: sbol3.Component, sequence: sbol3.Sequence, path: Path) -> Genome:
    identifier = str(component.name) if component.name else str(component.display_id)
    description = str(component.description) if component.description else ""
    bases = str(sequence.elements).upper().replace("U", "T")
    pieces = _pieces(component, path)
    if not any(piece.kind == "gene" for piece in pieces):
        return annotate(bases, identifier=identifier, description=description)
    return Genome(
        identifier,
        description,
        bases,
        _single(pieces, "leader", path),
        _single(pieces, "trailer", path),
        tuple(_genes(pieces, path)),
    )


def _pieces(component: sbol3.Component, path: Path) -> list[_Piece]:
    pieces: list[_Piece] = []
    for feature in component.features:
        if not isinstance(feature, sbol3.SequenceFeature):
            raise InputError(f"{path}: {feature.display_id} is not a sequence feature.")
        locations = list(feature.locations)
        if len(locations) != 1 or not isinstance(locations[0], sbol3.Range):
            raise InputError(f"{path}: {feature.display_id} has no single range.")
        kind = _kind(feature)
        if kind is None:
            raise InputError(f"{path}: {feature.display_id} has no genome role.")
        name = str(feature.name) if feature.name else str(feature.display_id)
        pieces.append(_Piece(kind, name, _span(locations[0], path, str(feature.display_id))))
    return pieces


def _kind(feature: sbol3.Feature) -> str | None:
    roles = {str(role) for role in feature.roles}
    if ROLE_LEADER in roles:
        return "leader"
    if ROLE_TRAILER in roles:
        return "trailer"
    if SO_CDS in roles:
        return "cds"
    if ROLE_GENE_START in roles:
        return "gene_start"
    if ROLE_GENE_END in roles:
        return "gene_end"
    if SO_GENE in roles:
        return "gene"
    return None


def _span(location: sbol3.Range, path: Path, name: str) -> Span:
    orientation = str(location.orientation or "")
    if orientation in {"", sbol3.SBOL_INLINE}:
        strand = "+"
    elif orientation == sbol3.SBOL_REVERSE_COMPLEMENT:
        strand = "-"
    else:
        raise InputError(f"{path}: {name} has an unknown orientation.")
    try:
        return Span(int(location.start) - 1, int(location.end), strand)
    except InputError as error:
        raise InputError(f"{path}: {name}: {error}") from error


def _single(pieces: Sequence[_Piece], kind: str, path: Path) -> Span | None:
    found = [piece for piece in pieces if piece.kind == kind]
    if len(found) > 1:
        raise InputError(f"{path}: expected one {kind}.")
    return None if not found else found[0].span


def _genes(pieces: Sequence[_Piece], path: Path) -> list[GeneSpan]:
    genes = sorted(
        (piece for piece in pieces if piece.kind == "gene"),
        key=lambda piece: (piece.span.start, piece.span.end),
    )
    symbols: Sequence[Gene | None] = list(Gene) if len(genes) == len(Gene) else [None] * len(genes)
    claimed: dict[int, str] = {}
    built: list[GeneSpan] = []
    for gene, symbol in zip(genes, symbols, strict=True):
        cds = _claim(pieces, "cds", gene, path, claimed)
        start = _claim(pieces, "gene_start", gene, path, claimed)
        end = _claim(pieces, "gene_end", gene, path, claimed)
        if len(start) > 1 or len(end) > 1:
            raise InputError(f"{path}: {gene.name} has more than one transcription signal.")
        chosen = max(cds, key=lambda piece: (len(piece.span), -piece.span.start), default=None)
        built.append(
            GeneSpan(
                gene.name,
                symbol,
                gene.span,
                None if chosen is None else chosen.span,
                gene.name if chosen is None else chosen.name,
                None if not start else start[0].span,
                None if not end else end[0].span,
            )
        )
    outside = [
        piece
        for piece in pieces
        if piece.kind in {"cds", "gene_start", "gene_end"} and id(piece) not in claimed
    ]
    if outside:
        raise InputError(f"{path}: {outside[0].name} is outside every gene.")
    return built


def _claim(
    pieces: Sequence[_Piece], kind: str, gene: _Piece, path: Path, claimed: dict[int, str]
) -> list[_Piece]:
    found = []
    for piece in pieces:
        if piece.kind != kind or not (gene.span.start <= piece.span.start < gene.span.end):
            continue
        owner = claimed.get(id(piece))
        if owner is not None:
            raise InputError(f"{path}: {piece.name} falls in {owner} and {gene.name}.")
        claimed[id(piece)] = gene.name
        found.append(piece)
    return found


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
