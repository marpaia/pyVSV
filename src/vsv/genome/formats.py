"""FASTA, GenBank, and SBOL genome readers.

FASTA is annotated from transcription signals. GenBank locations are simple spans. A ``join``
or ``complement`` location is rejected. Five ``gene`` features, in file order, are
nucleoprotein, phosphoprotein, matrix protein, glycoprotein, and large polymerase. The longest
coding sequence inside each gene is that gene's CDS, so the shorter overlapping ORF inside
phosphoprotein is not the phosphoprotein. With no gene features and exactly five CDS features,
those five are named in file order. An SBOL document is read by
:func:`vsv.genome.sbol.read.read_genomes`.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

from vsv.errors import InputError
from vsv.genes import Gene
from vsv.genome.annotate import annotate
from vsv.genome.record import GeneSpan, Genome, Span
from vsv.genome.sbol.read import read_genomes
from vsv.genome.sequence import read_fasta

FASTA_SUFFIXES = frozenset({".fasta", ".fa", ".fna", ".fas"})
GENBANK_SUFFIXES = frozenset({".gb", ".gbk", ".gbff", ".genbank"})
SBOL_SUFFIXES = frozenset({".ttl"})

_SPAN = re.compile(r"^(?P<lt><)?(?P<start>\d+)\.\.(?P<gt>>)?(?P<end>\d+)$")
_POINT = re.compile(r"^(?P<lt><)?(?P<point>\d+)(?P<gt>>)?$")
_LOCUS_LENGTH = re.compile(r"^LOCUS\s+\S+\s+(\d+)\s+bp\b")


def read_records(paths: Sequence[Path]) -> list[Genome]:
    """Read every genome in ``paths``, in input order, with VSV genes identified."""
    found: list[Genome] = []
    seen: set[str] = set()
    for path in paths:
        try:
            loaded = _read_one(path)
        except OSError as error:
            raise InputError(f"{path}: {error.strerror or error}") from error
        except InputError as error:
            if str(error).startswith(f"{path}:"):
                raise
            raise InputError(f"{path}: {error}") from error
        for genome in loaded:
            if genome.id in seen:
                raise InputError(f"Duplicate record identifier: {genome.id}")
            seen.add(genome.id)
            found.append(genome)
    if not found:
        raise InputError("No sequence records were found.")
    return found


def _read_one(path: Path) -> list[Genome]:
    suffix = path.suffix.lower()
    if suffix in FASTA_SUFFIXES:
        return [
            annotate(record.sequence, identifier=record.id, description=record.description)
            for record in read_fasta(path)
        ]
    if suffix in GENBANK_SUFFIXES:
        return _read_genbank(path)
    if suffix in SBOL_SUFFIXES:
        return read_genomes(path)
    raise InputError(
        f"Unrecognized file type {suffix or '(none)'}; expected FASTA "
        f"({', '.join(sorted(FASTA_SUFFIXES))}), GenBank "
        f"({', '.join(sorted(GENBANK_SUFFIXES))}), or SBOL (.ttl)."
    )


def _read_genbank(path: Path) -> list[Genome]:
    try:
        text = path.read_text(encoding="ascii")
    except UnicodeDecodeError as error:
        raise InputError(f"{path}: file is not ASCII text.") from error
    records = _genbank_records(text, path)
    if not records:
        raise InputError(f"{path}: not a readable GenBank file.")
    return records


class _Feature:
    def __init__(self, key: str, location: str) -> None:
        self.key = key
        self.location = location
        self.qualifiers: dict[str, str] = {}


def _genbank_records(text: str, path: Path) -> list[Genome]:
    records: list[Genome] = []
    identifier = ""
    description = ""
    locus_length: int | None = None
    version = ""
    accession = ""
    mode = "header"
    definition = False
    features: list[_Feature] = []
    sequence: list[str] = []
    saw_locus = False

    def flush() -> None:
        nonlocal identifier, description, locus_length, version, accession, mode, features, sequence
        nonlocal saw_locus, definition
        if not saw_locus:
            return
        bases = "".join(sequence).upper().replace("U", "T")
        if not bases:
            raise InputError(f"{path}: not a readable GenBank file.")
        if locus_length is not None and locus_length != len(bases):
            raise InputError(f"{path}: LOCUS length {locus_length} is not {len(bases)}.")
        record_id = version or accession or identifier
        records.append(_genome_from_features(record_id, description, bases, features, path))
        identifier = description = version = accession = ""
        locus_length = None
        mode = "header"
        definition = False
        features = []
        sequence = []
        saw_locus = False

    for raw in text.splitlines():
        line = raw.rstrip()
        if line.startswith("//"):
            flush()
            continue
        if mode == "origin":
            sequence.append(re.sub(r"[^A-Za-z]", "", line))
            continue
        if line.startswith("ORIGIN"):
            mode = "origin"
            definition = False
            continue
        if line.startswith("FEATURES"):
            mode = "features"
            definition = False
            continue
        if mode == "header":
            keyword = line.split()[0] if line[:1].isalpha() else ""
            if keyword == "LOCUS":
                saw_locus = True
                identifier = line.split()[1]
                matched = _LOCUS_LENGTH.match(line)
                locus_length = int(matched.group(1)) if matched else None
                definition = False
            elif keyword == "DEFINITION":
                description = line.split(None, 1)[1].strip().removesuffix(".")
                definition = True
            elif keyword == "ACCESSION":
                accession = line.split()[1]
                definition = False
            elif keyword == "VERSION":
                version = line.split()[1]
                definition = False
            elif definition and line[:1].isspace():
                description = f"{description} {line.strip()}".removesuffix(".")
            else:
                definition = False
            continue
        if len(line) > 5 and line.startswith("     ") and line[5] != " ":
            parts = line.strip().split(None, 1)
            if len(parts) == 2:
                features.append(_Feature(parts[0], parts[1].strip()))
            continue
        text_line = line.strip()
        if not text_line or not features:
            continue
        if text_line.startswith("/") and "=" in text_line:
            qual, _, value = text_line[1:].partition("=")
            features[-1].qualifiers[qual] = _unquote(value)
            continue
        if features[-1].qualifiers:
            qual = next(reversed(features[-1].qualifiers))
            features[-1].qualifiers[qual] += _unquote(text_line)
    if saw_locus:
        flush()
    return records


def _unquote(value: str) -> str:
    value = value.strip()
    if value.startswith('"') and value.endswith('"') and len(value) >= 2:
        return value[1:-1]
    if value.startswith('"'):
        return value[1:]
    if value.endswith('"'):
        return value[:-1]
    return value


def _location(feature: _Feature, path: Path) -> Span:
    location = feature.location.replace(" ", "")
    if any(token in location for token in ("join", "complement", "order")):
        name = feature.qualifiers.get("product") or feature.qualifiers.get("gene") or feature.key
        raise InputError(f"{path}: feature {name} uses a {location.split('(')[0]} location.")
    span = _SPAN.fullmatch(location)
    if span is not None:
        return Span(int(span.group("start")) - 1, int(span.group("end")))
    point = _POINT.fullmatch(location)
    if point is not None:
        start = int(point.group("point")) - 1
        return Span(start, start + 1)
    raise InputError(
        f"{path}: feature {feature.key} has an unreadable location {feature.location}."
    )


def _genome_from_features(
    identifier: str,
    description: str,
    bases: str,
    features: list[_Feature],
    path: Path,
) -> Genome:
    genes = [feature for feature in features if feature.key == "gene"]
    cdss = [feature for feature in features if feature.key == "CDS"]
    if len(genes) == len(Gene):
        return _from_genes(identifier, description, bases, genes, cdss, path)
    if not genes and len(cdss) == len(Gene):
        return _from_coding_sequences(identifier, description, bases, cdss, path)
    if not genes and not cdss:
        return annotate(bases, identifier=identifier, description=description)
    return Genome(identifier, description, bases)


def _from_genes(
    identifier: str,
    description: str,
    bases: str,
    genes: list[_Feature],
    cdss: list[_Feature],
    path: Path,
) -> Genome:
    coding = [(feature, _location(feature, path)) for feature in cdss]
    named: list[GeneSpan] = []
    for symbol, feature in zip(Gene, genes, strict=True):
        unit = _location(feature, path)
        inside = [
            (item, span)
            for item, span in coding
            if unit.start <= span.start and span.end <= unit.end
        ]
        chosen = max(inside, key=lambda item: len(item[1]), default=None)
        cds = None if chosen is None else chosen[1]
        product = (
            symbol.value if chosen is None else chosen[0].qualifiers.get("product", symbol.value)
        )
        named.append(GeneSpan(symbol.value, symbol, unit, cds, product))
    return Genome(identifier, description, bases, genes=tuple(named))


def _from_coding_sequences(
    identifier: str,
    description: str,
    bases: str,
    cdss: list[_Feature],
    path: Path,
) -> Genome:
    named: list[GeneSpan] = []
    for symbol, feature in zip(Gene, cdss, strict=True):
        span = _location(feature, path)
        product = feature.qualifiers.get("product", f"{symbol.value} CDS")
        named.append(GeneSpan(symbol.value, symbol, span, span, product))
    return Genome(identifier, description, bases, genes=tuple(named))
