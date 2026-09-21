"""SBOL3 model of VSV genomes.

A genome is an ``sbol3.Component`` that carries one DNA ``sbol3.Sequence``. Its annotations are
``sbol3.SequenceFeature`` objects with ``sbol3.Range`` locations. SBOL coordinates are one-based
and inclusive; :class:`Genome` exposes them as zero-based, half-open spans so they slice Python
strings directly. Sequences are stored as the DNA text of the antigenome, with ``T`` for thymine.

Standard Sequence Ontology roles mark what a feature is. VSV-specific meaning that the ontology
has no term for (leader, trailer, gene-start and gene-end signals) is carried by an additional
role in the ``vsv.bio`` namespace next to the closest SO role.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass

import sbol3

from vsv import NAMESPACE
from vsv.genes import Gene, Product
from vsv.sequence import InputError, reverse_complement, translate

_ROLE_NAMESPACE = f"{NAMESPACE}/role"

ROLE_LEADER = f"{_ROLE_NAMESPACE}/leader"
ROLE_TRAILER = f"{_ROLE_NAMESPACE}/trailer"
ROLE_GENE_START = f"{_ROLE_NAMESPACE}/gene-start"
ROLE_GENE_END = f"{_ROLE_NAMESPACE}/gene-end"

SO_REGION = "https://identifiers.org/SO:0000001"
SO_REGULATORY_REGION = "https://identifiers.org/SO:0005836"
SO_VIRAL_SEQUENCE = "https://identifiers.org/SO:0001041"
SO_GENE = sbol3.SO_GENE
SO_CDS = sbol3.SO_CDS
SO_TERMINATOR = sbol3.SO_TERMINATOR

_SO_TERM = re.compile(r"SO[:_](\d{7})$")
_INVALID_DISPLAY_ID = re.compile(r"[^A-Za-z0-9_]")
_ACCESSION_VERSION = re.compile(r"\.\d+$")
_REVERSE = frozenset({sbol3.SBOL_REVERSE_COMPLEMENT, sbol3.SO_REVERSE})


def role_key(role: str) -> str:
    """Normalize a role URI so that equivalent Sequence Ontology URIs compare equal.

    Sequence Ontology terms appear as ``https://identifiers.org/SO:0000704``,
    ``http://identifiers.org/so/SO:0000704``, and ``http://purl.obolibrary.org/obo/SO_0000704``;
    all normalize to ``SO:0000704``. Other URIs are returned unchanged.
    """
    match = _SO_TERM.search(role)
    return f"SO:{match.group(1)}" if match else role


def new_document() -> sbol3.Document:
    """An empty SBOL3 document whose new objects live in the ``vsv.bio`` namespace."""
    sbol3.set_namespace(NAMESPACE)
    return sbol3.Document()


def display_id(identifier: str) -> str:
    """The SBOL display ID for a sequence identifier.

    The accession version suffix is dropped and characters SBOL does not allow become ``_``, so
    ``NC_001560.1`` is ``NC_001560``.
    """
    cleaned = _INVALID_DISPLAY_ID.sub("_", _ACCESSION_VERSION.sub("", identifier))
    return f"_{cleaned}" if cleaned[:1].isdigit() else cleaned


def _range(feature: sbol3.SequenceFeature) -> sbol3.Range | None:
    return next((loc for loc in feature.locations if isinstance(loc, sbol3.Range)), None)


@dataclass(frozen=True, slots=True)
class GeneAnnotation:
    """The annotated features of one gene in a genome.

    ``gene_feature`` is the transcription unit and ``cds_feature`` its coding sequence. Either may
    be absent, as in curated records that annotate only coding sequences.
    """

    genome: Genome
    gene: Gene
    gene_feature: sbol3.SequenceFeature | None
    cds_feature: sbol3.SequenceFeature | None

    @property
    def gene_span(self) -> tuple[int, int] | None:
        """Zero-based, half-open span of the transcription unit."""
        return None if self.gene_feature is None else self.genome.span(self.gene_feature)

    @property
    def cds_span(self) -> tuple[int, int] | None:
        """Zero-based, half-open span of the coding sequence."""
        return None if self.cds_feature is None else self.genome.span(self.cds_feature)

    @property
    def gene_dna(self) -> str | None:
        """The DNA of the transcription unit."""
        return None if self.gene_feature is None else self.genome.feature_bases(self.gene_feature)

    @property
    def cds_dna(self) -> str | None:
        """The DNA of the coding sequence, including its stop codon."""
        return None if self.cds_feature is None else self.genome.feature_bases(self.cds_feature)

    @property
    def protein(self) -> str | None:
        """The encoded protein without its stop codon."""
        return None if self.cds_feature is None else self.genome.translate(self.cds_feature)

    def sequence(self, product: Product) -> str | None:
        """The requested sequence, or ``None`` if the genome does not annotate it."""
        match product:
            case Product.GENE:
                return self.gene_dna
            case Product.CDS:
                return self.cds_dna
            case Product.PROTEIN:
                return self.protein


class Genome:
    """A genome sequence and its annotations within an SBOL3 document."""

    def __init__(self, document: sbol3.Document, component: sbol3.Component) -> None:
        self.document = document
        self.component = component
        self.sequence = self._dna_sequence()
        self.bases: str = str(self.sequence.elements).upper()

    def _dna_sequence(self) -> sbol3.Sequence:
        for uri in self.component.sequences:
            sequence = self.document.find(uri)
            if (
                isinstance(sequence, sbol3.Sequence)
                and sequence.encoding == sbol3.IUPAC_DNA_ENCODING
            ):
                return sequence
        raise InputError(f"Component {self.component.identity} has no DNA sequence.")

    @property
    def id(self) -> str:
        """The sequence identifier, such as ``NC_001560.1``, or the SBOL display ID."""
        return str(self.component.name or self.component.display_id)

    @property
    def description(self) -> str:
        return str(self.component.description or "")

    def __len__(self) -> int:
        return len(self.bases)

    def features(self, role: str | None = None) -> list[sbol3.SequenceFeature]:
        """Located features ordered by position, optionally only those with ``role``."""
        found = [
            f
            for f in self.component.features
            if isinstance(f, sbol3.SequenceFeature) and _range(f) is not None
        ]
        if role is not None:
            wanted = role_key(role)
            found = [f for f in found if any(role_key(r) == wanted for r in f.roles)]
        return sorted(found, key=lambda f: (self.span(f)[0], -self.span(f)[1]))

    def span(self, feature: sbol3.SequenceFeature) -> tuple[int, int]:
        """Zero-based, half-open ``(start, end)`` of a feature's first range."""
        location = _range(feature)
        if location is None:
            raise ValueError(f"Feature {feature.identity} has no range.")
        return int(location.start) - 1, int(location.end)

    def is_reverse(self, feature: sbol3.SequenceFeature) -> bool:
        """Whether a feature lies on the reverse strand."""
        location = _range(feature)
        orientation = (location.orientation if location else None) or feature.orientation
        return orientation in _REVERSE

    def feature_bases(self, feature: sbol3.SequenceFeature) -> str:
        """The DNA of a feature read in its own direction."""
        start, end = self.span(feature)
        bases = self.bases[start:end]
        return reverse_complement(bases) if self.is_reverse(feature) else bases

    def translate(self, feature: sbol3.SequenceFeature) -> str:
        """The protein encoded by a coding feature, without the stop codon."""
        return translate(self.feature_bases(feature), to_stop=True)

    def named(self, role: str, name: str) -> sbol3.SequenceFeature | None:
        """The first feature with ``role`` and ``name``, if any."""
        return next((f for f in self.features(role) if f.name == name), None)

    def gene(self, name: str) -> sbol3.SequenceFeature | None:
        """The transcription unit named ``name``, if annotated."""
        return self.named(SO_GENE, name)

    def cds(self, name: str) -> sbol3.SequenceFeature | None:
        """The coding sequence of gene ``name``, if annotated."""
        return self.named(SO_CDS, name)

    def protein(self, name: str) -> str | None:
        """The translated product of gene ``name`` without the stop codon."""
        cds = self.cds(name)
        return None if cds is None else self.translate(cds)

    @property
    def genes(self) -> dict[Gene, GeneAnnotation]:
        """The annotated genes, in genome order. Genes with no gene or coding feature are absent."""
        found: dict[Gene, GeneAnnotation] = {}
        for gene in Gene:
            gene_feature, cds_feature = self.gene(gene), self.cds(gene)
            if gene_feature is not None or cds_feature is not None:
                found[gene] = GeneAnnotation(self, gene, gene_feature, cds_feature)
        return found

    def add_region(
        self, name: str, start: int, end: int, roles: Iterable[str]
    ) -> sbol3.SequenceFeature:
        """Add a forward-strand feature over the zero-based, half-open span ``[start, end)``."""
        if not 0 <= start < end <= len(self):
            raise ValueError(f"Span [{start}, {end}) is outside a sequence of length {len(self)}.")
        location = sbol3.Range(self.sequence, start + 1, end, orientation=sbol3.SBOL_INLINE)
        feature = sbol3.SequenceFeature(
            [location],
            roles=list(roles),
            orientation=sbol3.SBOL_INLINE,
            name=name,
        )
        self.component.features.append(feature)
        return feature

    def add_protein(self, cds: sbol3.SequenceFeature, name: str) -> sbol3.Component:
        """Record the protein encoded by ``cds`` as its own component.

        The protein is linked to the coding sequence by a genetic-production interaction.
        """
        base = f"{self.component.display_id}_{display_id(name)}_protein"
        sequence = sbol3.Sequence(
            f"{base}_seq",
            elements=self.translate(cds),
            encoding=sbol3.IUPAC_PROTEIN_ENCODING,
        )
        protein = sbol3.Component(
            base, sbol3.SBO_PROTEIN, name=f"{name} protein", sequences=[sequence]
        )
        self.document.add(sequence)
        self.document.add(protein)
        product = sbol3.SubComponent(protein, name=protein.name)
        self.component.features.append(product)
        self.component.interactions.append(
            sbol3.Interaction(
                [sbol3.SBO_GENETIC_PRODUCTION],
                participations=[
                    sbol3.Participation([sbol3.SBO_TEMPLATE], cds),
                    sbol3.Participation([sbol3.SBO_PRODUCT], product),
                ],
            )
        )
        return protein


def add_genome(document: sbol3.Document, identifier: str, description: str, bases: str) -> Genome:
    """Add a genome with the given sequence to ``document``."""
    key = display_id(identifier)
    sequence = sbol3.Sequence(
        f"{key}_seq", elements=bases.lower(), encoding=sbol3.IUPAC_DNA_ENCODING
    )
    component = sbol3.Component(
        key,
        [sbol3.SBO_DNA, sbol3.SO_LINEAR],
        roles=[SO_VIRAL_SEQUENCE],
        sequences=[sequence],
        name=identifier,
        description=description or None,
    )
    try:
        document.add(sequence)
        document.add(component)
    except ValueError as error:
        raise InputError(f"Duplicate record identifier: {identifier}") from error
    return Genome(document, component)


def genomes(document: sbol3.Document) -> list[Genome]:
    """Every component in ``document`` that has a DNA sequence, ordered by identity."""
    components = sorted(
        (o for o in document.objects if isinstance(o, sbol3.Component)),
        key=lambda c: str(c.identity),
    )
    found: list[Genome] = []
    for component in components:
        try:
            found.append(Genome(document, component))
        except InputError:
            continue
    return found
