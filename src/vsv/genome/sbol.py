"""SBOL3 model of VSV genomes and their parts.

A genome is an ``sbol3.Component`` that carries one DNA ``sbol3.Sequence``. Located annotations
are ``sbol3.SequenceFeature`` or ``sbol3.SubComponent`` objects with ``sbol3.Range`` locations.
SBOL coordinates are one-based and inclusive; :class:`Genome` exposes them as zero-based,
half-open spans so they slice Python strings directly. Sequences are stored as the DNA text of
the antigenome, with ``T`` for thymine.

A transcription unit, its coding sequence, and its protein are different parts. The transcription
unit carries the gene name (``Nucleoprotein``). The coding sequence is its own DNA component,
located on the genome as ``Nucleoprotein_cds``, and it keeps the product name
(``nucleocapsid (N) protein``) rather than repeating the gene name. The same split names a coding
sequence we infer from the antigenome ``Nucleoprotein CDS``. A role
``https://vsv.bio/gene/Nucleoprotein`` records which gene a part belongs to.

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
from vsv.errors import InputError
from vsv.genes import Gene, Product
from vsv.genome.sequence import reverse_complement, translate

sbol3.set_namespace(NAMESPACE)

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

type Part = sbol3.SequenceFeature | sbol3.SubComponent

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


def gene_symbol_role(symbol: str) -> str:
    """The role URI for a VSV gene, for example ``https://vsv.bio/gene/Nucleoprotein``."""
    return f"{NAMESPACE}/gene/{display_id(symbol)}"


def part_id(symbol: str, kind: str) -> str:
    """The display id of a gene's part, such as ``Nucleoprotein_cds``."""
    return f"{display_id(symbol)}_{kind}"


def pin_display_id(identified: sbol3.Identified, identifier: str) -> None:
    """Set the display id pySBOL3 keeps when this object is parented.

    The display id property has no setter. An owned object that already has one keeps it;
    otherwise pySBOL3 assigns a numbered id.
    """
    if display_id(identifier) != identifier:
        raise InputError(f"Invalid part id {identifier!r}.")
    identified._display_id = identifier  # pySBOL3 reads this private attribute when parenting


def ensure_role(feature: Part, role: str) -> None:
    """Add ``role`` unless the feature already carries an equivalent one."""
    if any(role_key(existing) == role_key(role) for existing in feature.roles):
        return
    feature.roles.append(role)


def _range(feature: Part) -> sbol3.Range | None:
    return next((loc for loc in feature.locations if isinstance(loc, sbol3.Range)), None)


def _part_roles(feature: Part, document: sbol3.Document) -> list[str]:
    """Roles of a located part, including those of the component it instances."""
    roles = [str(role) for role in feature.roles]
    if isinstance(feature, sbol3.SubComponent):
        instance = feature.instance_of
        if not isinstance(instance, sbol3.Component):
            found = document.find(str(instance))
            instance = found if isinstance(found, sbol3.Component) else None
        if isinstance(instance, sbol3.Component):
            roles.extend(str(role) for role in instance.roles)
    return roles


def _symbol_of(feature: Part, document: sbol3.Document) -> str | None:
    prefix = f"{NAMESPACE}/gene/"
    for role in _part_roles(feature, document):
        if role.startswith(prefix):
            return role[len(prefix) :]
    return None


@dataclass(frozen=True, slots=True)
class GeneAnnotation:
    """The annotated features of one gene in a genome.

    ``gene_feature`` is the transcription unit, named with the gene symbol. ``cds_feature`` is the
    coding sequence of that gene, a separate part. Either may be absent, as in curated records
    that annotate only coding sequences.
    """

    genome: Genome
    gene: Gene
    gene_feature: Part | None
    cds_feature: Part | None

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

    def features(self, role: str | None = None) -> list[Part]:
        """Located features ordered by position, optionally only those with ``role``.

        A coding sequence located as a subcomponent is included. A protein subcomponent has no
        range and is not a genomic feature.
        """
        found = [
            feature
            for feature in self.component.features
            if isinstance(feature, (sbol3.SequenceFeature, sbol3.SubComponent))
            and _range(feature) is not None
        ]
        if role is not None:
            wanted = role_key(role)
            found = [f for f in found if any(role_key(r) == wanted for r in f.roles)]
        return sorted(found, key=lambda f: (self.span(f)[0], -self.span(f)[1]))

    def span(self, feature: Part) -> tuple[int, int]:
        """Zero-based, half-open ``(start, end)`` of a feature's first range."""
        location = _range(feature)
        if location is None:
            raise ValueError(f"Feature {feature.identity} has no range.")
        return int(location.start) - 1, int(location.end)

    def is_reverse(self, feature: Part) -> bool:
        """Whether a feature lies on the reverse strand."""
        location = _range(feature)
        orientation = (location.orientation if location else None) or feature.orientation
        return orientation in _REVERSE

    def feature_bases(self, feature: Part) -> str:
        """The DNA of a feature read in its own direction."""
        start, end = self.span(feature)
        bases = self.bases[start:end]
        return reverse_complement(bases) if self.is_reverse(feature) else bases

    def translate(self, feature: Part) -> str:
        """The protein encoded by a coding feature, without the stop codon."""
        return translate(self.feature_bases(feature), to_stop=True)

    def named(self, role: str, name: str) -> Part | None:
        """The first feature with ``role`` and ``name``, if any."""
        return next((f for f in self.features(role) if f.name == name), None)

    def gene(self, name: str) -> Part | None:
        """The transcription unit of gene ``name``, if annotated."""
        symbol = display_id(name)
        for feature in self.features(SO_GENE):
            if (
                feature.display_id == symbol
                or feature.name == name
                or _symbol_of(feature, self.document) == symbol
            ):
                return feature
        return None

    def cds(self, name: str) -> Part | None:
        """The coding sequence of gene ``name``, if annotated.

        This is the part ``{name}_cds``, or the coding sequence marked with that gene's symbol.
        A feature merely named ``name`` is used only when no part is linked to the gene. The
        longest coding sequence inside the transcription unit wins when several overlap it.
        """
        symbol = display_id(name)
        coded = self.features(SO_CDS)
        identified = next(
            (feature for feature in coded if feature.display_id == part_id(symbol, "cds")),
            None,
        )
        if identified is not None:
            return identified
        marked = [feature for feature in coded if _symbol_of(feature, self.document) == symbol]
        if marked:
            return max(marked, key=lambda feature: self.span(feature)[1] - self.span(feature)[0])
        gene = self.gene(name)
        if gene is not None:
            start, end = self.span(gene)
            inside = [
                feature
                for feature in coded
                if start <= self.span(feature)[0] and self.span(feature)[1] <= end
            ]
            if inside:
                return max(
                    inside,
                    key=lambda feature: self.span(feature)[1] - self.span(feature)[0],
                )
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
        self, name: str, start: int, end: int, roles: Iterable[str], *, part: str | None = None
    ) -> sbol3.SequenceFeature:
        """Add a forward-strand feature over the zero-based, half-open span ``[start, end)``.

        ``part`` is the display id. Without it, pySBOL3 assigns a numbered id.
        """
        if not 0 <= start < end <= len(self):
            raise ValueError(f"Span [{start}, {end}) is outside a sequence of length {len(self)}.")
        location = sbol3.Range(self.sequence, start + 1, end, orientation=sbol3.SBOL_INLINE)
        feature = sbol3.SequenceFeature(
            [location],
            roles=list(roles),
            orientation=sbol3.SBOL_INLINE,
            name=name,
        )
        if part is not None:
            pin_display_id(feature, part)
        self.component.features.append(feature)
        return feature

    def identify_gene(self, feature: Part, symbol: str) -> Part:
        """Name a transcription unit with gene ``symbol`` and give it that part id.

        A previous name other than the feature key ``gene`` is kept as the description.
        """
        if feature.name and feature.name not in {symbol, "gene"} and not feature.description:
            feature.description = feature.name
        feature.name = symbol
        ensure_role(feature, gene_symbol_role(symbol))
        if isinstance(feature, sbol3.SequenceFeature):
            return self._reidentify(feature, display_id(symbol))
        return feature

    def add_cds(
        self, symbol: str, start: int, end: int, *, name: str | None = None
    ) -> sbol3.SubComponent:
        """Add gene ``symbol``'s coding sequence as a DNA component located on this genome.

        The located part's display id is ``{symbol}_cds``. ``name`` is the product name and
        defaults to ``{symbol} CDS``. It is not the gene symbol.
        """
        if not 0 <= start < end <= len(self):
            raise ValueError(f"Span [{start}, {end}) is outside a sequence of length {len(self)}.")
        product = name or f"{symbol} CDS"
        component = self._cds_component(symbol, product, self.bases[start:end])
        return self._locate_cds(component, symbol, product, start, end, sbol3.SBOL_INLINE)

    def define_cds(self, feature: Part, symbol: str) -> sbol3.SubComponent:
        """Make ``feature`` the located coding-sequence component of gene ``symbol``.

        The product name already on ``feature`` is kept. A second call returns the same part.
        """
        identifier = part_id(symbol, "cds")
        if isinstance(feature, sbol3.SubComponent) and feature.display_id == identifier:
            return feature
        start, end = self.span(feature)
        product = str(feature.name or f"{symbol} CDS")
        orientation = (
            sbol3.SBOL_REVERSE_COMPLEMENT if self.is_reverse(feature) else sbol3.SBOL_INLINE
        )
        component = self._cds_component(symbol, product, self.feature_bases(feature))
        located = self._locate_cds(component, symbol, product, start, end, orientation)
        if feature.description and not located.description:
            located.description = feature.description
        if feature in self.component.features:
            self.component.features.remove(feature)
        return located

    def constrain(self, subject: Part, target: Part, restriction: str) -> None:
        """Record ``restriction`` from ``subject`` to ``target`` unless it is already present."""
        key = (str(subject.identity), str(target.identity), restriction)
        existing = {
            (str(constraint.subject), str(constraint.object), str(constraint.restriction))
            for constraint in self.component.constraints
        }
        if key in existing:
            return
        self.component.constraints.append(sbol3.Constraint(restriction, subject, target))

    def _cds_component(self, symbol: str, product: str, bases: str) -> sbol3.Component:
        identifier = f"{self.component.display_id}_{part_id(symbol, 'cds')}"
        found = self.document.find(f"{NAMESPACE}/{identifier}")
        if isinstance(found, sbol3.Component):
            return found
        sequence = sbol3.Sequence(
            f"{identifier}_seq", elements=bases.lower(), encoding=sbol3.IUPAC_DNA_ENCODING
        )
        component = sbol3.Component(
            identifier,
            [sbol3.SBO_DNA],
            roles=[SO_CDS, gene_symbol_role(symbol)],
            name=product,
            sequences=[sequence],
        )
        try:
            self.document.add(sequence)
            self.document.add(component)
        except ValueError as error:
            raise InputError(f"Duplicate coding sequence: {symbol}") from error
        return component

    def _locate_cds(
        self,
        component: sbol3.Component,
        symbol: str,
        product: str,
        start: int,
        end: int,
        orientation: str,
    ) -> sbol3.SubComponent:
        identifier = part_id(symbol, "cds")
        existing = next(
            (
                feature
                for feature in self.component.features
                if isinstance(feature, sbol3.SubComponent) and feature.display_id == identifier
            ),
            None,
        )
        if existing is not None:
            return existing
        location = sbol3.Range(self.sequence, start + 1, end, orientation=orientation)
        # A subcomponent has a single role. The gene symbol is on the CDS component.
        located = sbol3.SubComponent(
            component,
            role_integration=SO_CDS,
            orientation=orientation,
            name=product,
            locations=[location],
        )
        pin_display_id(located, identifier)
        self.component.features.append(located)
        return located

    def _reidentify(self, feature: sbol3.SequenceFeature, identifier: str) -> sbol3.SequenceFeature:
        if feature.display_id == identifier or self.is_reverse(feature):
            return feature
        start, end = self.span(feature)
        replacement = self.add_region(
            str(feature.name or identifier),
            start,
            end,
            list(feature.roles),
            part=identifier,
        )
        if feature.description:
            replacement.description = feature.description
        self.component.features.remove(feature)
        return replacement

    def add_protein(self, cds: Part, name: str) -> sbol3.Component:
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
        if any(role_key(role) == role_key(SO_CDS) for role in component.roles):
            continue
        try:
            found.append(Genome(document, component))
        except InputError:
            continue
    return found
