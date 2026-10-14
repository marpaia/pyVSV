"""A VSV genome as a sequence plus the spans of its genes.

Coordinates are zero-based and half-open. The antigenome is positive sense.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from vsv.errors import InputError
from vsv.genes import Gene, Product
from vsv.genome.sbol.write import write as _write_document
from vsv.genome.sequence import IUPAC_NUCLEOTIDES, translate


@dataclass(frozen=True, slots=True)
class Span:
    """A half-open interval on the antigenome."""

    start: int
    end: int
    strand: str = "+"

    def __post_init__(self) -> None:
        if self.strand not in {"+", "-"}:
            raise InputError(f"Unknown strand {self.strand!r}.")
        if self.start < 0 or self.end < self.start:
            raise InputError(f"Invalid span {self.start}..{self.end}.")

    def __len__(self) -> int:
        return self.end - self.start


@dataclass(frozen=True, slots=True)
class GeneSpan:
    """One transcription unit and the coding sequence inside it.

    ``unit`` ends after the ``TATG`` of the gene-end signal. The poly(A) run that follows
    is ``gene_end`` and is not part of the gene. ``product`` is the name of the protein.
    """

    name: str
    gene: Gene | None
    unit: Span
    cds: Span | None
    product: str
    gene_start: Span | None = None
    gene_end: Span | None = None

    def sequence(self, bases: str, product: Product) -> str | None:
        """The transcription unit, its coding sequence, or the protein it encodes."""
        if product is Product.GENE:
            return bases[self.unit.start : self.unit.end]
        if self.cds is None:
            return None
        coding = bases[self.cds.start : self.cds.end]
        if product is Product.CDS:
            return coding
        return translate(coding, to_stop=True)


@dataclass(frozen=True, slots=True)
class Feature:
    """One annotated interval used by the annotation table."""

    kind: str
    name: str
    span: Span


@dataclass(frozen=True, slots=True)
class Genome:
    """An antigenome sequence and the genes identified on it."""

    id: str
    description: str
    bases: str
    leader: Span | None = None
    trailer: Span | None = None
    genes: tuple[GeneSpan, ...] = ()

    def __post_init__(self) -> None:
        if not self.id:
            raise InputError("A genome requires an identifier.")
        invalid = set(self.bases) - IUPAC_NUCLEOTIDES
        if not self.bases or invalid:
            symbols = "".join(sorted(invalid)) or "empty sequence"
            raise InputError(f"{self.id}: sequence is not IUPAC DNA ({symbols}).")
        length = len(self.bases)
        for span in self._spans():
            if span.end > length:
                raise InputError(f"{self.id}: a feature ends past the sequence.")
        object.__setattr__(self, "genes", tuple(self.genes))

    def __len__(self) -> int:
        return len(self.bases)

    def _spans(self) -> list[Span]:
        spans: list[Span] = []
        if self.leader is not None:
            spans.append(self.leader)
        if self.trailer is not None:
            spans.append(self.trailer)
        for gene in self.genes:
            spans.append(gene.unit)
            if gene.cds is not None:
                spans.append(gene.cds)
            if gene.gene_start is not None:
                spans.append(gene.gene_start)
            if gene.gene_end is not None:
                spans.append(gene.gene_end)
        return spans

    @property
    def named(self) -> dict[Gene, GeneSpan]:
        """The five VSV genes that were identified, in genome order."""
        return {gene.gene: gene for gene in self.genes if gene.gene is not None}

    def gene(self, name: Gene) -> GeneSpan:
        """The transcription unit of ``name``."""
        found = self.named.get(name)
        if found is None:
            raise InputError(f"{self.id}: no {name.value} annotated.")
        return found

    def cds(self, name: Gene) -> str:
        """The coding sequence of ``name``, including its stop codon."""
        coding = self.gene(name).sequence(self.bases, Product.CDS)
        if coding is None:
            raise InputError(f"{self.id}: no coding sequence annotated for {name.value}.")
        return coding

    def protein(self, name: Gene) -> str | None:
        """The translation of ``name`` without its stop codon."""
        found = self.named.get(name)
        if found is None:
            return None
        return found.sequence(self.bases, Product.PROTEIN)

    def features(self) -> tuple[Feature, ...]:
        """Leader, each gene and its signals and coding sequence, then the trailer."""
        rows: list[Feature] = []
        if self.leader is not None:
            rows.append(Feature("leader", "leader", self.leader))
        for gene in self.genes:
            rows.append(Feature("gene", gene.name, gene.unit))
            if gene.gene_start is not None:
                rows.append(Feature("gene_start", f"{gene.name} gene start", gene.gene_start))
            if gene.gene_end is not None:
                rows.append(Feature("gene_end", f"{gene.name} gene end", gene.gene_end))
            if gene.cds is not None:
                rows.append(Feature("CDS", gene.product, gene.cds))
        if self.trailer is not None:
            rows.append(Feature("trailer", "trailer", self.trailer))
        return tuple(rows)

    def write(self, path: str | Path) -> None:
        """Write this genome as an SBOL3 document. Callers do not handle SBOL objects."""
        _write_document(self, path)
