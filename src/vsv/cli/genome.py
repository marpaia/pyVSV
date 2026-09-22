"""Commands that read FASTA, GenBank, and SBOL3 genomes."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import click
import sbol3

from vsv.cli.render import cli_name, cli_text, table
from vsv.errors import InputError
from vsv.genes import Gene, Product
from vsv.genome import prepare
from vsv.genome.annotation import kind, label, to_gff3
from vsv.genome.compare import GenomeComparison, compare_all
from vsv.genome.formats import OUTPUT_FORMATS, read_genomes, write_document
from vsv.genome.sbol import Genome
from vsv.genome.sequence import composition, wrap_fasta

_FORMATS = ["table", "gff3", *OUTPUT_FORMATS]


class GeneName(click.ParamType[Gene]):
    """A VSV gene named ``n``, ``p``, ``m``, ``g``, or ``l``."""

    name = "{n,p,m,g,l}"

    def convert(
        self, value: object, param: click.Parameter | None, ctx: click.Context | None
    ) -> Gene:
        if isinstance(value, Gene):
            return value
        text = str(value)
        match = next((gene for gene in Gene if gene.short == text), None)
        if match is None:
            choices = ", ".join(gene.short for gene in Gene)
            self.fail(f"invalid gene {text!r}; expected one of {choices}", param, ctx)
        return match


class ProductName(click.ParamType[Product]):
    """``gene``, ``cds``, or ``protein``."""

    name = "product"

    def convert(
        self, value: object, param: click.Parameter | None, ctx: click.Context | None
    ) -> Product:
        if isinstance(value, Product):
            return value
        text = str(value)
        try:
            return Product(text)
        except ValueError:
            choices = ", ".join(product.value for product in Product)
            self.fail(f"invalid product {text!r}; expected one of {choices}", param, ctx)


def read(
    files: Sequence[Path], *, minimum: int, identify: bool, command: str
) -> tuple[sbol3.Document, list[Genome]]:
    """Read ``files`` and, when ``identify`` is set, name their VSV genes."""
    document, genomes = read_genomes(list(files))
    if len(genomes) < minimum:
        raise InputError(f"{command} requires at least {minimum} records.")
    if identify:
        for genome in genomes:
            for problem in prepare(genome):
                click.echo(f"vsv: warning: {genome.id}: {cli_text(problem)}", err=True)
    return document, genomes


def _summary(genomes: Sequence[Genome]) -> str:
    rows = []
    for genome in genomes:
        stats = composition(genome.bases)
        rows.append(
            [
                genome.id,
                str(stats.length),
                *(str(stats.counts[base]) for base in "ACGT"),
                str(stats.ambiguous),
                f"{stats.gc_content:.2%}",
                genome.description,
            ]
        )
    header = ["id", "length", "A", "C", "G", "T", "ambiguous", "GC", "description"]
    return table(header, rows)


def _annotation_table(genomes: Sequence[Genome]) -> str:
    rows = []
    for genome in genomes:
        for feature in genome.features():
            start, end = genome.span(feature)
            feature_kind = kind(feature)
            rows.append(
                [
                    genome.id,
                    feature_kind,
                    cli_name(label(feature)),
                    str(start + 1),
                    str(end),
                    str(end - start),
                    str(len(genome.translate(feature))) if feature_kind == "CDS" else "",
                ]
            )
    return table(["id", "feature", "name", "start", "end", "length", "aa"], rows)


def _span_cells(span: tuple[int, int] | None) -> list[str]:
    if span is None:
        return ["", "", ""]
    start, end = span
    return [str(start + 1), str(end), str(end - start)]


def _gene_table(genomes: Sequence[Genome]) -> str:
    rows = []
    for genome in genomes:
        found = genome.genes
        for gene in Gene:
            annotation = found.get(gene)
            protein = None if annotation is None else annotation.protein
            rows.append(
                [
                    genome.id,
                    gene.short,
                    *_span_cells(None if annotation is None else annotation.gene_span),
                    *_span_cells(None if annotation is None else annotation.cds_span),
                    "" if protein is None else str(len(protein)),
                ]
            )
    header = ["id", "gene", "start", "end", "length", "cds_start", "cds_end", "cds_length", "aa"]
    return table(header, rows)


def _extract(genomes: Sequence[Genome], genes: Sequence[Gene], product: Product) -> str:
    entries = []
    for genome in genomes:
        found = genome.genes
        for gene in genes:
            annotation = found.get(gene)
            sequence = None if annotation is None else annotation.sequence(product)
            if sequence is None:
                raise InputError(f"{genome.id}: no {product} annotated for gene {gene.short}.")
            entries.append(wrap_fasta(f"{genome.id} {gene.short} {product}", sequence))
    return "".join(entries).rstrip("\n")


def _comparison_report(comparisons: Sequence[GenomeComparison], k: int) -> str:
    genome_rows = []
    for comparison in comparisons:
        distance = comparison.kmers.mash_distance
        genome_rows.append(
            [
                comparison.id_a,
                comparison.id_b,
                str(comparison.length_a),
                str(comparison.length_b),
                f"{comparison.kmers.jaccard:.4f}",
                "NA" if distance is None else f"{distance:.4f}",
            ]
        )
    gene_rows = [
        [
            comparison.id_a,
            comparison.id_b,
            gene.gene.short,
            str(gene.alignment.length_a),
            str(gene.alignment.length_b),
            f"{gene.alignment.identity:.2%}",
        ]
        for comparison in comparisons
        for gene in comparison.genes
    ]
    return "\n".join(
        [
            f"Genome k-mer similarity (k={k})",
            table(["a", "b", "len_a", "len_b", "jaccard", "mash_distance"], genome_rows),
            "",
            "Protein identity by gene",
            table(["a", "b", "gene", "aa_a", "aa_b", "identity"], gene_rows),
        ]
    )


@click.command()
@click.argument("files", nargs=-1, required=True, type=click.Path(path_type=Path))
def summary(files: tuple[Path, ...]) -> None:
    """Report length, base composition, and GC content of each record."""
    _document, genomes = read(files, minimum=1, identify=False, command="summary")
    click.echo(_summary(genomes))


@click.command()
@click.argument("files", nargs=-1, required=True, type=click.Path(path_type=Path))
@click.option(
    "--format",
    "output_format",
    type=click.Choice(_FORMATS),
    default="table",
    show_default=True,
    help="A table, GFF3, or SBOL3 serialized as turtle, ntriples, rdfxml, or jsonld.",
)
def annotate(files: tuple[Path, ...], output_format: str) -> None:
    """Annotate leader, genes, signals, coding sequences, and trailer of VSV genomes.

    Genomes that already carry annotations keep them.
    """
    document, genomes = read(files, minimum=1, identify=True, command="annotate")
    if output_format == "gff3":
        click.echo(to_gff3(genomes).rstrip("\n"))
    elif output_format == "table":
        click.echo(_annotation_table(genomes))
    else:
        click.echo(write_document(document, output_format).rstrip("\n"))


@click.command()
@click.argument("files", nargs=-1, required=True, type=click.Path(path_type=Path))
def genes(files: tuple[Path, ...]) -> None:
    """List the n, p, m, g, and l genes and coding sequences of VSV genomes."""
    _document, genomes = read(files, minimum=1, identify=True, command="genes")
    click.echo(_gene_table(genomes))


@click.command()
@click.argument("files", nargs=-1, required=True, type=click.Path(path_type=Path))
@click.option(
    "-g",
    "--gene",
    "genes",
    multiple=True,
    type=GeneName(),
    metavar="{n,p,m,g,l}",
    help="Gene to extract. Repeat for several. The default is all five.",
)
@click.option(
    "--kind",
    type=ProductName(),
    default=Product.CDS,
    show_default=True,
    help="The transcription unit, its coding sequence, or the protein.",
)
def extract(files: tuple[Path, ...], genes: tuple[Gene, ...], kind: Product) -> None:
    """Write gene, coding-sequence, or protein sequences of VSV genomes as FASTA."""
    _document, genomes = read(files, minimum=1, identify=True, command="extract")
    click.echo(_extract(genomes, genes or tuple(Gene), kind))


@click.command()
@click.argument("files", nargs=-1, required=True, type=click.Path(path_type=Path))
@click.option(
    "-k",
    type=click.IntRange(min=1),
    default=12,
    show_default=True,
    help="k-mer size.",
)
def compare(files: tuple[Path, ...], k: int) -> None:
    """Compare genomes by k-mer similarity and per-gene protein identity."""
    _document, genomes = read(files, minimum=2, identify=True, command="compare")
    click.echo(_comparison_report(compare_all(genomes, k), k))
