"""Commands that read VSV genomes and write sequences as SBOL."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import click

from vsv.cli.render import cli_name, cli_text, table
from vsv.errors import InputError
from vsv.genes import Gene, Product
from vsv.genome import gene_problems
from vsv.genome.formats import read_records
from vsv.genome.record import Genome
from vsv.genome.sbol.write import write_sequences
from vsv.genome.sequence import translate
from vsv.sbol import DNA_ENCODING, PROTEIN_ENCODING


class GeneName(click.ParamType[Gene]):
    """A VSV gene: ``n`` nucleoprotein, ``p`` phosphoprotein, ``m`` matrix protein,
    ``g`` glycoprotein, or ``l`` large polymerase.
    """

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


def _read(files: Sequence[Path], *, minimum: int, identify: bool, command: str) -> list[Genome]:
    """Read ``files`` and, when ``identify`` is set, warn about unnamed VSV genes."""
    genomes = read_records(list(files))
    if len(genomes) < minimum:
        raise InputError(f"{command} requires at least {minimum} records.")
    if identify:
        for genome in genomes:
            for problem in gene_problems(genome):
                click.echo(f"vsv: warning: {genome.id}: {cli_text(problem)}", err=True)
    return genomes


def _annotation_table(genomes: Sequence[Genome]) -> str:
    rows = []
    for genome in genomes:
        for feature in genome.features():
            start, end = feature.span.start, feature.span.end
            amino = ""
            if feature.kind == "CDS":
                protein = translate(genome.bases[start:end], to_stop=True)
                amino = str(len(protein))
            rows.append(
                [
                    genome.id,
                    feature.kind,
                    cli_name(feature.name),
                    str(start + 1),
                    str(end),
                    str(end - start),
                    amino,
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
        found = genome.named
        for gene in Gene:
            annotation = found.get(gene)
            protein = (
                None if annotation is None else annotation.sequence(genome.bases, Product.PROTEIN)
            )
            rows.append(
                [
                    genome.id,
                    gene.short,
                    *_span_cells(
                        None if annotation is None else (annotation.unit.start, annotation.unit.end)
                    ),
                    *_span_cells(
                        None
                        if annotation is None or annotation.cds is None
                        else (annotation.cds.start, annotation.cds.end)
                    ),
                    "" if protein is None else str(len(protein)),
                ]
            )
    header = ["id", "gene", "start", "end", "length", "cds_start", "cds_end", "cds_length", "aa"]
    return table(header, rows)


def _extract(genomes: Sequence[Genome], genes: Sequence[Gene], product: Product) -> str:
    encoding = PROTEIN_ENCODING if product is Product.PROTEIN else DNA_ENCODING
    entries = []
    for genome in genomes:
        for gene in genes:
            annotation = genome.named.get(gene)
            sequence = None if annotation is None else annotation.sequence(genome.bases, product)
            if sequence is None:
                raise InputError(
                    f"{genome.id}: no {product.value} annotated for gene {gene.short}."
                )
            entries.append((f"{genome.id} {gene.short} {product.value}", sequence, encoding))
    return write_sequences(entries).rstrip("\n")


@click.command()
@click.argument("files", nargs=-1, required=True, type=click.Path(path_type=Path))
def annotate(files: tuple[Path, ...]) -> None:
    """Annotate leader, genes, signals, coding sequences, and trailer of VSV genomes."""
    click.echo(_annotation_table(_read(files, minimum=1, identify=True, command="annotate")))


@click.command()
@click.argument("files", nargs=-1, required=True, type=click.Path(path_type=Path))
def genes(files: tuple[Path, ...]) -> None:
    """List each VSV gene and its coding sequence.

    Nucleoprotein, phosphoprotein, matrix protein, glycoprotein, and large polymerase.
    """
    click.echo(_gene_table(_read(files, minimum=1, identify=True, command="genes")))


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
    """Write gene, coding-sequence, or protein sequences of VSV genomes as SBOL."""
    genomes = _read(files, minimum=1, identify=True, command="extract")
    click.echo(_extract(genomes, genes or tuple(Gene), kind))
