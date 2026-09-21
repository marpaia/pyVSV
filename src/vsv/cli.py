"""Command-line entry point for VSV."""

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from vsv import __version__
from vsv.annotation import kind, label, to_gff3
from vsv.api import prepare
from vsv.compare import GenomeComparison, compare_all
from vsv.formats import OUTPUT_FORMATS, read_genomes, write_document
from vsv.genes import Gene, Product
from vsv.sbol import Genome
from vsv.sequence import InputError, composition, wrap_fasta


def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    widths = [max(len(cell) for cell in column) for column in zip(header, *rows, strict=True)]
    return "\n".join(
        "  ".join(cell.ljust(width) for cell, width in zip(row, widths, strict=True)).rstrip()
        for row in (header, *rows)
    )


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
    return _table(header, rows)


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
                    label(feature),
                    str(start + 1),
                    str(end),
                    str(end - start),
                    str(len(genome.translate(feature))) if feature_kind == "CDS" else "",
                ]
            )
    return _table(["id", "feature", "name", "start", "end", "length", "aa"], rows)


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
                    gene,
                    *_span_cells(None if annotation is None else annotation.gene_span),
                    *_span_cells(None if annotation is None else annotation.cds_span),
                    "" if protein is None else str(len(protein)),
                ]
            )
    header = ["id", "gene", "start", "end", "length", "cds_start", "cds_end", "cds_length", "aa"]
    return _table(header, rows)


def _extract(genomes: Sequence[Genome], genes: Sequence[Gene], product: Product) -> str:
    entries = []
    for genome in genomes:
        found = genome.genes
        for gene in genes:
            annotation = found.get(gene)
            sequence = None if annotation is None else annotation.sequence(product)
            if sequence is None:
                raise InputError(f"{genome.id}: no {product} annotated for gene {gene}.")
            entries.append(wrap_fasta(f"{genome.id} {gene} {product}", sequence))
    return "".join(entries).rstrip("\n")


def _comparison_report(comparisons: Sequence[GenomeComparison], k: int) -> str:
    genome_rows = []
    for c in comparisons:
        distance = c.kmers.mash_distance
        genome_rows.append(
            [
                c.id_a,
                c.id_b,
                str(c.length_a),
                str(c.length_b),
                f"{c.kmers.jaccard:.4f}",
                "NA" if distance is None else f"{distance:.4f}",
            ]
        )
    gene_rows = [
        [
            c.id_a,
            c.id_b,
            g.gene,
            str(g.alignment.length_a),
            str(g.alignment.length_b),
            f"{g.alignment.identity:.2%}",
        ]
        for c in comparisons
        for g in c.genes
    ]
    return "\n".join(
        [
            f"Genome k-mer similarity (k={k})",
            _table(["a", "b", "len_a", "len_b", "jaccard", "mash_distance"], genome_rows),
            "",
            "Protein identity by gene",
            _table(["a", "b", "gene", "aa_a", "aa_b", "identity"], gene_rows),
        ]
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vsv",
        description="Computational genetics tools for vesicular stomatitis virus research.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", metavar="command")

    def add_command(name: str, help_text: str, *, minimum: int = 1) -> argparse.ArgumentParser:
        command = commands.add_parser(name, help=help_text, description=help_text)
        command.add_argument(
            "inputs",
            nargs="+",
            type=Path,
            metavar="file",
            help="FASTA, GenBank, or SBOL3 files to read",
        )
        command.set_defaults(minimum=minimum)
        return command

    add_command("summary", "Report length, base composition, and GC content of each record.")
    annotate_parser = add_command(
        "annotate",
        "Annotate leader, genes, signals, coding sequences, and trailer of VSV genomes. "
        "Genomes that already carry annotations keep them.",
    )
    annotate_parser.add_argument(
        "--format",
        choices=["table", "gff3", *OUTPUT_FORMATS],
        default="table",
        help="output format: a table, GFF3, or SBOL3 serialized as turtle, ntriples, rdfxml, "
        "or jsonld (default: table)",
    )
    add_command("genes", "List the N, P, M, G, and L genes and coding sequences of VSV genomes.")
    extract_parser = add_command(
        "extract", "Write gene, coding-sequence, or protein sequences of VSV genomes as FASTA."
    )
    extract_parser.add_argument(
        "-g",
        "--gene",
        dest="genes",
        action="append",
        type=Gene,
        choices=list(Gene),
        help="gene to extract; repeat for several (default: all five)",
    )
    extract_parser.add_argument(
        "--kind",
        type=Product,
        choices=list(Product),
        default=Product.CDS,
        help="the transcription unit, its coding sequence, or the protein (default: cds)",
    )
    compare_parser = add_command(
        "compare",
        "Compare genomes by k-mer similarity and per-gene protein identity.",
        minimum=2,
    )
    compare_parser.add_argument("-k", type=int, default=12, help="k-mer size (default: 12)")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI and return its exit status."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if args.command == "compare" and args.k < 1:
        parser.error("-k must be positive.")
    try:
        document, genomes = read_genomes(args.inputs)
        if len(genomes) < args.minimum:
            raise InputError(f"{args.command} requires at least {args.minimum} records.")
        if args.command == "summary":
            output = _summary(genomes)
        else:
            for genome in genomes:
                for problem in prepare(genome):
                    print(f"vsv: warning: {genome.id}: {problem}", file=sys.stderr)
            if args.command == "compare":
                output = _comparison_report(compare_all(genomes, args.k), args.k)
            elif args.command == "genes":
                output = _gene_table(genomes)
            elif args.command == "extract":
                output = _extract(genomes, args.genes or list(Gene), args.kind)
            elif args.format == "gff3":
                output = to_gff3(genomes).rstrip("\n")
            elif args.format == "table":
                output = _annotation_table(genomes)
            else:
                output = write_document(document, args.format).rstrip("\n")
    except InputError as error:
        print(f"vsv: error: {error}", file=sys.stderr)
        return 1
    print(output)
    return 0
