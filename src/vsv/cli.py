"""Command-line entry point for VSV."""

import argparse
import sys
import textwrap
from collections.abc import Sequence
from pathlib import Path

from vsv import __version__
from vsv.annotation import kind, label, to_gff3
from vsv.api import prepare
from vsv.architecture import Architecture, Derivation, LocusStatus, Module, ParentStatus
from vsv.compare import GenomeComparison, compare_all
from vsv.evidence import EvidenceLevel, Observation, Phenotype
from vsv.formats import OUTPUT_FORMATS, read_genomes, write_document
from vsv.genes import Gene, Product
from vsv.ingest import DEFAULT_QUERY, INGEST_ROOT, ingest
from vsv.library import LIBRARY_ROOT, load_design, load_library
from vsv.query import DesignRecord, select
from vsv.rules import KNOWLEDGE_GAPS, Host, check
from vsv.sbol import Genome
from vsv.sequence import InputError, composition, wrap_fasta

_EVIDENCE = {
    "none": EvidenceLevel.NONE,
    "in-vitro": EvidenceLevel.IN_VITRO,
    "rodent": EvidenceLevel.RODENT,
    "preclinical": EvidenceLevel.RODENT,
    "nonhuman-primate": EvidenceLevel.NONHUMAN_PRIMATE,
    "human": EvidenceLevel.HUMAN,
}


def _cli_name(text: str) -> str:
    """Show a VSV gene as its one-letter CLI name."""
    for gene in sorted(Gene, key=lambda item: len(item.value), reverse=True):
        if text == gene.value:
            return gene.short
        prefix = f"{gene.value} "
        if text.startswith(prefix):
            return f"{gene.short}{text[len(gene.value) :]}"
    return text


def _cli_text(text: str) -> str:
    for gene in sorted(Gene, key=lambda item: len(item.value), reverse=True):
        text = text.replace(gene.value, gene.short)
    return text


def _parse_gene(value: str) -> Gene:
    match = next((gene for gene in Gene if gene.short == value), None)
    if match is None:
        choices = ", ".join(gene.short for gene in Gene)
        raise argparse.ArgumentTypeError(f"invalid gene {value!r}; expected one of {choices}")
    return match


def _row(cells: Sequence[str], widths: Sequence[int]) -> str:
    return "  ".join(cell.ljust(width) for cell, width in zip(cells, widths, strict=True)).rstrip()


def _widths(header: Sequence[str], rows: Sequence[Sequence[str]]) -> list[int]:
    return [max(len(cell) for cell in column) for column in zip(header, *rows, strict=True)]


def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    widths = _widths(header, rows)
    return "\n".join(_row(row, widths) for row in (header, *rows))


def _indent(text: str, spaces: int = 2) -> str:
    pad = " " * spaces
    return "\n".join(pad + line if line else "" for line in text.splitlines())


def _wrap(text: str, indent: str, width: int = 78) -> list[str]:
    return textwrap.wrap(text, width=width, initial_indent=indent, subsequent_indent=indent)


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
                    _cli_name(label(feature)),
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
                    gene.short,
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
                raise InputError(f"{genome.id}: no {product} annotated for gene {gene.short}.")
            entries.append(wrap_fasta(f"{genome.id} {gene.short} {product}", sequence))
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
            g.gene.short,
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
    add_command("genes", "List the n, p, m, g, and l genes and coding sequences of VSV genomes.")
    extract_parser = add_command(
        "extract", "Write gene, coding-sequence, or protein sequences of VSV genomes as FASTA."
    )
    extract_parser.add_argument(
        "-g",
        "--gene",
        dest="genes",
        action="append",
        type=_parse_gene,
        metavar="{n,p,m,g,l}",
        help="gene to extract (n, p, m, g, or l); repeat for several (default: all five)",
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
    designs = commands.add_parser(
        "designs",
        help="List curated VSV architectures and the evidence attached to them.",
    )
    designs.add_argument("--antigen-class", default=None, help="antigen class slug")
    designs.add_argument("--phenotype", default=None, help="observation phenotype")
    designs.add_argument("--host", default=None, help="observation host")
    designs.add_argument(
        "--evidence",
        default=None,
        help="minimum evidence: in-vitro, rodent, preclinical, nonhuman-primate, or human",
    )
    designs.add_argument("--library", type=Path, default=LIBRARY_ROOT, help="design library")
    check = commands.add_parser(
        "check", help="Report design-rule violations for one curated design."
    )
    check.add_argument("design", type=Path, help="design directory")
    show = commands.add_parser(
        "show",
        help="Show one curated design, its implementations, and its observations.",
    )
    show.add_argument("design", help="design id or directory")
    fetch = commands.add_parser("ingest", help="Fetch new VSV vector papers into the review queue.")
    fetch.add_argument("--dest", type=Path, default=INGEST_ROOT, help="ingest directory")
    fetch.add_argument(
        "--since",
        default=None,
        help="only papers first published on or after this date",
    )
    fetch.add_argument("--query", default=DEFAULT_QUERY, help="Europe PMC query")
    return parser


def _phenotype(value: object) -> Phenotype | None:
    if value is None:
        return None
    try:
        return Phenotype(str(value))
    except ValueError as error:
        raise InputError(f"unknown phenotype {value!r}.") from error


def _host(value: object) -> Host | None:
    if value is None:
        return None
    try:
        return Host(str(value))
    except ValueError as error:
        raise InputError(f"unknown host {value!r}.") from error


def _evidence(value: object) -> EvidenceLevel | None:
    if value is None:
        return None
    try:
        return _EVIDENCE[str(value)]
    except KeyError as error:
        raise InputError(f"unknown evidence {value!r}.") from error


def _design_rows(records: Sequence[DesignRecord]) -> str:
    rows = []
    for record in records:
        card = record.card
        rows.append(
            [
                record.architecture.id,
                card.replication.value,
                ",".join(card.antigen_classes) or "-",
                card.evidence_level.value,
                ",".join(card.hosts) or "-",
                _gap_text(record),
            ]
        )
    if not rows:
        return "no designs match"
    return _table(["id", "replication", "antigens", "evidence", "hosts", "gaps"], rows)


def _gap_text(record: DesignRecord) -> str:
    names = [
        item.rule
        for item in check(record.architecture, record.observations)
        if item.rule in KNOWLEDGE_GAPS
    ]
    return ",".join(names) or "-"


def _parent_text(architecture: Architecture) -> str:
    if architecture.parent_status is ParentStatus.KNOWN and architecture.parent_id:
        return architecture.parent_id
    return architecture.parent_status.value


def _derivation_text(derivation: Derivation) -> str:
    target = _cli_name(derivation.target)
    module = None if derivation.module_id is None else _cli_name(derivation.module_id)
    if derivation.activity == "delete-module":
        text = f"delete {target}"
    elif derivation.activity == "insert-module":
        text = f"insert after {target}"
        if module:
            text += f" with {module}"
    else:
        text = f"replace {target}"
        if module and module != target:
            text += f" with {module}"
    if derivation.replication is not None:
        text += f", {derivation.replication.value}"
    return text


def _fact_lines(pairs: Sequence[tuple[str, str]]) -> list[str]:
    width = max(len(key) for key, _value in pairs)
    return [f"{key:<{width}}  {value}" for key, value in pairs]


def _locus(module: Module) -> str:
    if module.locus_status is LocusStatus.NOT_TRANSCRIBED:
        return "-"
    if module.locus_status is LocusStatus.UNKNOWN or module.locus is None:
        return "unknown"
    return str(module.locus)


def _module_table(modules: Sequence[Module]) -> str:
    rows = [
        [_locus(module), _cli_name(module.id), module.antigen_class or ""] for module in modules
    ]
    return _table(["locus", "module", "class"], rows)


def _observation_lines(observations: Sequence[Observation], citation: str) -> list[str]:
    if not observations:
        return ["  no observations"]
    header = ["phenotype", "host", "assay", "outcome"]
    rows = [
        [item.phenotype.value, item.host.value, item.assay.value, item.outcome.value]
        for item in observations
    ]
    widths = _widths(header, rows)
    lines = ["  " + _row(header, widths)]
    for item, row in zip(observations, rows, strict=True):
        lines.append("  " + _row(row, widths))
        if item.citation != citation:
            lines.append(f"    {item.citation}")
        if item.note:
            lines.extend(_wrap(item.note, "    "))
        lines.append("")
    return lines[:-1]


def _show(record: DesignRecord) -> str:
    architecture = record.architecture
    facts = [
        ("backbone", architecture.backbone),
        ("replication", architecture.replication.value),
        ("parent", _parent_text(architecture)),
    ]
    if architecture.derivation is not None:
        facts.append(("derivation", _derivation_text(architecture.derivation)))
    if record.realization:
        facts.append(("realization", record.realization))
    if record.model_source:
        facts.append(("model", record.model_source))
    lines = [
        architecture.id,
        architecture.name,
        "",
        *_fact_lines(facts),
        "",
        "modules",
        _indent(_module_table(architecture.modules)),
    ]
    known = {item.id for item in record.implementations}
    grouped: dict[str, list[Observation]] = {item_id: [] for item_id in known}
    for observation in record.observations:
        if observation.implementation_id in known:
            grouped[observation.implementation_id].append(observation)
    if not record.implementations:
        lines.extend(["", "implementations", "  none"])
    for implementation in record.implementations:
        lines.extend(
            [
                "",
                f"{implementation.stage.value}  {implementation.id}",
                f"  {implementation.citation}",
            ]
        )
        if implementation.note:
            lines.extend(_wrap(implementation.note, "  "))
        lines.append("")
        lines.extend(_observation_lines(grouped[implementation.id], implementation.citation))
    if record.unmodeled:
        key_width = max(len(field.key) for field in record.unmodeled)
        lines.extend(["", "unmodeled"])
        for field in record.unmodeled:
            lines.append(f"  {field.key:<{key_width}}  {field.value}")
            if field.source:
                lines.extend(_wrap(field.source, " " * (key_width + 4)))
    loose = [
        observation.id
        for observation in record.observations
        if observation.implementation_id not in known
    ]
    if loose:
        lines.extend(["", "observations without an implementation", "  " + ", ".join(loose)])
    gaps = [
        item.message
        for item in check(architecture, record.observations)
        if item.rule in KNOWLEDGE_GAPS
    ]
    if gaps:
        lines.extend(["", "gaps", *[f"  {gap}" for gap in gaps]])
    return "\n".join(lines)


def _resolve_show(identifier: str) -> DesignRecord:
    path = Path(identifier)
    if path.is_dir():
        return load_design(path)
    matches = [record for record in load_library() if record.architecture.id == identifier]
    if len(matches) != 1:
        raise InputError(f"{identifier}: design not found.")
    return matches[0]


def _atlas(args: argparse.Namespace) -> tuple[int, str]:
    if args.command == "designs":
        records = select(
            load_library(args.library),
            antigen_class=args.antigen_class,
            phenotype=_phenotype(args.phenotype),
            host=_host(args.host),
            min_evidence=_evidence(args.evidence),
        )
        return 0, _design_rows(records)
    if args.command == "show":
        return 0, _show(_resolve_show(args.design))
    if args.command == "ingest":
        report = ingest(args.dest, query=args.query, since=args.since)
        return 0, (
            f"papers: {report.papers}\nproposals: {report.proposals}\ncursor: {report.cursor}"
        )
    record = load_design(args.design)
    violations = check(record.architecture, record.observations)
    structural = [item for item in violations if item.rule not in KNOWLEDGE_GAPS]
    gaps = [item for item in violations if item.rule in KNOWLEDGE_GAPS]
    if not structural and not gaps:
        return 0, f"{record.architecture.id}: ok"
    blocks = []
    if structural:
        blocks.append(
            _table(["rule", "message"], [(item.rule, item.message) for item in structural])
        )
    if gaps:
        blocks.append(_table(["gap", "message"], [(item.rule, item.message) for item in gaps]))
    return (1 if structural else 0), "\n".join(blocks)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI and return its exit status."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if args.command in {"designs", "check", "show", "ingest"}:
        try:
            status, output = _atlas(args)
        except InputError as error:
            print(f"vsv: error: {error}", file=sys.stderr)
            return 1
        if output:
            print(output)
        return status
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
                    print(f"vsv: warning: {genome.id}: {_cli_text(problem)}", file=sys.stderr)
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
