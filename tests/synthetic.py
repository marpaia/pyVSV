"""A synthetic VSV-like genome whose feature coordinates are known exactly."""

from dataclasses import dataclass

from vsv.genes import Gene
from vsv.genome.annotation import GENE_END_SIGNAL, TRANSCRIPT_END

_LEADER = "ACGCTTCCGCTCCTCCTCGCCTCCCTCTCCGCTCCTCCGC"
_TRAILER = "GCGGAGGAGCGGAGAGGGAGGCGAGGAGGAGCGGAAGCGT"
_START = "AACAGTTATC"
_INTERGENIC = "GA"


@dataclass(frozen=True)
class SyntheticGenome:
    """A five-gene genome and the zero-based, half-open spans of everything in it."""

    identifier: str
    description: str
    sequence: str
    leader: tuple[int, int]
    trailer: tuple[int, int]
    genes: dict[str, tuple[int, int]]
    cds: dict[str, tuple[int, int]]
    proteins: dict[str, str]

    @property
    def fasta(self) -> str:
        return f">{self.identifier} {self.description}\n{self.sequence}\n"


def build_genome(
    identifier: str = "SYN0001.1", codons: dict[str, str] | None = None
) -> SyntheticGenome:
    """Assemble a genome from a leader, five gene units, and a trailer.

    Each unit is a gene-start signal, two spacer bases, an open reading frame, and a gene-end
    signal; units are separated by an intergenic dinucleotide.
    """
    codons = codons or {name: "GCT" * (10 + 5 * i) for i, name in enumerate(Gene)}
    sequence = _LEADER
    genes: dict[str, tuple[int, int]] = {}
    cds: dict[str, tuple[int, int]] = {}
    proteins: dict[str, str] = {}
    for i, name in enumerate(Gene):
        if i:
            sequence += _INTERGENIC
        start = len(sequence)
        orf = "ATG" + codons[name] + "TAA"
        cds_start = start + len(_START) + 2
        cds[name] = (cds_start, cds_start + len(orf))
        proteins[name] = "M" + "A" * (len(codons[name]) // 3)
        sequence += _START + "GG" + orf + "CTCT" + GENE_END_SIGNAL
        genes[name] = (start, len(sequence) - len(GENE_END_SIGNAL) + len(TRANSCRIPT_END))
    trailer = (len(sequence), len(sequence) + len(_TRAILER))
    return SyntheticGenome(
        identifier,
        "synthetic vesicular stomatitis virus, complete genome",
        sequence + _TRAILER,
        (0, len(_LEADER)),
        trailer,
        genes,
        cds,
        proteins,
    )


_PRODUCTS = {
    Gene.NUCLEOPROTEIN: "nucleocapsid protein",
    Gene.PHOSPHOPROTEIN: "phosphoprotein",
    Gene.MATRIX_PROTEIN: "matrix protein",
    Gene.GLYCOPROTEIN: "glycoprotein",
    Gene.LARGE_POLYMERASE: "L polymerase protein",
}


def genbank_text(genome: SyntheticGenome, *, genes: bool = False) -> str:
    """Render a synthetic genome as a GenBank record with free-text product names."""
    accession = genome.identifier.split(".")[0]
    lines = [
        f"LOCUS       {accession:<16}{len(genome.sequence):>7} bp    cRNA    linear   VRL",
        f"DEFINITION  {genome.description}.",
        f"ACCESSION   {accession}",
        f"VERSION     {genome.identifier}",
        "SOURCE      Vesicular stomatitis virus",
        "  ORGANISM  Vesicular stomatitis virus",
        "            Viruses.",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{len(genome.sequence)}",
        '                     /organism="Vesicular stomatitis virus"',
        '                     /mol_type="viral cRNA"',
    ]
    for name in Gene:
        if genes:
            start, end = genome.genes[name]
            lines += [
                f"     gene            {start + 1}..{end}",
                f'                     /gene="{name}"',
            ]
        start, end = genome.cds[name]
        lines += [
            f"     CDS             {start + 1}..{end}",
            f'                     /product="{_PRODUCTS[name]}"',
        ]
    lines.append("ORIGIN")
    sequence = genome.sequence.lower()
    for offset in range(0, len(sequence), 60):
        chunk = sequence[offset : offset + 60]
        groups = " ".join(chunk[i : i + 10] for i in range(0, len(chunk), 10))
        lines.append(f"{offset + 1:>9} {groups}")
    lines.append("//")
    return "\n".join(lines) + "\n"
