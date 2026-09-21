"""Genome annotation for vesicular stomatitis virus (VSV) antigenome sequences.

VSV is a non-segmented negative-strand virus. Its five genes (N, P, M, G, L) are transcribed
as separate mRNAs, each delimited by a conserved gene-end signal (``TATG`` followed by a run of
seven adenines, in antigenome sense) and a conserved gene-start signal (``AACAG``, two variable
bases, ``ATC``). Sequences are expected in antigenome (positive) sense.

Annotations are SBOL3 features on a :class:`vsv.sbol.Genome`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sbol3

from vsv import sbol
from vsv.genes import Gene
from vsv.sbol import Genome
from vsv.sequence import STOP_CODONS

TRANSCRIPT_END = "TATG"
GENE_END_SIGNAL = TRANSCRIPT_END + "A" * 7
GENE_START_SIGNAL = "AACAGNNATC"
_START_MAX_MISMATCHES = 1
_START_SEARCH_WINDOW = 80

_KINDS = (
    (sbol.ROLE_LEADER, "leader"),
    (sbol.ROLE_TRAILER, "trailer"),
    (sbol.ROLE_GENE_START, "gene_start"),
    (sbol.ROLE_GENE_END, "gene_end"),
    (sbol.SO_GENE, "gene"),
    (sbol.SO_CDS, "CDS"),
    (sbol3.SO_MRNA, "mRNA"),
    ("https://identifiers.org/SO:0000418", "signal_peptide"),
    ("https://identifiers.org/SO:0000419", "mature_protein_region"),
    ("https://identifiers.org/SO:0001060", "sequence_variant"),
)


@dataclass(frozen=True, slots=True)
class _Unit:
    """A transcription unit: half-open span of its gene-start signal, gene-end signal, and ORF.

    The gene ends after the ``TATG`` of its gene-end signal; the adenine run that follows templates
    the poly(A) tail and is part of the signal but not of the gene.
    """

    start: int
    end_signal: int
    end: int
    orf: tuple[int, int] | None

    @property
    def gene_end(self) -> int:
        return self.end_signal + len(TRANSCRIPT_END)


def kind(feature: sbol3.SequenceFeature) -> str:
    """A short label for what a feature is, such as ``gene``, ``CDS``, or ``leader``."""
    roles = {sbol.role_key(role) for role in feature.roles}
    for role, name in _KINDS:
        if sbol.role_key(role) in roles:
            return name
    return "region"


def label(feature: sbol3.SequenceFeature) -> str:
    """The name of a feature, falling back to its display ID."""
    return str(feature.name or feature.display_id)


def find_gene_ends(sequence: str) -> list[int]:
    """Start positions of every exact gene-end signal, in increasing order."""
    positions: list[int] = []
    at = sequence.find(GENE_END_SIGNAL)
    while at != -1:
        positions.append(at)
        at = sequence.find(GENE_END_SIGNAL, at + 1)
    return positions


def find_gene_start(sequence: str, begin: int, end: int) -> int | None:
    """Position of the best gene-start signal starting in ``[begin, end)``.

    Up to one mismatch to the consensus is tolerated. The hit with the fewest mismatches wins,
    the earliest among equals. Returns ``None`` if no window matches.
    """
    width = len(GENE_START_SIGNAL)
    best: tuple[int, int] | None = None
    for position in range(begin, min(end, len(sequence) - width + 1)):
        window = sequence[position : position + width]
        mismatches = sum(
            1
            for want, got in zip(GENE_START_SIGNAL, window, strict=True)
            if want != "N" and want != got
        )
        if mismatches <= _START_MAX_MISMATCHES and (best is None or mismatches < best[1]):
            best = (position, mismatches)
            if mismatches == 0:
                break
    return None if best is None else best[0]


def find_orf(sequence: str, begin: int, end: int) -> tuple[int, int] | None:
    """The open reading frame from the first ``ATG`` in ``[begin, end)`` to its stop codon.

    Returns half-open coordinates including the stop codon, or ``None`` when there is no start
    codon or the frame is not closed before ``end``.
    """
    start = sequence.find("ATG", begin, end)
    if start == -1:
        return None
    for position in range(start, end - 2, 3):
        if sequence[position : position + 3] in STOP_CODONS:
            return start, position + 3
    return None


def _find_units(sequence: str) -> list[_Unit]:
    units: list[_Unit] = []
    previous_end = 0
    for end_position in find_gene_ends(sequence):
        window_end = min(previous_end + _START_SEARCH_WINDOW, end_position)
        start_position = find_gene_start(sequence, previous_end, window_end)
        if start_position is not None:
            gene_end = end_position + len(GENE_END_SIGNAL)
            orf = find_orf(sequence, start_position + len(GENE_START_SIGNAL), end_position)
            units.append(_Unit(start_position, end_position, gene_end, orf))
            previous_end = gene_end
    return units


def is_annotated(genome: Genome) -> bool:
    """Whether the genome already has gene or coding-sequence features."""
    return bool(genome.features(sbol.SO_GENE) or genome.features(sbol.SO_CDS))


def annotate(genome: Genome) -> None:
    """Add the leader, genes, signals, coding sequences, and trailer of a VSV genome.

    Each gene-end signal closes a transcription unit whose start signal is searched for shortly
    after the previous unit. Units are named N, P, M, G, L in genome order only when exactly five
    are found; otherwise they are named ``gene1``, ``gene2``, and so on. Every coding sequence is
    also recorded as a protein component linked to it.
    """
    units = _find_units(genome.bases)
    names = (
        [str(gene) for gene in Gene]
        if len(units) == len(Gene)
        else [f"gene{i}" for i in range(1, len(units) + 1)]
    )
    if units:
        genome.add_region("leader", 0, units[0].start, [sbol.SO_REGION, sbol.ROLE_LEADER])
    for name, unit in zip(names, units, strict=True):
        genome.add_region(name, unit.start, unit.gene_end, [sbol.SO_GENE])
        genome.add_region(
            f"{name} gene start",
            unit.start,
            unit.start + len(GENE_START_SIGNAL),
            [sbol.SO_REGULATORY_REGION, sbol.ROLE_GENE_START],
        )
        genome.add_region(
            f"{name} gene end",
            unit.end_signal,
            unit.end,
            [sbol.SO_TERMINATOR, sbol.ROLE_GENE_END],
        )
        if unit.orf is not None:
            cds = genome.add_region(name, unit.orf[0], unit.orf[1], [sbol.SO_CDS])
            genome.add_protein(cds, name)
    if units and units[-1].end < len(genome):
        genome.add_region(
            "trailer", units[-1].end, len(genome), [sbol.SO_REGION, sbol.ROLE_TRAILER]
        )


def name_curated_genes(genome: Genome) -> None:
    """Name the genes and coding sequences of a curated genome N, P, M, G, L by position.

    Curated records label features with free-text products such as ``phosphoprotein`` or
    ``NS protein``, and often leave genes unnamed. VSV genes always occur in the order N, P, M, G,
    L, so five genes are named in that order, and each takes the longest coding sequence within
    it. Without gene features, exactly five coding sequences are named in order instead. Genomes
    that fit neither pattern are left as they are. A name that is replaced is kept as the
    feature description.
    """
    genes = genome.features(sbol.SO_GENE)
    cdss = genome.features(sbol.SO_CDS)
    if len(genes) == len(Gene):
        for name, gene in zip(Gene, genes, strict=True):
            start, end = genome.span(gene)
            inside = [c for c in cdss if start <= genome.span(c)[0] and genome.span(c)[1] <= end]
            _rename(gene, name)
            if inside:
                _rename(max(inside, key=lambda c: genome.span(c)[1] - genome.span(c)[0]), name)
    elif not genes and len(cdss) == len(Gene):
        for name, cds in zip(Gene, cdss, strict=True):
            _rename(cds, name)


def _rename(feature: sbol3.SequenceFeature, name: str) -> None:
    if feature.name == name:
        return
    if feature.name and not feature.description:
        feature.description = feature.name
    feature.name = name


def to_gff3(genomes: Sequence[Genome]) -> str:
    """Render genome annotations as GFF3 with one-based, inclusive coordinates."""
    lines = ["##gff-version 3"]
    for genome in genomes:
        lines.append(f"##sequence-region {genome.id} 1 {len(genome)}")
    for genome in genomes:
        for feature in genome.features():
            start, end = genome.span(feature)
            feature_kind = kind(feature)
            attributes = f"ID={genome.id}:{feature_kind}:{label(feature)};Name={label(feature)}"
            lines.append(
                "\t".join(
                    [
                        genome.id,
                        "vsv",
                        feature_kind,
                        str(start + 1),
                        str(end),
                        ".",
                        "-" if genome.is_reverse(feature) else "+",
                        "0" if feature_kind == "CDS" else ".",
                        attributes,
                    ]
                )
            )
    return "\n".join(lines) + "\n"
