"""Find VSV transcription units in an antigenome string.

VSV transcribes five genes, in order: Nucleoprotein, Phosphoprotein, MatrixProtein,
Glycoprotein, LargePolymerase. Each unit begins with a gene-start (``AACAG``, two variable
bases, ``ATC``) and ends at ``TATG`` followed by seven adenines. Sequences are antigenome sense.
"""

from __future__ import annotations

from dataclasses import dataclass

from vsv.genes import Gene
from vsv.genome.record import GeneSpan, Genome, Span
from vsv.genome.sequence import STOP_CODONS

TRANSCRIPT_END = "TATG"
GENE_END_SIGNAL = TRANSCRIPT_END + "A" * 7
GENE_START_SIGNAL = "AACAGNNATC"
_START_MAX_MISMATCHES = 1
_START_SEARCH_WINDOW = 80


@dataclass(frozen=True, slots=True)
class _Unit:
    """Half-open coordinates of one transcription unit.

    The gene ends after ``TATG``. The adenine run is part of the gene-end signal only.
    """

    start: int
    end_signal: int
    end: int
    orf: tuple[int, int] | None

    @property
    def gene_end(self) -> int:
        return self.end_signal + len(TRANSCRIPT_END)


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
            signal_end = end_position + len(GENE_END_SIGNAL)
            orf = find_orf(sequence, start_position + len(GENE_START_SIGNAL), end_position)
            units.append(_Unit(start_position, end_position, signal_end, orf))
            previous_end = signal_end
    return units


def _span(start: int, end: int) -> Span:
    return Span(start, end)


def annotate(bases: str, *, identifier: str, description: str = "") -> Genome:
    """Leader, genes, signals, coding sequences, and trailer of a VSV antigenome.

    Exactly five units take the gene names in genome order. Any other count is named
    ``gene1``, ``gene2``, and so on. The coding sequence is the first closed ``ATG`` after
    the gene-start, and its product name is ``{name} CDS``.
    """
    units = _find_units(bases)
    names: list[tuple[str, Gene | None]] = (
        [(gene.value, gene) for gene in Gene]
        if len(units) == len(Gene)
        else [(f"gene{index}", None) for index in range(1, len(units) + 1)]
    )
    genes: list[GeneSpan] = []
    for (name, gene), unit in zip(names, units, strict=True):
        cds = None if unit.orf is None else _span(unit.orf[0], unit.orf[1])
        genes.append(
            GeneSpan(
                name,
                gene,
                _span(unit.start, unit.gene_end),
                cds,
                f"{name} CDS",
                _span(unit.start, unit.start + len(GENE_START_SIGNAL)),
                _span(unit.end_signal, unit.end),
            )
        )
    leader = _span(0, units[0].start) if units and units[0].start > 0 else None
    trailer = _span(units[-1].end, len(bases)) if units and units[-1].end < len(bases) else None
    return Genome(identifier, description, bases, leader, trailer, tuple(genes))
