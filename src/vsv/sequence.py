"""Nucleotide and protein sequence records, FASTA input and output, and basic statistics."""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

IUPAC_NUCLEOTIDES = frozenset("ACGTRYSWKMBDHVN")
_UNAMBIGUOUS = "ACGT"
_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)

_BASES = "TCAG"
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE: dict[str, str] = {
    a + b + c: _AMINO_ACIDS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}

STOP_CODONS = frozenset(codon for codon, aa in CODON_TABLE.items() if aa == "*")


class InputError(ValueError):
    """Raised when an input file cannot be read as genome data."""


class FastaError(InputError):
    """Raised when FASTA input is malformed."""


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    """A named nucleotide sequence.

    ``sequence`` is upper-case IUPAC nucleotide text using ``T`` for thymine.
    """

    id: str
    description: str
    sequence: str

    def __post_init__(self) -> None:
        if not self.id:
            raise FastaError("A sequence record requires a non-empty identifier.")
        invalid = set(self.sequence) - IUPAC_NUCLEOTIDES
        if invalid:
            raise FastaError(
                f"Record {self.id} contains non-nucleotide symbols: {''.join(sorted(invalid))}"
            )

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def header(self) -> str:
        """The FASTA header line without the leading ``>``."""
        return f"{self.id} {self.description}".rstrip()


@dataclass(frozen=True, slots=True)
class Composition:
    """Nucleotide composition of a sequence."""

    length: int
    counts: dict[str, int]
    ambiguous: int

    @property
    def unambiguous(self) -> int:
        return self.length - self.ambiguous

    @property
    def gc_content(self) -> float:
        """Fraction of unambiguous bases that are G or C; ``nan`` if there are none."""
        if self.unambiguous == 0:
            return float("nan")
        return (self.counts["G"] + self.counts["C"]) / self.unambiguous


def parse_fasta(lines: Iterable[str]) -> Iterator[SequenceRecord]:
    """Yield records from FASTA text lines.

    Sequence text is upper-cased, ``U`` is read as ``T``, and whitespace is ignored.
    """
    header: str | None = None
    chunks: list[str] = []
    for number, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if header is not None:
                yield _build_record(header, chunks)
            header, chunks = line[1:].strip(), []
        elif header is None:
            raise FastaError(f"Line {number}: sequence data before the first '>' header.")
        else:
            chunks.append(re.sub(r"\s+", "", line))
    if header is not None:
        yield _build_record(header, chunks)


def _build_record(header: str, chunks: list[str]) -> SequenceRecord:
    identifier, _, description = header.partition(" ")
    sequence = "".join(chunks).upper().replace("U", "T")
    if not sequence:
        raise FastaError(f"Record {identifier or '(unnamed)'} has no sequence.")
    return SequenceRecord(identifier, description.strip(), sequence)


def read_fasta(path: str | Path) -> list[SequenceRecord]:
    """Read every record from a FASTA file."""
    path = Path(path)
    try:
        with path.open(encoding="ascii") as handle:
            return list(parse_fasta(handle))
    except FastaError as error:
        raise FastaError(f"{path}: {error}") from error
    except UnicodeDecodeError as error:
        raise FastaError(f"{path}: file is not ASCII text.") from error


def wrap_fasta(header: str, sequence: str, width: int = 70) -> str:
    """Render one FASTA entry, of nucleotides or protein, with lines wrapped at ``width``."""
    if width < 1:
        raise ValueError("width must be positive.")
    lines = [f">{header}", *(sequence[i : i + width] for i in range(0, len(sequence), width))]
    return "\n".join(lines) + "\n"


def format_fasta(records: Iterable[SequenceRecord], width: int = 70) -> str:
    """Render records as FASTA text with sequence lines wrapped at ``width``."""
    return "".join(wrap_fasta(record.header, record.sequence, width) for record in records)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an IUPAC nucleotide string."""
    return sequence.translate(_COMPLEMENT)[::-1]


def translate(sequence: str, *, to_stop: bool = False) -> str:
    """Translate codons in frame from the first base with the standard genetic code.

    Codons containing ambiguous bases translate to ``X``. A trailing partial codon is ignored.
    With ``to_stop`` the protein ends before the first stop codon; otherwise stops are ``*``.
    """
    protein: list[str] = []
    for i in range(0, len(sequence) - len(sequence) % 3, 3):
        amino_acid = CODON_TABLE.get(sequence[i : i + 3], "X")
        if amino_acid == "*" and to_stop:
            break
        protein.append(amino_acid)
    return "".join(protein)


def composition(sequence: str) -> Composition:
    """Count nucleotides in a sequence."""
    tally = Counter(sequence)
    counts = {base: tally[base] for base in _UNAMBIGUOUS}
    ambiguous = len(sequence) - sum(counts.values())
    return Composition(length=len(sequence), counts=counts, ambiguous=ambiguous)
