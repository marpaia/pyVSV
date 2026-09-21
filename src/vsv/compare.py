"""Comparative measures between sequences: k-mer similarity and protein alignment identity."""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations

from vsv.genes import Gene
from vsv.sbol import Genome

_MATCH = 2
_MISMATCH = -1
_GAP = -2


@dataclass(frozen=True, slots=True)
class KmerSimilarity:
    """Similarity of two sequences by their sets of distinct k-mers."""

    k: int
    shared: int
    union: int

    @property
    def jaccard(self) -> float:
        return self.shared / self.union if self.union else 0.0

    @property
    def mash_distance(self) -> float | None:
        """Estimated substitutions per site (Mash distance); ``None`` with no shared k-mers."""
        if self.shared == 0:
            return None
        jaccard = self.jaccard
        return -math.log(2 * jaccard / (1 + jaccard)) / self.k


@dataclass(frozen=True, slots=True)
class ProteinAlignment:
    """A global pairwise protein alignment."""

    matches: int
    columns: int
    length_a: int
    length_b: int

    @property
    def identity(self) -> float:
        """Identical columns over alignment columns."""
        return self.matches / self.columns if self.columns else 0.0


@dataclass(frozen=True, slots=True)
class GeneComparison:
    """Protein-level comparison of one gene between two annotated genomes."""

    gene: Gene
    alignment: ProteinAlignment


@dataclass(frozen=True, slots=True)
class GenomeComparison:
    """Comparison of two annotated genomes."""

    id_a: str
    id_b: str
    length_a: int
    length_b: int
    kmers: KmerSimilarity
    genes: tuple[GeneComparison, ...]


def kmer_set(sequence: str, k: int) -> set[str]:
    """Distinct k-mers of ``sequence`` that contain only A, C, G, and T."""
    if k < 1:
        raise ValueError("k must be positive.")
    return {
        sequence[i : i + k]
        for i in range(len(sequence) - k + 1)
        if not set(sequence[i : i + k]) - set("ACGT")
    }


def kmer_similarity(a: str, b: str, k: int = 12) -> KmerSimilarity:
    """Compare the distinct k-mer sets of two sequences."""
    set_a, set_b = kmer_set(a, k), kmer_set(b, k)
    shared = len(set_a & set_b)
    return KmerSimilarity(k=k, shared=shared, union=len(set_a) + len(set_b) - shared)


def align_proteins(a: str, b: str) -> ProteinAlignment:
    """Globally align two protein sequences (match +2, mismatch -1, gap -2) and count identity."""
    if not a or not b:
        return ProteinAlignment(0, max(len(a), len(b)), len(a), len(b))
    rows, cols = len(a) + 1, len(b) + 1
    # Traceback directions: 0 diagonal, 1 up (gap in b), 2 left (gap in a).
    trace = [bytearray(cols) for _ in range(rows)]
    for j in range(1, cols):
        trace[0][j] = 2
    for i in range(1, rows):
        trace[i][0] = 1
    previous = [j * _GAP for j in range(cols)]
    for i in range(1, rows):
        residue = a[i - 1]
        current = [i * _GAP] + [0] * (cols - 1)
        row = trace[i]
        for j in range(1, cols):
            diagonal = previous[j - 1] + (_MATCH if residue == b[j - 1] else _MISMATCH)
            up = previous[j] + _GAP
            left = current[j - 1] + _GAP
            if diagonal >= up and diagonal >= left:
                current[j] = diagonal
            elif up >= left:
                current[j] = up
                row[j] = 1
            else:
                current[j] = left
                row[j] = 2
        previous = current
    i, j, matches, columns = len(a), len(b), 0, 0
    while i > 0 or j > 0:
        direction = trace[i][j]
        if direction == 0:
            matches += a[i - 1] == b[j - 1]
            i, j = i - 1, j - 1
        elif direction == 1:
            i -= 1
        else:
            j -= 1
        columns += 1
    return ProteinAlignment(matches, columns, len(a), len(b))


def compare_genomes(a: Genome, b: Genome, k: int = 12) -> GenomeComparison:
    """Compare two annotated genomes by k-mer similarity and per-gene protein identity."""
    genes: list[GeneComparison] = []
    for name in Gene:
        protein_a, protein_b = a.protein(name), b.protein(name)
        if protein_a is not None and protein_b is not None:
            genes.append(GeneComparison(name, align_proteins(protein_a, protein_b)))
    return GenomeComparison(
        id_a=a.id,
        id_b=b.id,
        length_a=len(a),
        length_b=len(b),
        kmers=kmer_similarity(a.bases, b.bases, k),
        genes=tuple(genes),
    )


def compare_all(genomes: list[Genome], k: int = 12) -> list[GenomeComparison]:
    """Compare every unordered pair of annotated genomes, in input order."""
    return [compare_genomes(a, b, k) for a, b in combinations(genomes, 2)]
