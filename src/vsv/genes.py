"""The five VSV genes and the sequences that can be read from each.

VSV transcribes its genes in a fixed order from the 3' end of the genome: N, P, M, G, L.
"""

from enum import StrEnum


class Gene(StrEnum):
    """A VSV gene, declared in genome order."""

    N = "N"
    P = "P"
    M = "M"
    G = "G"
    L = "L"


class Product(StrEnum):
    """A sequence that can be read from a gene."""

    GENE = "gene"
    CDS = "cds"
    PROTEIN = "protein"
