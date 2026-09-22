"""The five VSV genes and the sequences that can be read from each.

VSV transcribes its genes in a fixed order from the 3' end of the genome: Nucleoprotein,
Phosphoprotein, MatrixProtein, Glycoprotein, LargePolymerase. Those are the classical N, P, M, G,
and L genes.
"""

from enum import StrEnum


class Gene(StrEnum):
    """A VSV gene, declared in genome order."""

    NUCLEOPROTEIN = "Nucleoprotein"
    PHOSPHOPROTEIN = "Phosphoprotein"
    MATRIX_PROTEIN = "MatrixProtein"
    GLYCOPROTEIN = "Glycoprotein"
    LARGE_POLYMERASE = "LargePolymerase"

    @property
    def short(self) -> str:
        """The one-letter CLI name: ``n``, ``p``, ``m``, ``g``, or ``l``."""
        return _SHORT_NAMES[self]


_SHORT_NAMES: dict[Gene, str] = {
    Gene.NUCLEOPROTEIN: "n",
    Gene.PHOSPHOPROTEIN: "p",
    Gene.MATRIX_PROTEIN: "m",
    Gene.GLYCOPROTEIN: "g",
    Gene.LARGE_POLYMERASE: "l",
}


class Product(StrEnum):
    """A sequence that can be read from a gene."""

    GENE = "gene"
    CDS = "cds"
    PROTEIN = "protein"
