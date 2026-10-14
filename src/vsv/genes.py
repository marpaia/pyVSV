"""The five VSV genes and the sequences that can be read from each.

VSV transcribes its genes in a fixed order from the 3' end of the genome: nucleoprotein,
phosphoprotein, matrix protein, glycoprotein, and large polymerase.
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
        """Command-line name: ``n`` nucleoprotein, ``p`` phosphoprotein, ``m`` matrix protein,
        ``g`` glycoprotein, or ``l`` large polymerase.
        """
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
