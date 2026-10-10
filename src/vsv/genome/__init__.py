"""VSV genomes, their annotations, and the files they are read from."""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from os import PathLike
from pathlib import Path

from vsv.errors import InputError
from vsv.genes import Gene
from vsv.genome.formats import read_records
from vsv.genome.record import Genome

__all__ = [
    "GeneWarning",
    "Genome",
    "gene_problems",
    "load",
    "load_all",
    "read_record",
]


class GeneWarning(UserWarning):
    """A genome whose five VSV genes could not all be identified."""


def gene_problems(genome: Genome) -> list[str]:
    """What keeps ``genome`` from having all five genes with coding sequences, if anything."""
    found = genome.named
    absent = [gene for gene in Gene if gene not in found]
    uncoded = [gene for gene, annotation in found.items() if annotation.cds is None]
    problems = []
    if absent:
        problems.append(f"genes not identified: {', '.join(absent)}")
    if uncoded:
        problems.append(f"genes without a coding sequence: {', '.join(uncoded)}")
    return problems


def load_all(paths: Sequence[str | PathLike[str]]) -> list[Genome]:
    """Read every genome from FASTA, GenBank, or SBOL files, with genes identified.

    Each genome that lacks one of the five VSV genes raises a :class:`GeneWarning`.
    """
    genomes = read_records([Path(path) for path in paths])
    for genome in genomes:
        for problem in gene_problems(genome):
            warnings.warn(GeneWarning(f"{genome.id}: {problem}"), stacklevel=2)
    return genomes


def load(path: str | PathLike[str]) -> Genome:
    """Read the single genome in a FASTA, GenBank, or SBOL file, with genes identified.

    Raises :class:`~vsv.errors.InputError` if the file holds more than one record.
    """
    genomes = load_all([path])
    if len(genomes) != 1:
        raise InputError(f"{path}: expected one record, found {len(genomes)}; use load_all.")
    return genomes[0]


def read_record(path: str | PathLike[str]) -> Genome:
    """Read one FASTA, GenBank, or SBOL record and identify its VSV genes.

    Raises :class:`~vsv.errors.InputError` if the file holds more than one record.
    A genome that lacks one of the five genes does not warn; :func:`load` does.
    """
    genomes = read_records([Path(path)])
    if len(genomes) != 1:
        raise InputError(f"{path}: expected one record, found {len(genomes)}; use load_all.")
    return genomes[0]
