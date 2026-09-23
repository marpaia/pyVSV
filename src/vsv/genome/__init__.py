"""VSV genomes, their annotations, and the files they are read from."""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from os import PathLike
from pathlib import Path

import sbol3

from vsv.errors import InputError
from vsv.genes import Gene
from vsv.genome.annotation import annotate, is_annotated, name_curated_genes
from vsv.genome.formats import read_genomes
from vsv.genome.payload import AntigenClass, AntigenFeature
from vsv.genome.sbol import GeneAnnotation, Genome

__all__ = [
    "AntigenClass",
    "AntigenFeature",
    "GeneAnnotation",
    "GeneWarning",
    "Genome",
    "gene_problems",
    "load",
    "load_all",
    "prepare",
    "read_record",
]


class GeneWarning(UserWarning):
    """A genome whose five VSV genes could not all be identified."""


def gene_problems(genome: Genome) -> list[str]:
    """What keeps ``genome`` from having all five genes with coding sequences, if anything."""
    found = genome.genes
    absent = [gene for gene in Gene if gene not in found]
    uncoded = [gene for gene, annotation in found.items() if annotation.cds_feature is None]
    problems = []
    if absent:
        problems.append(f"genes not identified: {', '.join(absent)}")
    if uncoded:
        problems.append(f"genes without a coding sequence: {', '.join(uncoded)}")
    return problems


def prepare(genome: Genome) -> list[str]:
    """Identify the genes of ``genome`` in place and return the remaining :func:`gene_problems`.

    A genome with gene or coding-sequence features keeps them and has its genes named by position.
    Any other genome is annotated from its sequence.
    """
    if is_annotated(genome):
        name_curated_genes(genome)
    else:
        annotate(genome)
    return gene_problems(genome)


def load_all(paths: Sequence[str | PathLike[str]]) -> tuple[sbol3.Document, list[Genome]]:
    """Read every genome from FASTA, GenBank, or SBOL3 files, with genes identified.

    Each genome that lacks one of the five VSV genes raises a :class:`GeneWarning`.
    """
    document, genomes = read_genomes([Path(path) for path in paths])
    for genome in genomes:
        for problem in prepare(genome):
            warnings.warn(GeneWarning(f"{genome.id}: {problem}"), stacklevel=2)
    return document, genomes


def load(path: str | PathLike[str]) -> Genome:
    """Read the single genome in a FASTA, GenBank, or SBOL3 file, with genes identified.

    Raises :class:`~vsv.errors.InputError` if the file holds more than one record.
    """
    _, genomes = load_all([path])
    if len(genomes) != 1:
        raise InputError(f"{path}: expected one record, found {len(genomes)}; use load_all.")
    return genomes[0]


def read_record(path: str | PathLike[str]) -> Genome:
    """Read one FASTA, GenBank, or SBOL3 record and leave its features unnamed as VSV genes.

    A GenBank or SBOL record keeps the features it arrived with. A FASTA record has none until
    they are added. Raises :class:`~vsv.errors.InputError` if the file holds more than one record.
    """
    _, genomes = read_genomes([Path(path)], name_genes=False)
    if len(genomes) != 1:
        raise InputError(f"{path}: expected one record, found {len(genomes)}; use load_all.")
    return genomes[0]
