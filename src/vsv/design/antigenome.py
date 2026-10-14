"""Splice one coding sequence into the Indiana antigenome NC_001560.1.

A glycoprotein replacement leaves both flanking ``CT`` dinucleotides on the backbone. An insert
between the glycoprotein and the large polymerase places ``CT`` before the new gene-start. An
insert after the leader places ``CT`` after the new gene-end, so the gene-end is not joined to
the nucleoprotein gene-start. The papers' infectious clones stay on ``Implementation.plasmid``.
NC_001560.1 is the reference antigenome, not a rescue plasmid.
"""

from __future__ import annotations

from functools import cache
from pathlib import Path

from vsv.design.evidence import Violation
from vsv.errors import InputError
from vsv.genome.annotate import GENE_END_SIGNAL, annotate
from vsv.genome.record import Genome
from vsv.genome.sbol.read import read_genomes
from vsv.genome.sequence import Nucleotides

_ACCESSION = "NC_001560.1"
_GENE_START = "AACAGAGATC"
_INTERGENIC = "CT"
_L_GENE_START = "AACAGCAATC"
_N_GENE_START = "AACAGTAATC"
_G_START = 3048
_G_END = 4720
_GL_CUT = 4720
_LEADER_CUT = 50
SLOT_CITATION = {
    "g": "https://doi.org/10.1128/jvi.78.10.5458-5465.2004",
    "g-l": "https://doi.org/10.1128/jvi.71.8.5982-5989.1997",
    "leader": "https://doi.org/10.1128/jvi.00019-14",
}


def cassette(slot: str, coding: str) -> str:
    """Bases spliced into ``slot``.

    ``g`` is the glycoprotein transcription unit. ``g-l`` inserts between the glycoprotein and
    the large polymerase and prefixes ``CT``. ``leader`` inserts before the nucleoprotein and
    suffixes ``CT``.
    """
    Nucleotides(coding)
    body = _GENE_START + coding + GENE_END_SIGNAL
    if slot == "g":
        return body
    if slot == "g-l":
        return _INTERGENIC + body
    if slot == "leader":
        return body + _INTERGENIC
    raise InputError(f"Unknown slot {slot}.")


def splice(
    coding: str,
    slot: str,
    *,
    identifier: str = _ACCESSION,
    description: str = "",
) -> Genome:
    """Place ``coding`` in ``slot`` of NC_001560.1 and identify the transcription units.

    A replacement leaves five genes, so ``cds`` names nucleoprotein, phosphoprotein, matrix
    protein, glycoprotein, and large polymerase. An insert leaves
    six units, and those units keep the numbered names from the annotator.
    """
    bases = _indiana_sequence()
    _require_accession(bases)
    inserted = cassette(slot, coding)
    if slot == "g":
        start, end = _G_START, _G_END
    elif slot == "g-l":
        start, end = _GL_CUT, _GL_CUT
    elif slot == "leader":
        start, end = _LEADER_CUT, _LEADER_CUT
    else:
        raise InputError(f"Unknown slot {slot}.")
    return annotate(
        bases[:start] + inserted + bases[end:],
        identifier=identifier,
        description=description,
    )


def cassette_violations(identifier: str, coding: str, slot: str) -> tuple[Violation, ...]:
    """Gene-end and junction failures for one inserted coding sequence."""
    violations: list[Violation] = []
    unit = _GENE_START + coding + GENE_END_SIGNAL
    if not unit.endswith(GENE_END_SIGNAL) or unit.count(GENE_END_SIGNAL) != 1:
        violations.append(
            Violation(
                "gene-end",
                f"{identifier}: a gene-end signal occurs inside the transcription unit.",
            )
        )
    problem = _junction_problem(splice(coding, slot).bases, slot)
    if problem is not None:
        violations.append(Violation("junction", f"{identifier}: {problem}"))
    return tuple(violations)


def _junction_problem(bases: str, slot: str) -> str | None:
    """A sentence when a new junction lacks its intergenic dinucleotide."""
    if GENE_END_SIGNAL + "AACAG" in bases:
        return "a gene-end is joined directly to a gene-start."
    required = {
        "g": (
            GENE_END_SIGNAL + _INTERGENIC + _GENE_START,
            GENE_END_SIGNAL + _INTERGENIC + _L_GENE_START,
        ),
        "g-l": (
            GENE_END_SIGNAL + _INTERGENIC + _GENE_START,
            GENE_END_SIGNAL + _INTERGENIC + _L_GENE_START,
        ),
        "leader": (GENE_END_SIGNAL + _INTERGENIC + _N_GENE_START,),
    }.get(slot, ())
    if any(junction not in bases for junction in required):
        return "a gene junction is missing its intergenic dinucleotide."
    return None


def _require_accession(bases: str) -> None:
    expected = {
        (3035, 3048): "TATGAAAAAAACT",
        (3048, 3058): _GENE_START,
        (4709, 4720): GENE_END_SIGNAL,
        (4720, 4732): "CTAACAGCAATC",
        (50, 60): _N_GENE_START,
    }
    for (start, end), sequence in expected.items():
        if bases[start:end] != sequence:
            raise InputError(f"{_ACCESSION}: the cited signal at {start} is not on this sequence.")
    if len(bases) != 11161:
        raise InputError(f"{_ACCESSION}: expected 11161 nt, found {len(bases)}.")


@cache
def _indiana_sequence() -> str:
    path = _antigenome_path()
    genomes = read_genomes(path)
    if len(genomes) != 1 or genomes[0].id != _ACCESSION:
        raise InputError(f"{path}: expected one {_ACCESSION} record.")
    return genomes[0].bases


def _antigenome_path() -> Path:
    root = Path(__file__).resolve().parents[3]
    checkout = root / "data" / "raw" / "NC_001560.ttl"
    if (root / "pyproject.toml").is_file() and checkout.is_file():
        return checkout
    packaged = Path(__file__).resolve().parents[1] / "NC_001560.ttl"
    if packaged.is_file():
        return packaged
    raise InputError("The Indiana antigenome is not installed.")
