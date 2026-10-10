"""Load the curated VSV vaccine library.

Each vaccine is a directory with ``design.ttl`` and a provenance note. The file describes
the ordered genome and the papers. A nucleotide sequence is present only when an antigen
feature has an ``sbol:Sequence``. The shelf cards have none.
"""

from __future__ import annotations

from pathlib import Path

from vsv.design.classes import antigen_matches
from vsv.design.evidence import EvidenceLevel, Host, Phenotype, Timing, evidence_level, host_matches
from vsv.design.sbol import read
from vsv.design.vaccine import GlycoproteinForm, Vaccine
from vsv.errors import InputError

LIBRARY_ROOT = Path(__file__).resolve().parents[3] / "data" / "designs"


def load(path: Path) -> Vaccine:
    """Load one design directory, or a ``design.ttl`` file."""
    turtle = path / "design.ttl" if path.is_dir() else path
    if not turtle.is_file():
        raise InputError(f"{path}: missing design.ttl.")
    return read(turtle)


def library(root: Path | None = None) -> list[Vaccine]:
    """Load every design directory under ``root``, in directory-name order."""
    directory = LIBRARY_ROOT if root is None else root
    if not directory.is_dir():
        raise InputError(f"{directory}: design library not found.")
    records = [load(path) for path in sorted(directory.iterdir()) if path.is_dir()]
    identifiers = [record.id for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise InputError(f"{directory}: duplicate design ids.")
    return records


def select(
    records: list[Vaccine] | tuple[Vaccine, ...],
    *,
    antigen_class: str | None = None,
    phenotype: Phenotype | None = None,
    host: Host | None = None,
    min_evidence: EvidenceLevel | None = None,
    strain: str | None = None,
    form: GlycoproteinForm | None = None,
    timing: Timing | None = None,
) -> list[Vaccine]:
    """Vaccines that meet every filter that is set.

    Evidence is judged from every observation, before phenotype, host, and timing filters.
    ``antigen_class`` matches that class and its subclasses.
    """
    found: list[Vaccine] = []
    for vaccine in records:
        level = evidence_level([observation.host for observation in vaccine.observations])
        if min_evidence is not None and level.rank < min_evidence.rank:
            continue
        classes = [antigen.antigen_class for antigen in vaccine.antigens]
        if antigen_class is not None and not any(
            antigen_matches(slug, antigen_class) for slug in classes
        ):
            continue
        if strain is not None and strain not in {antigen.strain for antigen in vaccine.antigens}:
            continue
        if form is not None and form not in {antigen.form for antigen in vaccine.antigens}:
            continue
        narrowed = phenotype is not None or host is not None or timing is not None
        matched = any(
            (phenotype is None or observation.phenotype is phenotype)
            and (host is None or host_matches(observation.host, host))
            and (timing is None or observation.timing is timing)
            for observation in vaccine.observations
        )
        if narrowed and not matched:
            continue
        found.append(vaccine)
    return found
