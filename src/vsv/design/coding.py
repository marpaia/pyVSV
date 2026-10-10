"""Rules for one antigen coding sequence.

These rules read the antigen's class, form, and replication. They do not cut the
Indiana antigenome.
"""

from __future__ import annotations

import re

from vsv.design.classes import class_label, editing_length
from vsv.design.evidence import Violation
from vsv.design.vaccine import GlycoproteinForm, Replication
from vsv.genome.sequence import STOP_CODONS, hydrophobic_window, translate

_ADENOSINE_RUN = re.compile(r"A{7,8}")


def coding_violations(
    identifier: str,
    coding: str,
    *,
    replication: Replication,
    form: GlycoproteinForm,
    antigen_class: str,
    replaced: bool,
) -> tuple[Violation, ...]:
    """ORF, entry, anchor, and editing-site failures. An empty result passes."""
    violations: list[Violation] = []
    closed = _closed_orf(coding)
    if not closed:
        violations.append(
            Violation(
                "open-reading-frame",
                f"{identifier}: the coding sequence is not a closed open reading frame.",
            )
        )
    competent_replacement = replication is Replication.COMPETENT and replaced
    if competent_replacement and form is GlycoproteinForm.SOLUBLE:
        violations.append(
            Violation(
                "soluble-entry",
                f"{identifier}: a soluble glycoprotein does not supply entry.",
            )
        )
    if closed and competent_replacement and form.needs_anchor():
        protein = translate(coding, to_stop=True)
        if not hydrophobic_window(protein):
            violations.append(
                Violation(
                    "transmembrane-anchor",
                    f"{identifier}: a replication-competent glycoprotein replacement "
                    "has no transmembrane anchor.",
                )
            )
    expected = editing_length(antigen_class, form.value)
    if expected is not None and _editing_run(coding) != str(expected):
        violations.append(
            Violation(
                "editing-site",
                f"{identifier}: {class_label(antigen_class)} editing site "
                f"is a run of {expected} A's.",
            )
        )
    return tuple(violations)


def _editing_run(coding: str) -> str | None:
    """The single run of seven or eight A's, or ``None`` when the ORF does not have one.

    A run of nine or more is one greedy match of eight.
    """
    runs = _ADENOSINE_RUN.findall(coding)
    if len(runs) != 1:
        return None
    return str(len(runs[0]))


def _closed_orf(coding: str) -> bool:
    if len(coding) < 6 or len(coding) % 3 != 0 or not coding.startswith("ATG"):
        return False
    if coding[-3:] not in STOP_CODONS:
        return False
    body = range(0, len(coding) - 3, 3)
    return all(coding[index : index + 3] not in STOP_CODONS for index in body)
