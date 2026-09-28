"""Assign each corpus paper a single class from its title, abstract, and indexing."""

from __future__ import annotations

import re

from vsv.literature.analysis.records import Classification, PaperKind
from vsv.literature.records import Paper

_RECOMBINANT = re.compile(
    r"\brvsv\b|"
    r"recombinant (?:vesicular stomatitis|vsv)\b|"
    r"\bvsv[- ]based (?:vaccine|vector)\b|"
    r"\bvsv (?:vaccine )?vector\b|"
    r"\bvector(?:ized|ised)? vsv\b|"
    r"vsv(?:\W{0,3})?(?:Δ|∆|delta)\W{0,3}g\b",
    re.IGNORECASE,
)
_VSV_CARGO = re.compile(
    r"\bvsv\b.{0,60}\b(?:expressing|encoding)\b|"
    r"\b(?:expressing|encoding)\b.{0,40}\bvsv\b",
    re.IGNORECASE,
)
_ONCOLYTIC = re.compile(r"\boncoly", re.IGNORECASE)
_REPORTER = re.compile(r"\b(?:egfp|gfp|luciferase|reporter virus)\b", re.IGNORECASE)
_PSEUDOTYPE = re.compile(
    r"pseudotyp|"
    r"\bvsv-g\b.{0,40}\b(?:lentivir|pseudovir|envelope|murine leukemia)\b|"
    r"\blentiviral\b",
    re.IGNORECASE,
)
_PATHOGEN = re.compile(
    r"\b(?:ebola|zebov|marburg|marv|lassa|lasv|nipah|hendra|sars-cov-2|sars-cov2|"
    r"influenza|hiv|respiratory syncytial|measles|crimean-congo|hantavirus|zika|"
    r"chikungunya|rabies|norovirus)\b",
    re.IGNORECASE,
)
_EPIDEMIOLOGY = re.compile(
    r"\b(?:outbreaks?|livestock|cattle|equine|horses?|seroprevalence|serosurvey|"
    r"dairy|veterinary)\b",
    re.IGNORECASE,
)
_STRUCTURE = re.compile(
    r"cryo-electron|cryo-em|crystal structure|atomic (?:model|structure)|ångstrom|\bangstrom\b|"
    r"\bmorphology\b|surface structure",
    re.IGNORECASE,
)
_TITLE_SUBJECT = re.compile(r"vesicular stomatitis", re.IGNORECASE)
_VIROLOGY = re.compile(
    r"\b(?:polymerase|nucleocapsid|nucleoprotein|transcription|interferon|"
    r"replication cycle|glycoprotein structure)\b",
    re.IGNORECASE,
)


def classify(paper: Paper) -> Classification:
    """Classify ``paper`` from its title, abstract, publication types, and MeSH terms.

    A review keeps the subject class when the text supports one. A paper that matches no
    rule is ``other`` and ``unresolved``.
    """
    text = f"{paper.title}\n{paper.abstract}"
    review = any(item.casefold() == "review" for item in paper.pub_types)
    kind = _kind(text)
    if kind is PaperKind.OTHER and _TITLE_SUBJECT.search(paper.title):
        kind = PaperKind.BASIC_VIROLOGY
    if kind is PaperKind.OTHER and review:
        kind = PaperKind.REVIEW
    status = "unresolved" if kind is PaperKind.OTHER else "classified"
    return Classification(kind, review, status)


def _kind(text: str) -> PaperKind:
    if not text.strip():
        return PaperKind.OTHER
    recombinant = _recombinant(text)
    pathogen = _PATHOGEN.search(text) is not None
    if recombinant and pathogen:
        return PaperKind.VACCINE_VECTOR
    if recombinant and _ONCOLYTIC.search(text):
        return PaperKind.ONCOLYTIC
    if recombinant and _REPORTER.search(text):
        return PaperKind.REPORTER
    if recombinant:
        return PaperKind.VACCINE_VECTOR
    if _PSEUDOTYPE.search(text):
        return PaperKind.PSEUDOTYPE
    if _ONCOLYTIC.search(text):
        return PaperKind.ONCOLYTIC
    if _EPIDEMIOLOGY.search(text):
        return PaperKind.EPIDEMIOLOGY
    if _STRUCTURE.search(text):
        return PaperKind.STRUCTURE
    if _VIROLOGY.search(text):
        return PaperKind.BASIC_VIROLOGY
    return PaperKind.OTHER


def _recombinant(text: str) -> bool:
    """Whether the text describes a recombinant VSV genome, not another recombinant virus."""
    if _RECOMBINANT.search(text):
        return True
    return _VSV_CARGO.search(text) is not None
