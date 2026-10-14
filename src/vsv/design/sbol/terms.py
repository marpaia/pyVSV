"""Properties stored in an SBOL3 vaccine-plan file."""

from __future__ import annotations

import sbol3

from vsv.design.evidence import DoseUnit
from vsv.errors import InputError
from vsv.sbol import NAMESPACE

_SO_VIRAL_SEQUENCE = "https://identifiers.org/SO:0001041"
_ACTIVITY = f"{NAMESPACE}/activity"
_ANTIGEN_CLASS = f"{NAMESPACE}/antigen-class"
_BACKBONE = f"{NAMESPACE}/backbone/VSV-Indiana"
_INDEX = f"{NAMESPACE}/index"
_STAGE = f"{NAMESPACE}/stage"
_PHENOTYPE = f"{NAMESPACE}/phenotype"
_HOST = f"{NAMESPACE}/host"
_ASSAY = f"{NAMESPACE}/assay"
_OUTCOME = f"{NAMESPACE}/outcome"
_ROUTE = f"{NAMESPACE}/route"
_TIMING = f"{NAMESPACE}/timing"
_CHALLENGE_VIRUS = f"{NAMESPACE}/challenge-virus"
_CHALLENGE_STRAIN = f"{NAMESPACE}/challenge-strain"
_CHALLENGE_ROUTE = f"{NAMESPACE}/challenge-route"
_CHALLENGE_AMOUNT = f"{NAMESPACE}/challenge-amount"
_CHALLENGE_UNIT = f"{NAMESPACE}/challenge-unit"
_CHALLENGE_DOSE = f"{NAMESPACE}/challenge-dose"
_DAY = f"{NAMESPACE}/day"
_VACCINE_AMOUNT = f"{NAMESPACE}/vaccine-amount"
_VACCINE_UNIT = f"{NAMESPACE}/vaccine-unit"
_ACCESSION = f"{NAMESPACE}/accession"
_ALIAS = f"{NAMESPACE}/alias"
_FORM = f"{NAMESPACE}/glycoprotein-form"
_STRAIN = f"{NAMESPACE}/strain"
_ISOLATE_ACCESSION = f"{NAMESPACE}/isolate-accession"
_PLASMID = f"{NAMESPACE}/plasmid"
_CLONING_SITES = f"{NAMESPACE}/cloning-sites"
_SLOT = f"{NAMESPACE}/slot"
_MODULE = f"{NAMESPACE}/module"


def _attachments(document: sbol3.Document, owner: sbol3.TopLevel) -> list[sbol3.Attachment]:
    found = []
    for uri in owner.attachments:
        attachment = document.find(str(uri))
        if not isinstance(attachment, sbol3.Attachment):
            raise InputError(f"{owner.display_id}: missing attachment {uri}.")
        found.append(attachment)
    return found


def _attachment_sources(document: sbol3.Document, owner: sbol3.TopLevel) -> tuple[str, ...]:
    return tuple(str(attachment.source) for attachment in _attachments(document, owner))


def _values(obj: sbol3.SBOLObject, predicate: str) -> list[str]:
    return [str(value) for value in obj._properties.get(predicate, [])]


def _text(obj: sbol3.SBOLObject, predicate: str) -> str:
    values = _values(obj, predicate)
    if len(values) > 1:
        raise InputError(f"{obj.display_id}: expected one {predicate}.")
    return values[0] if values else ""


def _write_text(obj: sbol3.SBOLObject, predicate: str, value: str) -> None:
    if value:
        sbol3.TextProperty(obj, predicate, 1, 1, initial_value=value)


def _optional_int(obj: sbol3.SBOLObject, predicate: str) -> int | None:
    values = _values(obj, predicate)
    if not values:
        return None
    if len(values) != 1:
        raise InputError(f"{obj.display_id}: expected one {predicate}.")
    return int(values[0])


def _set_index(obj: sbol3.TopLevel, index: int) -> None:
    sbol3.IntProperty(obj, _INDEX, 1, 1, initial_value=index)


def _index_of(obj: sbol3.SBOLObject) -> int:
    values = _values(obj, _INDEX)
    if len(values) != 1:
        raise InputError(f"{obj.display_id}: expected an index.")
    return int(values[0])


def _dose_unit(obj: sbol3.SBOLObject, predicate: str) -> DoseUnit:
    text = _text(obj, predicate)
    if not text:
        return DoseUnit.UNSPECIFIED
    try:
        return DoseUnit(text)
    except ValueError as error:
        raise InputError(f"{obj.display_id}: unknown dose unit {text!r}.") from error


def _suffixed(values: list[str], namespace: str, label: str) -> str:
    matched = [value[len(namespace) + 1 :] for value in values if value.startswith(f"{namespace}/")]
    if len(matched) != 1:
        raise InputError(f"Expected one {label}.")
    return matched[0]


def _optional_suffix(values: list[str], namespace: str, default: str) -> str:
    matched = [value[len(namespace) + 1 :] for value in values if value.startswith(f"{namespace}/")]
    if len(matched) > 1:
        raise InputError(f"Expected at most one {namespace}.")
    return matched[0] if matched else default


def _citation_key(attachment: sbol3.Attachment) -> tuple[int, str]:
    name = str(attachment.display_id)
    suffix = name.rsplit("_citation_", 1)[-1]
    return (int(suffix) if suffix.isdigit() else 10**9, name)
