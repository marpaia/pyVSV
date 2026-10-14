"""Write a vaccine plan as an SBOL3 design file."""

from __future__ import annotations

import math
from itertools import pairwise
from pathlib import Path
from typing import TYPE_CHECKING
from urllib.parse import quote

import sbol3

from vsv.design.evidence import (
    Challenge,
    DoseUnit,
    Implementation,
    Observation,
    Route,
    Timing,
)
from vsv.design.modules import Antigen, Replace, Viral
from vsv.design.sbol.terms import (
    _ACCESSION,
    _ACTIVITY,
    _ALIAS,
    _ANTIGEN_CLASS,
    _ASSAY,
    _BACKBONE,
    _CHALLENGE_AMOUNT,
    _CHALLENGE_ROUTE,
    _CHALLENGE_STRAIN,
    _CHALLENGE_UNIT,
    _CHALLENGE_VIRUS,
    _CLONING_SITES,
    _DAY,
    _FORM,
    _HOST,
    _ISOLATE_ACCESSION,
    _MODULE,
    _OUTCOME,
    _PHENOTYPE,
    _PLASMID,
    _ROUTE,
    _SLOT,
    _SO_VIRAL_SEQUENCE,
    _STAGE,
    _STRAIN,
    _TIMING,
    _VACCINE_AMOUNT,
    _VACCINE_UNIT,
    _set_index,
    _write_text,
)
from vsv.design.vector import GlycoproteinForm
from vsv.errors import InputError
from vsv.sbol import NAMESPACE, ROLE, ROLE_LEADER, ROLE_TRAILER, SO_GENE, SO_REGION, pin

if TYPE_CHECKING:
    from vsv.design.vaccine import Vaccine


def write(vaccine: Vaccine, path: str | Path) -> None:
    """Write ``vaccine`` as an SBOL3 design file."""
    _document(vaccine).write(str(path))


def write_implementation(
    document: sbol3.Document,
    component: sbol3.Component,
    owner_id: str,
    implementation: Implementation,
    index: int,
) -> sbol3.Implementation:
    """Add one implementation whose ``sbol:built`` property is ``component``."""
    local = _local_id(implementation.id)
    stored = sbol3.Implementation(
        f"{owner_id}__impl_{local}",
        built=component,
        name=implementation.id,
        description=implementation.note or None,
    )
    sbol3.URIProperty(
        stored, sbol3.SBOL_ROLE, 1, 1, initial_value=f"{_STAGE}/{implementation.stage.value}"
    )
    _set_index(stored, index)
    if implementation.plasmid:
        sbol3.TextProperty(stored, _PLASMID, 1, 1, initial_value=implementation.plasmid)
    if implementation.cloning_sites:
        sbol3.TextProperty(stored, _CLONING_SITES, 1, 1, initial_value=implementation.cloning_sites)
    document.add(stored)
    _cite(document, stored, f"{owner_id}__impl_{local}_citation", implementation.citation)
    return stored


def write_observation(
    document: sbol3.Document,
    owner_id: str,
    observation: Observation,
    implementation: sbol3.Implementation,
    index: int,
) -> None:
    """Add one observation as an Activity that uses ``implementation``."""
    local = _local_id(observation.id)
    types = [
        f"{_PHENOTYPE}/{observation.phenotype.value}",
        f"{_HOST}/{observation.host.value}",
        f"{_ASSAY}/{observation.assay.value}",
        f"{_OUTCOME}/{observation.outcome.value}",
    ]
    if observation.route is not Route.UNKNOWN:
        types.append(f"{_ROUTE}/{observation.route.value}")
    if observation.timing is not Timing.UNKNOWN:
        types.append(f"{_TIMING}/{observation.timing.value}")
    activity = sbol3.Activity(
        f"{owner_id}__obs_{local}",
        types=types,
        name=observation.id,
        description=observation.note or None,
        usage=[sbol3.Usage(implementation.identity)],
    )
    if observation.day is not None:
        sbol3.IntProperty(activity, _DAY, 1, 1, initial_value=observation.day)
    if observation.vaccine_amount:
        sbol3.TextProperty(
            activity, _VACCINE_AMOUNT, 1, 1, initial_value=observation.vaccine_amount
        )
    if observation.vaccine_unit is not DoseUnit.UNSPECIFIED:
        sbol3.TextProperty(
            activity, _VACCINE_UNIT, 1, 1, initial_value=observation.vaccine_unit.value
        )
    _write_challenge(activity, observation.challenge)
    _set_index(activity, index)
    document.add(activity)
    _cite(document, activity, f"{owner_id}__obs_{local}_citation", observation.citation)


def _document(vaccine: Vaccine) -> sbol3.Document:
    sbol3.set_namespace(NAMESPACE)
    document = sbol3.Document()
    component = sbol3.Component(
        vaccine.id,
        [sbol3.SBO_DNA],
        roles=[_SO_VIRAL_SEQUENCE, _BACKBONE, f"{ROLE}/{vaccine.replication.value}"],
        name=vaccine.name,
    )
    if vaccine.parent is not None:
        component.derived_from = [f"{NAMESPACE}/{vaccine.parent}"]
    if vaccine.aliases:
        sbol3.TextProperty(component, _ALIAS, 0, math.inf, initial_value=list(vaccine.aliases))
    document.add(component)
    features = [_terminus(component, "leader", ROLE_LEADER)]
    for module in vaccine.modules:
        feature = _module_feature(document, component, module)
        features.append(feature)
    features.append(_terminus(component, "trailer", ROLE_TRAILER))
    for earlier, later in pairwise(features):
        component.constraints.append(sbol3.Constraint(sbol3.SBOL_STRICTLY_PRECEDES, earlier, later))
    for index, citation in enumerate(vaccine.citations, start=1):
        attachment = sbol3.Attachment(f"{vaccine.id}_citation_{index}", citation)
        document.add(attachment)
        component.attachments.append(attachment.identity)
    if vaccine.realization:
        attachment = sbol3.Attachment(
            "realization",
            f"{_ACCESSION}/{quote(vaccine.realization, safe='')}",
            name=vaccine.realization,
        )
        document.add(attachment)
        component.attachments.append(attachment.identity)
    _write_edit(document, component, vaccine)
    stored: list[sbol3.Implementation] = []
    for index, implementation in enumerate(vaccine.implementations):
        stored.append(write_implementation(document, component, vaccine.id, implementation, index))
    by_name = {str(item.name): item for item in stored}
    for index, observation in enumerate(vaccine.observations):
        built = by_name.get(observation.implementation_id)
        if built is None:
            raise InputError(
                f"{observation.id} cites unknown implementation {observation.implementation_id}."
            )
        write_observation(document, vaccine.id, observation, built, index)
    return document


def _terminus(component: sbol3.Component, identifier: str, role: str) -> sbol3.LocalSubComponent:
    feature = sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[SO_REGION, role], name=identifier)
    pin(feature, identifier)
    component.features.append(feature)
    return feature


def _module_feature(
    document: sbol3.Document, component: sbol3.Component, module: Viral | Antigen
) -> sbol3.LocalSubComponent:
    if isinstance(module, Viral):
        roles = [SO_GENE]
        name = module.gene.value
    else:
        roles = [SO_GENE, f"{_ANTIGEN_CLASS}/{module.antigen_class}"]
        name = module.name
    feature = sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=roles, name=name)
    pin(feature, module.id)
    if isinstance(module, Antigen):
        _write_text(feature, _STRAIN, module.strain)
        _write_text(feature, _ISOLATE_ACCESSION, module.accession)
        if module.form is not GlycoproteinForm.UNSPECIFIED:
            _write_text(feature, _FORM, module.form.value)
        if module.coding:
            sequence = sbol3.Sequence(
                f"{component.display_id}__{module.id}_sequence",
                elements=module.coding,
                encoding=sbol3.IUPAC_DNA_ENCODING,
            )
            document.add(sequence)
            sbol3.URIProperty(feature, sbol3.SBOL_SEQUENCES, 0, 1, initial_value=sequence.identity)
    component.features.append(feature)
    return feature


def _write_edit(document: sbol3.Document, component: sbol3.Component, vaccine: Vaccine) -> None:
    edit = vaccine.edit
    if edit is None:
        return
    if isinstance(edit, Replace):
        activity_name = "replace-module"
        slot = "g"
        module = edit.module
    else:
        activity_name = "insert-module"
        slot = edit.after
        module = edit.module
    activity = sbol3.Activity(f"{vaccine.id}_derivation", types=[f"{_ACTIVITY}/{activity_name}"])
    sbol3.TextProperty(activity, _SLOT, 1, 1, initial_value=slot)
    sbol3.TextProperty(activity, _MODULE, 1, 1, initial_value=module)
    document.add(activity)
    component.generated_by = [activity.identity]
    if vaccine.citations:
        _cite(document, activity, f"{vaccine.id}_derivation_citation", vaccine.citations[0])


def _write_challenge(activity: sbol3.Activity, challenge: Challenge | None) -> None:
    if challenge is None:
        return
    _write_text(activity, _CHALLENGE_VIRUS, challenge.virus)
    _write_text(activity, _CHALLENGE_STRAIN, challenge.strain)
    if challenge.route is not Route.UNKNOWN:
        _write_text(activity, _CHALLENGE_ROUTE, challenge.route.value)
    _write_text(activity, _CHALLENGE_AMOUNT, challenge.amount)
    if challenge.unit is not DoseUnit.UNSPECIFIED:
        _write_text(activity, _CHALLENGE_UNIT, challenge.unit.value)


def _cite(document: sbol3.Document, owner: sbol3.TopLevel, display: str, source: str) -> None:
    attachment = sbol3.Attachment(display, source)
    document.add(attachment)
    owner.attachments = [*list(owner.attachments), attachment]


def _local_id(identifier: str) -> str:
    cleaned = "".join(
        character if character.isalnum() or character == "_" else "_" for character in identifier
    )
    if not cleaned or cleaned[0].isdigit():
        cleaned = f"_{cleaned}"
    return cleaned
