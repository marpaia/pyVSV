"""Read and write a vaccine as an SBOL3 design file.

The file is the shareable record. Callers pass and receive a :class:`~vsv.design.vaccine.Vaccine`.
Leader and trailer are features in the file so the ordered genome is complete, and they are
not modules. A coding sequence is an ``sbol:Sequence`` only when the antigen has one.
"""

from __future__ import annotations

import math
from itertools import pairwise
from pathlib import Path
from urllib.parse import quote

import sbol3

from vsv.design.evidence import (
    Assay,
    Challenge,
    DoseUnit,
    Host,
    Implementation,
    Observation,
    Outcome,
    Phenotype,
    Route,
    Stage,
    Timing,
)
from vsv.design.vaccine import (
    Antigen,
    GlycoproteinForm,
    Insert,
    Replace,
    Replication,
    Vaccine,
    Viral,
)
from vsv.errors import InputError
from vsv.genes import Gene
from vsv.sbol import NAMESPACE, ROLE, ROLE_LEADER, ROLE_TRAILER, SO_GENE, SO_REGION, pin

sbol3.set_namespace(NAMESPACE)

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


def read(path: str | Path) -> Vaccine:
    """Read one SBOL3 design file. The result is a :class:`Vaccine`."""
    document = _open(path)
    return _vaccine_of(document)


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


def read_implementations(
    document: sbol3.Document, component: sbol3.Component
) -> list[Implementation]:
    """Implementations built from ``component``, in stored order."""
    found: list[tuple[int, Implementation]] = []
    for obj in document.objects:
        if not isinstance(obj, sbol3.Implementation) or str(obj.built) != str(component.identity):
            continue
        if not obj.name:
            raise InputError(f"{obj.display_id}: an implementation needs a name.")
        stage = _suffixed(_values(obj, sbol3.SBOL_ROLE), _STAGE, "stage")
        citations = _attachment_sources(document, obj)
        if len(citations) != 1:
            raise InputError(f"{obj.display_id}: an implementation needs one citation.")
        found.append(
            (
                _index_of(obj),
                Implementation(
                    str(obj.name),
                    Stage(stage),
                    citations[0],
                    str(obj.description or ""),
                    _text(obj, _PLASMID),
                    _text(obj, _CLONING_SITES),
                ),
            )
        )
    found.sort(key=lambda item: item[0])
    return [implementation for _index, implementation in found]


def read_observations(
    document: sbol3.Document,
    component: sbol3.Component,
    implementations: list[Implementation],
) -> list[Observation]:
    """Observations that use an implementation of ``component``, in stored order."""
    stored = [
        obj
        for obj in document.objects
        if isinstance(obj, sbol3.Implementation) and str(obj.built) == str(component.identity)
    ]
    by_identity = {str(item.identity): item for item in stored}
    names = {item.id for item in implementations}
    found: list[tuple[int, Observation]] = []
    for activity in document.objects:
        if not isinstance(activity, sbol3.Activity) or not activity.usage or not activity.name:
            continue
        if any(str(kind).startswith(f"{_ACTIVITY}/") for kind in activity.types):
            continue
        implementation = by_identity.get(str(activity.usage[0].entity))
        if implementation is None or str(implementation.name) not in names:
            continue
        terms = [str(kind) for kind in activity.types]
        citations = _attachment_sources(document, activity)
        if len(citations) != 1:
            raise InputError(f"{activity.display_id}: an observation needs one citation.")
        found.append(
            (
                _index_of(activity),
                Observation(
                    str(activity.name),
                    str(implementation.name),
                    Phenotype(_suffixed(terms, _PHENOTYPE, "phenotype")),
                    Host(_suffixed(terms, _HOST, "host")),
                    Assay(_suffixed(terms, _ASSAY, "assay")),
                    Outcome(_suffixed(terms, _OUTCOME, "outcome")),
                    citations[0],
                    str(activity.description or ""),
                    Route(_optional_suffix(terms, _ROUTE, Route.UNKNOWN.value)),
                    Timing(_optional_suffix(terms, _TIMING, Timing.UNKNOWN.value)),
                    _optional_int(activity, _DAY),
                    _text(activity, _VACCINE_AMOUNT),
                    _dose_unit(activity, _VACCINE_UNIT),
                    _challenge_from(activity),
                ),
            )
        )
    found.sort(key=lambda item: item[0])
    return [observation for _index, observation in found]


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


def _vaccine_of(document: sbol3.Document) -> Vaccine:
    component = _plan_component(document)
    modules = tuple(_module_of(document, feature) for feature in _ordered(component))
    implementations = read_implementations(document, component)
    return Vaccine(
        str(component.display_id),
        str(component.name or ""),
        _replication(component),
        modules,
        parent=_parent(component),
        edit=_edit(document, component),
        citations=_citations(document, component),
        aliases=tuple(_values(component, _ALIAS)),
        implementations=tuple(implementations),
        observations=tuple(read_observations(document, component, implementations)),
        realization=_realization(document, component),
    )


def _plan_component(document: sbol3.Document) -> sbol3.Component:
    found = [
        obj
        for obj in document.objects
        if isinstance(obj, sbol3.Component) and sbol3.SBO_DNA in {str(kind) for kind in obj.types}
    ]
    if len(found) != 1:
        raise InputError(f"Expected one vaccine, found {len(found)}.")
    component = found[0]
    if not component.name:
        raise InputError(f"{component.display_id}: a vaccine needs a name.")
    return component


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


def _module_of(document: sbol3.Document, feature: sbol3.Feature) -> Viral | Antigen:
    roles = [str(role) for role in feature.roles]
    classes = [role for role in roles if role.startswith(f"{_ANTIGEN_CLASS}/")]
    if classes:
        if len(classes) != 1:
            raise InputError(f"{feature.display_id}: expected one antigen class.")
        return Antigen(
            str(feature.display_id),
            str(feature.name or ""),
            classes[0].rsplit("/", 1)[-1],
            _text(feature, _STRAIN),
            _text(feature, _ISOLATE_ACCESSION),
            _form_of(feature),
            _coding(document, feature),
        )
    try:
        gene = Gene(str(feature.display_id))
    except ValueError as error:
        raise InputError(f"{feature.display_id}: expected a VSV gene or an antigen.") from error
    return Viral(gene)


def _ordered(component: sbol3.Component) -> list[sbol3.Feature]:
    features = [
        feature for feature in component.features if isinstance(feature, sbol3.LocalSubComponent)
    ]
    by_identity = {str(feature.identity): feature for feature in features}
    outgoing: dict[str, str] = {}
    incoming: set[str] = set()
    for constraint in component.constraints:
        if str(constraint.restriction) != sbol3.SBOL_STRICTLY_PRECEDES:
            continue
        subject, target = str(constraint.subject), str(constraint.object)
        if subject not in by_identity or target not in by_identity:
            continue
        if subject in outgoing or target in incoming:
            raise InputError(f"{component.display_id}: module order is not a single chain.")
        outgoing[subject] = target
        incoming.add(target)
    linked = set(outgoing) | incoming
    if not linked:
        ordered = list(features)
    else:
        starts = [
            feature
            for feature in features
            if str(feature.identity) in linked and str(feature.identity) not in incoming
        ]
        if len(starts) != 1:
            raise InputError(f"{component.display_id}: module order is not a single chain.")
        ordered = []
        seen: set[str] = set()
        current: sbol3.Feature | None = starts[0]
        while current is not None:
            identity = str(current.identity)
            if identity in seen:
                raise InputError(f"{component.display_id}: module order cycles.")
            seen.add(identity)
            ordered.append(current)
            nxt = outgoing.get(identity)
            current = by_identity.get(nxt) if nxt is not None else None
        if len(ordered) != len(linked):
            raise InputError(f"{component.display_id}: module order does not include every module.")
    return [feature for feature in ordered if not _is_terminus(feature)]


def _is_terminus(feature: sbol3.Feature) -> bool:
    keys = {str(role) for role in feature.roles}
    return ROLE_LEADER in keys or ROLE_TRAILER in keys


def _edit(document: sbol3.Document, component: sbol3.Component) -> Replace | Insert | None:
    activities = []
    for uri in component.generated_by:
        activity = document.find(str(uri))
        if isinstance(activity, sbol3.Activity) and any(
            str(kind).startswith(f"{_ACTIVITY}/") for kind in activity.types
        ):
            activities.append(activity)
    if not activities:
        return None
    if len(activities) != 1:
        raise InputError(f"{component.display_id}: expected one derivation.")
    activity = activities[0]
    kind = next(
        str(item).rsplit("/", 1)[-1]
        for item in activity.types
        if str(item).startswith(f"{_ACTIVITY}/")
    )
    module = _text(activity, _MODULE)
    slot = _text(activity, _SLOT)
    if not module:
        raise InputError(f"{component.display_id}: a derivation names the module.")
    if kind == "replace-module":
        return Replace(module)
    if kind == "insert-module" and slot == "g":
        return Insert("g", module)
    if kind == "insert-module" and slot == "leader":
        return Insert("leader", module)
    raise InputError(f"{component.display_id}: unknown derivation {kind} at {slot!r}.")


def _replication(component: sbol3.Component) -> Replication:
    prefix = f"{ROLE}/"
    matched = [
        str(role)[len(prefix) :]
        for role in component.roles
        if str(role).startswith(prefix)
        and str(role)[len(prefix) :] in {item.value for item in Replication}
    ]
    if len(matched) != 1:
        raise InputError(f"{component.display_id}: expected one replication class.")
    return Replication(matched[0])


def _parent(component: sbol3.Component) -> str | None:
    values = [str(value) for value in component.derived_from]
    if not values:
        return None
    if len(values) != 1 or not values[0].startswith(f"{NAMESPACE}/"):
        raise InputError(f"{component.display_id}: expected one parent in {NAMESPACE}.")
    return values[0].rsplit("/", 1)[-1]


def _citations(document: sbol3.Document, component: sbol3.Component) -> tuple[str, ...]:
    attachments = [
        attachment
        for attachment in _attachments(document, component)
        if str(attachment.display_id) != "realization"
    ]
    attachments.sort(key=_citation_key)
    return tuple(str(attachment.source) for attachment in attachments)


def _realization(document: sbol3.Document, component: sbol3.Component) -> str | None:
    found = [
        attachment
        for attachment in _attachments(document, component)
        if str(attachment.display_id) == "realization"
    ]
    if not found:
        return None
    if len(found) != 1 or not found[0].name:
        raise InputError(f"{component.display_id}: expected one realization.")
    return str(found[0].name)


def _coding(document: sbol3.Document, feature: sbol3.Feature) -> str:
    uris = _values(feature, sbol3.SBOL_SEQUENCES)
    if not uris:
        return ""
    if len(uris) != 1:
        raise InputError(f"{feature.display_id}: expected one coding sequence.")
    sequence = document.find(uris[0])
    if not isinstance(sequence, sbol3.Sequence) or not sequence.elements:
        raise InputError(f"{feature.display_id}: missing coding sequence.")
    return str(sequence.elements).upper()


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


def _challenge_from(activity: sbol3.Activity) -> Challenge | None:
    virus = _text(activity, _CHALLENGE_VIRUS)
    strain = _text(activity, _CHALLENGE_STRAIN)
    route_text = _text(activity, _CHALLENGE_ROUTE)
    amount = _text(activity, _CHALLENGE_AMOUNT)
    if not amount:
        amount = _text(activity, _CHALLENGE_DOSE)
    unit = _dose_unit(activity, _CHALLENGE_UNIT)
    route = Route(route_text) if route_text in Route else Route.UNKNOWN
    stated = bool(virus or strain or amount or route_text) or unit is not DoseUnit.UNSPECIFIED
    if not stated:
        return None
    return Challenge(virus, strain, route, amount, unit)


def _open(path: str | Path) -> sbol3.Document:
    sbol3.set_namespace(NAMESPACE)
    document = sbol3.Document()
    try:
        document.read(str(path))
    except (OSError, ValueError) as error:
        raise InputError(f"{path}: {error}") from error
    return document


def _cite(document: sbol3.Document, owner: sbol3.TopLevel, display: str, source: str) -> None:
    attachment = sbol3.Attachment(display, source)
    document.add(attachment)
    owner.attachments = [*list(owner.attachments), attachment]


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


def _form_of(feature: sbol3.Feature) -> GlycoproteinForm:
    text = _text(feature, _FORM)
    if not text:
        return GlycoproteinForm.UNSPECIFIED
    try:
        return GlycoproteinForm(text)
    except ValueError as error:
        raise InputError(f"{feature.display_id}: unknown glycoprotein form {text!r}.") from error


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


def _local_id(identifier: str) -> str:
    cleaned = "".join(
        character if character.isalnum() or character == "_" else "_" for character in identifier
    )
    if not cleaned or cleaned[0].isdigit():
        cleaned = f"_{cleaned}"
    return cleaned
