"""Read a vaccine from an SBOL3 design file.

The file is the shareable record. Callers pass and receive a :class:`~vsv.design.vaccine.Vaccine`.
Leader and trailer are features in the file so the ordered genome is complete, and they are
not modules. A coding sequence is an ``sbol:Sequence`` only when the antigen has one.
"""

from __future__ import annotations

from pathlib import Path

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
from vsv.design.modules import Antigen, Insert, Replace, Viral
from vsv.design.sbol.terms import (
    _ACTIVITY,
    _ALIAS,
    _ANTIGEN_CLASS,
    _ASSAY,
    _CHALLENGE_AMOUNT,
    _CHALLENGE_DOSE,
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
    _STAGE,
    _STRAIN,
    _TIMING,
    _VACCINE_AMOUNT,
    _VACCINE_UNIT,
    _attachment_sources,
    _attachments,
    _citation_key,
    _dose_unit,
    _index_of,
    _optional_int,
    _optional_suffix,
    _suffixed,
    _text,
    _values,
)
from vsv.design.vaccine import GlycoproteinForm, Replication, Vaccine
from vsv.errors import InputError
from vsv.genes import Gene
from vsv.sbol import NAMESPACE, ROLE, ROLE_LEADER, ROLE_TRAILER

sbol3.set_namespace(NAMESPACE)


def read(path: str | Path) -> Vaccine:
    """Read one SBOL3 design file. The result is a :class:`Vaccine`."""
    document = _open(path)
    return _vaccine_of(document)


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


def _form_of(feature: sbol3.Feature) -> GlycoproteinForm:
    text = _text(feature, _FORM)
    if not text:
        return GlycoproteinForm.UNSPECIFIED
    try:
        return GlycoproteinForm(text)
    except ValueError as error:
        raise InputError(f"{feature.display_id}: unknown glycoprotein form {text!r}.") from error
