"""An architecture as an SBOL3 document with no nucleotide sequence.

Module order is a chain of ``strictlyPrecedes`` constraints. A transcribed module with no locus
stays off that chain. The edit that produced a design is an Activity whose type is
``replace-module``, ``insert-module``, or ``delete-module`` and whose description is
``target=…;module=…;replication=…``. Implementations, observations, a realization accession, and
a linked model live in the same document.
"""

from __future__ import annotations

import math
from itertools import pairwise
from urllib.parse import quote, unquote

import sbol3

from vsv import NAMESPACE
from vsv.design.architecture import Architecture, ParentStatus, Replication
from vsv.design.edits import Delete, Edit, Insert, Replace
from vsv.design.evidence import (
    Assay,
    Design,
    Host,
    Implementation,
    Observation,
    Outcome,
    Phenotype,
    Stage,
    Unmodeled,
)
from vsv.design.model import attach_model
from vsv.design.modules import (
    Antigen,
    Leader,
    Module,
    ModuleKind,
    Reporter,
    Trailer,
    ViralGene,
)
from vsv.errors import InputError
from vsv.genome.sbol import (
    ROLE_LEADER,
    ROLE_TRAILER,
    SO_GENE,
    SO_REGION,
    SO_VIRAL_SEQUENCE,
    new_document,
    pin_display_id,
    role_key,
)

sbol3.set_namespace(NAMESPACE)

_ACTIVITY_NAMESPACE = f"{NAMESPACE}/activity"
_ANTIGEN_CLASS_NAMESPACE = f"{NAMESPACE}/antigen-class"
_BACKBONE_NAMESPACE = f"{NAMESPACE}/backbone"
_ROLE_NAMESPACE = f"{NAMESPACE}/role"
_ROLE_REPORTER = f"{_ROLE_NAMESPACE}/reporter"
_ROLE_REPLICATION_COMPETENT = f"{_ROLE_NAMESPACE}/replication-competent"
_ROLE_REPLICATION_SINGLE_CYCLE = f"{_ROLE_NAMESPACE}/single-cycle"
_REPLICATION_ROLE = {
    Replication.COMPETENT: _ROLE_REPLICATION_COMPETENT,
    Replication.SINGLE_CYCLE: _ROLE_REPLICATION_SINGLE_CYCLE,
}
_ROLE_REPLICATION = {role: replication for replication, role in _REPLICATION_ROLE.items()}
_ROLE_PARENT_UNKNOWN = f"{_ROLE_NAMESPACE}/parent-unknown"
_INDEX = f"{NAMESPACE}/index"
_STAGE = f"{NAMESPACE}/stage"
_PHENOTYPE = f"{NAMESPACE}/phenotype"
_HOST = f"{NAMESPACE}/host"
_ASSAY = f"{NAMESPACE}/assay"
_OUTCOME = f"{NAMESPACE}/outcome"
_ACCESSION = f"{NAMESPACE}/accession"
_ALIAS = f"{NAMESPACE}/alias"
_UNMODELED_VALUE = f"{NAMESPACE}/unmodeled/value/"
_UNMODELED_SOURCE = f"{NAMESPACE}/unmodeled-source/"
_DERIVATION_KEYS = frozenset({"target", "module", "replication"})


def _module_roles(module: Module) -> list[str]:
    match module:
        case Leader():
            return [SO_REGION, ROLE_LEADER]
        case Trailer():
            return [SO_REGION, ROLE_TRAILER]
        case Antigen(antigen_class=antigen_class):
            return [SO_GENE, f"{_ANTIGEN_CLASS_NAMESPACE}/{antigen_class}"]
        case Reporter():
            return [SO_GENE, _ROLE_REPORTER]
        case ViralGene():
            return [SO_GENE]


def _leaders(architecture: Architecture) -> list[Module]:
    return [module for module in architecture.modules if isinstance(module, Leader)]


def _trailers(architecture: Architecture) -> list[Module]:
    return [module for module in architecture.modules if isinstance(module, Trailer)]


def _require_genome_order(architecture: Architecture) -> tuple[list[Module], list[Module]]:
    """Refuse to serialize a module list whose known loci disagree with its order.

    The constraint chain is the leader, the transcribed modules whose loci are known, and the
    trailer. Modules with an unknown locus sit between that chain's genes and the trailer, in
    display-id order. Constraints record order only, so a list that is not already arranged
    that way would come back in a different order.
    """
    known: list[Module] = [
        module for module in architecture.transcribed if module.locus is not None
    ]
    unknown: list[Module] = [module for module in architecture.transcribed if module.locus is None]
    if [module.locus for module in known] != list(range(len(known))):
        raise InputError(f"{architecture.id}: modules are not in locus order.")
    if [module.id for module in unknown] != sorted(module.id for module in unknown):
        raise InputError(f"{architecture.id}: modules are not in genome order.")
    leaders = _leaders(architecture)
    trailers = _trailers(architecture)
    ordered = leaders + known + unknown + trailers
    if [module.id for module in architecture.modules] != [module.id for module in ordered]:
        raise InputError(f"{architecture.id}: modules are not in genome order.")
    return leaders + known + trailers, unknown


def _activity(edit: Edit) -> str:
    match edit:
        case Replace():
            return "replace-module"
        case Insert():
            return "insert-module"
        case Delete():
            return "delete-module"


def _derivation_text(edit: Edit) -> str:
    match edit:
        case Delete(target=target, replication=replication):
            parts = [f"target={target}"]
        case Replace(target=target, module=module, replication=replication):
            parts = [f"target={target}", f"module={module}"]
        case Insert(after=after, module=module, replication=replication):
            parts = [f"target={after}", f"module={module}"]
    if replication is not None:
        parts.append(f"replication={replication.value}")
    return ";".join(parts)


def _parse_derivation(activity: str, description: str) -> Edit:
    fields: dict[str, str] = {}
    for part in description.split(";"):
        if not part:
            continue
        key, separator, value = part.partition("=")
        if separator != "=" or key not in _DERIVATION_KEYS or not value or key in fields:
            raise InputError(f"Invalid derivation {description!r}.")
        fields[key] = value
    if "target" not in fields:
        raise InputError(f"Invalid derivation {description!r}.")
    replication = None
    if "replication" in fields:
        try:
            replication = Replication(fields["replication"])
        except ValueError as error:
            raise InputError(f"Invalid derivation {description!r}.") from error
    module = fields.get("module")
    match activity:
        case "delete-module":
            if module is not None:
                raise InputError("A deletion derivation has no module.")
            return Delete(fields["target"], replication)
        case "replace-module" | "insert-module":
            if module is None:
                raise InputError(f"A {activity} derivation names the resulting module.")
            if activity == "replace-module":
                return Replace(fields["target"], module, replication)
            return Insert(fields["target"], module, replication)
        case _:
            raise InputError(f"Unknown derivation {activity!r}.")


def to_document(architecture: Architecture) -> sbol3.Document:
    """Serialize an architecture as an SBOL3 document with no nucleotide sequence."""
    chain, _unknown = _require_genome_order(architecture)
    document = new_document()
    roles = [
        SO_VIRAL_SEQUENCE,
        f"{_BACKBONE_NAMESPACE}/{architecture.backbone}",
    ]
    if architecture.replication is not Replication.UNKNOWN:
        roles.append(_REPLICATION_ROLE[architecture.replication])
    if architecture.parent is ParentStatus.UNKNOWN:
        roles.append(_ROLE_PARENT_UNKNOWN)
    component = sbol3.Component(
        architecture.id,
        [sbol3.SBO_DNA],
        roles=roles,
        name=architecture.name,
    )
    if architecture.parent_id is not None:
        component.derived_from = [f"{NAMESPACE}/{architecture.parent_id}"]
    document.add(component)
    features: dict[str, sbol3.LocalSubComponent] = {}
    for module in architecture.modules:
        feature = sbol3.LocalSubComponent(
            [sbol3.SBO_DNA], roles=_module_roles(module), name=module.name
        )
        pin_display_id(feature, module.id)
        component.features.append(feature)
        features[module.id] = feature
    for earlier, later in pairwise(chain):
        component.constraints.append(
            sbol3.Constraint(sbol3.SBOL_STRICTLY_PRECEDES, features[earlier.id], features[later.id])
        )
    for index, citation in enumerate(architecture.citations, start=1):
        attachment = sbol3.Attachment(f"{architecture.id}_citation_{index}", citation)
        document.add(attachment)
        component.attachments.append(attachment.identity)
    if architecture.edit is not None:
        activity = sbol3.Activity(
            f"{architecture.id}_derivation",
            types=[f"{_ACTIVITY_NAMESPACE}/{_activity(architecture.edit)}"],
            description=_derivation_text(architecture.edit),
        )
        parent = architecture.parent_id
        if parent is None:
            raise InputError(f"{architecture.id}: a derivation requires a parent.")
        usage = sbol3.Usage(f"{NAMESPACE}/{parent}")
        pin_display_id(usage, "parent")
        activity.usage = [usage]
        document.add(activity)
        component.generated_by = [activity.identity]
    return document


def _kind_of(feature: sbol3.Feature) -> tuple[ModuleKind, str | None]:
    roles = [str(role) for role in feature.roles]
    keys = {role_key(role) for role in roles}
    if role_key(ROLE_LEADER) in keys:
        return ModuleKind.LEADER, None
    if role_key(ROLE_TRAILER) in keys:
        return ModuleKind.TRAILER, None
    antigen = [role for role in roles if role.startswith(f"{_ANTIGEN_CLASS_NAMESPACE}/")]
    if len(antigen) > 1:
        raise InputError(f"Feature {feature.identity} has more than one antigen class.")
    if antigen:
        return ModuleKind.ANTIGEN, antigen[0].rsplit("/", 1)[-1]
    if role_key(_ROLE_REPORTER) in keys:
        return ModuleKind.REPORTER, None
    return ModuleKind.VIRAL_GENE, None


def _ordered_features(
    component: sbol3.Component,
) -> tuple[list[sbol3.Feature], list[sbol3.Feature]]:
    """The precedence chain and the features that have no locus constraint."""
    features = [
        feature for feature in component.features if isinstance(feature, sbol3.LocalSubComponent)
    ]
    by_identity = {str(feature.identity): feature for feature in features}
    outgoing: dict[str, str] = {}
    incoming: set[str] = set()
    for constraint in component.constraints:
        if constraint.restriction != sbol3.SBOL_STRICTLY_PRECEDES:
            continue
        subject, target = str(constraint.subject), str(constraint.object)
        if subject not in by_identity or target not in by_identity:
            continue
        if subject in outgoing or target in incoming:
            raise InputError(f"{component.display_id}: module order is not a single chain.")
        outgoing[subject] = target
        incoming.add(target)
    linked = set(outgoing) | incoming
    free = [feature for feature in features if str(feature.identity) not in linked]
    if not linked:
        return [], free
    starts = [
        feature
        for feature in features
        if str(feature.identity) in linked and str(feature.identity) not in incoming
    ]
    if len(starts) != 1:
        raise InputError(f"{component.display_id}: module order is not a single chain.")
    ordered: list[sbol3.Feature] = []
    seen: set[str] = set()
    current: sbol3.Feature | None = starts[0]
    while current is not None:
        identity = str(current.identity)
        if identity in seen:
            raise InputError(f"{component.display_id}: module order cycles.")
        seen.add(identity)
        ordered.append(current)
        nxt = outgoing.get(identity)
        current = by_identity[nxt] if nxt is not None else None
    if len(ordered) != len(linked):
        raise InputError(f"{component.display_id}: module order does not include every module.")
    return ordered, free


def _component_attachments(
    document: sbol3.Document, component: sbol3.Component
) -> list[sbol3.Attachment]:
    attachments = []
    for uri in component.attachments:
        attachment = document.find(str(uri))
        if not isinstance(attachment, sbol3.Attachment):
            raise InputError(f"{component.display_id}: missing citation {uri}.")
        attachments.append(attachment)
    return attachments


def _citation_sources(document: sbol3.Document, component: sbol3.Component) -> tuple[str, ...]:
    attachments = [
        attachment
        for attachment in _component_attachments(document, component)
        if str(attachment.display_id) != "realization"
    ]
    attachments.sort(key=lambda attachment: str(attachment.display_id))
    return tuple(str(attachment.source) for attachment in attachments)


def _read_edit(document: sbol3.Document, component: sbol3.Component) -> Edit | None:
    activities = []
    for uri in component.generated_by:
        activity = document.find(str(uri))
        if isinstance(activity, sbol3.Activity) and any(
            str(kind).startswith(f"{_ACTIVITY_NAMESPACE}/") for kind in activity.types
        ):
            activities.append(activity)
    if not activities:
        return None
    if len(activities) != 1:
        raise InputError(f"{component.display_id}: expected one derivation.")
    activity = activities[0]
    activity_name = next(
        str(kind).rsplit("/", 1)[-1]
        for kind in activity.types
        if str(kind).startswith(f"{_ACTIVITY_NAMESPACE}/")
    )
    return _parse_derivation(activity_name, str(activity.description or ""))


def _module_from_feature(feature: sbol3.Feature, locus: int | None) -> Module:
    kind, antigen_class = _kind_of(feature)
    name = str(feature.name or feature.display_id)
    identity = str(feature.identity)
    module_id = str(feature.display_id)
    match kind:
        case ModuleKind.LEADER:
            return Leader(name, module_id, identity)
        case ModuleKind.TRAILER:
            return Trailer(name, module_id, identity)
        case ModuleKind.ANTIGEN:
            if antigen_class is None:
                raise InputError(f"Feature {feature.identity} has no antigen class.")
            return Antigen(module_id, name, antigen_class, locus, identity)
        case ModuleKind.REPORTER:
            return Reporter(module_id, name, locus, identity)
        case ModuleKind.VIRAL_GENE:
            return ViralGene(module_id, name, locus, identity)


def _architecture_from_component(
    document: sbol3.Document, component: sbol3.Component
) -> Architecture:
    replications = [
        _ROLE_REPLICATION[str(role)] for role in component.roles if str(role) in _ROLE_REPLICATION
    ]
    backbones = [
        str(role).rsplit("/", 1)[-1]
        for role in component.roles
        if str(role).startswith(f"{_BACKBONE_NAMESPACE}/")
    ]
    if len(replications) > 1 or len(backbones) != 1:
        raise InputError(
            f"{component.display_id}: an architecture needs one backbone and replication class."
        )
    replication = replications[0] if replications else Replication.UNKNOWN
    parent_unknown = any(str(role) == _ROLE_PARENT_UNKNOWN for role in component.roles)
    derived = [str(uri) for uri in component.derived_from]
    if len(derived) > 1 or (derived and parent_unknown):
        raise InputError(f"{component.display_id}: expected at most one parent.")
    parent: str | ParentStatus
    if derived:
        parent = derived[0].rstrip("/").rsplit("/", 1)[-1]
    elif parent_unknown:
        parent = ParentStatus.UNKNOWN
    else:
        parent = ParentStatus.ROOT
    ordered, free = _ordered_features(component)
    modules: list[Module] = []
    locus = 0
    trailers: list[sbol3.Feature] = []
    reached_trailer = False
    for feature in ordered:
        kind, _antigen = _kind_of(feature)
        if kind is ModuleKind.TRAILER:
            reached_trailer = True
            trailers.append(feature)
            continue
        if reached_trailer:
            raise InputError(f"{component.display_id}: the trailer is not last.")
        module_locus = None if kind is ModuleKind.LEADER else locus
        if module_locus is not None:
            locus += 1
        modules.append(_module_from_feature(feature, module_locus))
    for feature in sorted(free, key=lambda item: str(item.display_id)):
        kind, _antigen = _kind_of(feature)
        if kind in {ModuleKind.LEADER, ModuleKind.TRAILER}:
            raise InputError(f"{component.display_id}: the leader and trailer are ordered.")
        modules.append(_module_from_feature(feature, None))
    modules.extend(_module_from_feature(feature, None) for feature in trailers)
    return Architecture(
        str(component.display_id),
        str(component.name or component.display_id),
        backbones[0],
        replication,
        tuple(modules),
        parent=parent,
        citations=_citation_sources(document, component),
        edit=_read_edit(document, component),
    )


def architectures(document: sbol3.Document) -> list[Architecture]:
    """Every architecture component in ``document``, ordered by identity."""
    components = sorted(
        (obj for obj in document.objects if isinstance(obj, sbol3.Component)),
        key=lambda component: str(component.identity),
    )
    found = []
    for component in components:
        if any(str(role).startswith(f"{_BACKBONE_NAMESPACE}/") for role in component.roles):
            found.append(_architecture_from_component(document, component))
    return found


def architecture_of(document: sbol3.Document) -> Architecture:
    """The single architecture in ``document``."""
    found = architectures(document)
    if len(found) != 1:
        raise InputError(f"Expected one architecture, found {len(found)}.")
    return found[0]


_SEQUENCE_KEYS = frozenset(
    {
        "sequence",
        "sequences",
        "cds",
        "protein",
        "bases",
        "nucleotides",
        "amino_acids",
        "translation",
    }
)


def _local_id(identifier: str) -> str:
    cleaned = "".join(
        character if character.isalnum() or character == "_" else "_" for character in identifier
    )
    if not cleaned or cleaned[0].isdigit():
        cleaned = f"_{cleaned}"
    return cleaned


def _values(obj: sbol3.SBOLObject, predicate: str) -> list[str]:
    return [str(value) for value in obj._properties.get(predicate, [])]


def _set_index(obj: sbol3.TopLevel, index: int) -> None:
    sbol3.IntProperty(obj, _INDEX, 1, 1, initial_value=index)


def _index_of(obj: sbol3.SBOLObject) -> int:
    values = _values(obj, _INDEX)
    if len(values) != 1:
        raise InputError(f"{obj.display_id}: expected an index.")
    return int(values[0])


def _attachment_sources(document: sbol3.Document, owner: sbol3.TopLevel) -> tuple[str, ...]:
    sources = []
    for uri in owner.attachments:
        attachment = document.find(str(uri))
        if not isinstance(attachment, sbol3.Attachment):
            raise InputError(f"{owner.display_id}: missing attachment {uri}.")
        sources.append(str(attachment.source))
    return tuple(sources)


def _suffixed(values: list[str], namespace: str, label: str) -> str:
    matched = [value[len(namespace) + 1 :] for value in values if value.startswith(f"{namespace}/")]
    if len(matched) != 1:
        raise InputError(f"Expected one {label}.")
    return matched[0]


def _write_realization(
    document: sbol3.Document, component: sbol3.Component, realization: str
) -> None:
    attachment = sbol3.Attachment(
        "realization",
        f"{_ACCESSION}/{quote(realization, safe='')}",
        name=realization,
    )
    document.add(attachment)
    component.attachments.append(attachment.identity)


def _realization(document: sbol3.Document, component: sbol3.Component) -> str | None:
    found = [
        attachment
        for attachment in _component_attachments(document, component)
        if str(attachment.display_id) == "realization"
    ]
    if not found:
        return None
    if len(found) != 1 or not found[0].name:
        raise InputError(f"{component.display_id}: expected one realization.")
    return str(found[0].name)


def _model_source(document: sbol3.Document, component: sbol3.Component) -> str | None:
    models = list(component.models)
    if not models:
        return None
    if len(models) != 1:
        raise InputError(f"{component.display_id}: expected one model.")
    model = models[0]
    if not isinstance(model, sbol3.Model):
        found = document.find(str(model))
        if not isinstance(found, sbol3.Model):
            raise InputError(f"{component.display_id}: missing model.")
        model = found
    return str(model.source)


def _write_unmodeled(component: sbol3.Component, fields: tuple[Unmodeled, ...]) -> None:
    for field in fields:
        if field.key in _SEQUENCE_KEYS:
            raise InputError(f"sequence payload is not allowed ({field.key}).")
        token = quote(field.key, safe="")
        sbol3.TextProperty(component, f"{_UNMODELED_VALUE}{token}", 1, 1, initial_value=field.value)
        if field.source:
            sbol3.TextProperty(
                component, f"{_UNMODELED_SOURCE}{token}", 1, 1, initial_value=field.source
            )


def _unmodeled(component: sbol3.Component) -> tuple[Unmodeled, ...]:
    values: dict[str, str] = {}
    sources: dict[str, str] = {}
    for predicate, stored in component._properties.items():
        if predicate.startswith(_UNMODELED_VALUE):
            key = unquote(predicate.removeprefix(_UNMODELED_VALUE))
            if key in _SEQUENCE_KEYS:
                raise InputError(f"sequence payload is not allowed ({key}).")
            values[key] = str(stored[-1])
        elif predicate.startswith(_UNMODELED_SOURCE):
            sources[unquote(predicate.removeprefix(_UNMODELED_SOURCE))] = str(stored[-1])
    if set(sources).difference(values):
        raise InputError(f"{component.display_id}: an unmodeled source has no value.")
    return tuple(Unmodeled(key, values[key], sources.get(key, "")) for key in sorted(values))


def _write_implementation(
    document: sbol3.Document,
    component: sbol3.Component,
    architecture_id: str,
    implementation: Implementation,
    index: int,
) -> sbol3.Implementation:
    local = _local_id(implementation.id)
    stored = sbol3.Implementation(
        f"{architecture_id}__{local}",
        built=component,
        name=implementation.id,
        description=implementation.note or None,
    )
    sbol3.URIProperty(
        stored, sbol3.SBOL_ROLE, 1, 1, initial_value=f"{_STAGE}/{implementation.stage.value}"
    )
    _set_index(stored, index)
    document.add(stored)
    attachment = sbol3.Attachment(f"{architecture_id}__{local}_citation", implementation.citation)
    document.add(attachment)
    stored.attachments = [attachment]
    return stored


def _implementations(document: sbol3.Document, component: sbol3.Component) -> list[Implementation]:
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
                    str(obj.name), Stage(stage), citations[0], str(obj.description or "")
                ),
            )
        )
    found.sort(key=lambda item: item[0])
    return [implementation for _index, implementation in found]


def _write_observation(
    document: sbol3.Document,
    architecture_id: str,
    observation: Observation,
    implementation: sbol3.Implementation,
    index: int,
) -> None:
    local = _local_id(observation.id)
    activity = sbol3.Activity(
        f"{architecture_id}__{local}",
        types=[
            f"{_PHENOTYPE}/{observation.phenotype.value}",
            f"{_HOST}/{observation.host.value}",
            f"{_ASSAY}/{observation.assay.value}",
            f"{_OUTCOME}/{observation.outcome.value}",
        ],
        name=observation.id,
        description=observation.note or None,
        usage=[sbol3.Usage(implementation.identity)],
    )
    _set_index(activity, index)
    document.add(activity)
    data = sbol3.ExperimentalData(f"{architecture_id}__{local}_data", generated_by=[activity])
    document.add(data)
    attachment = sbol3.Attachment(f"{architecture_id}__{local}_citation", observation.citation)
    document.add(attachment)
    data.attachments = [attachment]


def _observations(
    document: sbol3.Document,
    implementations: list[Implementation],
    stored: list[sbol3.Implementation],
) -> list[Observation]:
    by_identity = {str(item.identity): item for item in stored}
    names = {item.id for item in implementations}
    found: list[tuple[int, Observation]] = []
    for data in document.objects:
        if not isinstance(data, sbol3.ExperimentalData):
            continue
        activities = [
            activity
            for uri in data.generated_by
            if isinstance(activity := document.find(str(uri)), sbol3.Activity)
        ]
        if len(activities) != 1:
            raise InputError(f"{data.display_id}: expected one generating activity.")
        activity = activities[0]
        if len(activity.usage) != 1 or not activity.name:
            raise InputError(f"{activity.display_id}: an observation names one implementation.")
        implementation = by_identity.get(str(activity.usage[0].entity))
        if implementation is None or str(implementation.name) not in names:
            raise InputError(
                f"{activity.name} cites unknown implementation {activity.usage[0].entity}."
            )
        terms = [str(kind) for kind in activity.types]
        citations = _attachment_sources(document, data)
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
                ),
            )
        )
    found.sort(key=lambda item: item[0])
    return [observation for _index, observation in found]


def design_to_document(design: Design) -> sbol3.Document:
    """Serialize a design, including evidence, as an SBOL3 document with no sequence."""
    document = to_document(design.architecture)
    component = document.find(f"{NAMESPACE}/{design.architecture.id}")
    if not isinstance(component, sbol3.Component):
        raise InputError(f"{design.architecture.id}: missing architecture component.")
    if design.aliases:
        sbol3.TextProperty(component, _ALIAS, 0, math.inf, initial_value=list(design.aliases))
    if design.realization:
        _write_realization(document, component, design.realization)
    if design.model_source:
        model = attach_model(document, design.architecture.id)
        model.source = design.model_source
    _write_unmodeled(component, design.unmodeled)
    stored: list[sbol3.Implementation] = []
    for index, implementation in enumerate(design.implementations):
        stored.append(
            _write_implementation(
                document, component, design.architecture.id, implementation, index
            )
        )
    by_name = {str(item.name): item for item in stored}
    for index, observation in enumerate(design.observations):
        built = by_name.get(observation.implementation_id)
        if built is None:
            raise InputError(
                f"{observation.id} cites unknown implementation {observation.implementation_id}."
            )
        _write_observation(document, design.architecture.id, observation, built, index)
    if any(isinstance(obj, sbol3.Sequence) for obj in document.objects):
        raise InputError("sequence payload is not allowed (sequence).")
    return document


def design_of(document: sbol3.Document) -> Design:
    """The single design in ``document``."""
    if any(isinstance(obj, sbol3.Sequence) for obj in document.objects):
        raise InputError("sequence payload is not allowed (sequence).")
    architecture = architecture_of(document)
    component = document.find(f"{NAMESPACE}/{architecture.id}")
    if not isinstance(component, sbol3.Component):
        raise InputError(f"{architecture.id}: missing architecture component.")
    implementations = _implementations(document, component)
    stored = [
        obj
        for obj in document.objects
        if isinstance(obj, sbol3.Implementation) and str(obj.built) == str(component.identity)
    ]
    return Design(
        architecture,
        tuple(_observations(document, implementations, stored)),
        tuple(implementations),
        _realization(document, component),
        _model_source(document, component),
        _unmodeled(component),
        tuple(sorted(_values(component, _ALIAS))),
    )
