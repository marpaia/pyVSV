"""An architecture as an SBOL3 document with no nucleotide sequence.

Module order is a chain of ``strictlyPrecedes`` constraints. The edit that produced a design is
an Activity whose type is ``replace-module``, ``insert-module``, or ``delete-module`` and whose
description is ``target=…;module=…;replication=…``.
"""

from __future__ import annotations

from itertools import pairwise

import sbol3

from vsv import NAMESPACE
from vsv.design.architecture import Architecture, ParentStatus, Replication
from vsv.design.edits import Delete, Edit, Insert, Replace
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


def _require_genome_order(architecture: Architecture) -> None:
    """Refuse to serialize a module list whose loci disagree with its order.

    Constraints record order only. Writing a list that is not already in locus order would
    assign new loci on the way back in.
    """
    body = list(architecture.transcribed)
    if any(module.locus is None for module in body):
        raise InputError(f"{architecture.id}: every transcribed module needs a locus.")
    if [module.locus for module in body] != list(range(len(body))):
        raise InputError(f"{architecture.id}: modules are not in locus order.")
    leaders = [module for module in architecture.modules if isinstance(module, Leader)]
    trailers = [module for module in architecture.modules if isinstance(module, Trailer)]
    ordered = leaders + body + trailers
    if [module.id for module in architecture.modules] != [module.id for module in ordered]:
        raise InputError(f"{architecture.id}: modules are not in genome order.")


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
    if architecture.replication is Replication.UNKNOWN:
        raise InputError(f"{architecture.id}: replication is unknown.")
    _require_genome_order(architecture)
    document = new_document()
    component = sbol3.Component(
        architecture.id,
        [sbol3.SBO_DNA],
        roles=[
            SO_VIRAL_SEQUENCE,
            _REPLICATION_ROLE[architecture.replication],
            f"{_BACKBONE_NAMESPACE}/{architecture.backbone}",
        ],
        name=architecture.name,
    )
    if architecture.parent_id is not None:
        component.derived_from = [f"{NAMESPACE}/{architecture.parent_id}"]
    document.add(component)
    features: list[sbol3.LocalSubComponent] = []
    for module in architecture.modules:
        feature = sbol3.LocalSubComponent(
            [sbol3.SBO_DNA], roles=_module_roles(module), name=module.name
        )
        pin_display_id(feature, module.id)
        component.features.append(feature)
        features.append(feature)
    for earlier, later in pairwise(features):
        component.constraints.append(sbol3.Constraint(sbol3.SBOL_STRICTLY_PRECEDES, earlier, later))
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


def _ordered_features(component: sbol3.Component) -> list[sbol3.Feature]:
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
    if not features:
        return []
    starts = [feature for feature in features if str(feature.identity) not in incoming]
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
    if len(ordered) != len(features):
        raise InputError(f"{component.display_id}: module order does not include every module.")
    return ordered


def _citation_sources(document: sbol3.Document, component: sbol3.Component) -> tuple[str, ...]:
    attachments = []
    for uri in component.attachments:
        attachment = document.find(str(uri))
        if not isinstance(attachment, sbol3.Attachment):
            raise InputError(f"{component.display_id}: missing citation {uri}.")
        attachments.append(attachment)
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
    replication = next(
        (
            _ROLE_REPLICATION[str(role)]
            for role in component.roles
            if str(role) in _ROLE_REPLICATION
        ),
        None,
    )
    backbones = [
        str(role).rsplit("/", 1)[-1]
        for role in component.roles
        if str(role).startswith(f"{_BACKBONE_NAMESPACE}/")
    ]
    if replication is None or len(backbones) != 1:
        raise InputError(
            f"{component.display_id}: an architecture needs one backbone and replication class."
        )
    derived = [str(uri) for uri in component.derived_from]
    if len(derived) > 1:
        raise InputError(f"{component.display_id}: expected at most one parent.")
    parent: str | ParentStatus
    parent = ParentStatus.ROOT if not derived else derived[0].rstrip("/").rsplit("/", 1)[-1]
    modules: list[Module] = []
    locus = 0
    for feature in _ordered_features(component):
        kind, _antigen = _kind_of(feature)
        module_locus = None if kind in {ModuleKind.LEADER, ModuleKind.TRAILER} else locus
        if module_locus is not None:
            locus += 1
        modules.append(_module_from_feature(feature, module_locus))
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
