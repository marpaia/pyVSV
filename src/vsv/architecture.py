"""Sequence-free VSV architectures.

A :class:`Genome` is a nucleotide sequence and its features. An :class:`Architecture` is the
design those features realize: backbone, replication class, and modules in genome order. Modules
are ordered with SBOL ``strictlyPrecedes`` constraints, so a design can be stored and compared
without a sequence. Projecting a canonical five-gene genome produces the chassis those designs
derive from.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import StrEnum
from itertools import pairwise

import sbol3

from vsv import NAMESPACE
from vsv.genes import Gene
from vsv.sbol import (
    ROLE_LEADER,
    ROLE_TRAILER,
    SO_GENE,
    SO_REGION,
    SO_VIRAL_SEQUENCE,
    Genome,
    display_id,
    new_document,
    pin_display_id,
    role_key,
)
from vsv.sequence import InputError

_ACTIVITY_NAMESPACE = f"{NAMESPACE}/activity"
_ANTIGEN_CLASS_NAMESPACE = f"{NAMESPACE}/antigen-class"
_BACKBONE_NAMESPACE = f"{NAMESPACE}/backbone"
_ROLE_NAMESPACE = f"{NAMESPACE}/role"

ROLE_REPORTER = f"{_ROLE_NAMESPACE}/reporter"
ROLE_REPLICATION_COMPETENT = f"{_ROLE_NAMESPACE}/replication-competent"
ROLE_REPLICATION_SINGLE_CYCLE = f"{_ROLE_NAMESPACE}/single-cycle"

_BACKBONE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")
_SLUG = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*\Z")
_ACTIVITIES = frozenset({"replace-module", "insert-module", "delete-module"})
_DERIVATION_KEYS = frozenset({"target", "module", "replication"})


class ModuleKind(StrEnum):
    """What a module contributes to the vector."""

    LEADER = "leader"
    TRAILER = "trailer"
    VIRAL_GENE = "viral-gene"
    ANTIGEN = "antigen"
    REPORTER = "reporter"


class Replication(StrEnum):
    """Whether the vector completes a cycle on its own."""

    COMPETENT = "replication-competent"
    SINGLE_CYCLE = "single-cycle"
    UNKNOWN = "unknown"


class LocusStatus(StrEnum):
    """Whether a module's transcription index is known."""

    KNOWN = "known"
    UNKNOWN = "unknown"
    NOT_TRANSCRIBED = "not-transcribed"


class ParentStatus(StrEnum):
    """Whether the parent design is the chassis root, named, or unpublished in the record."""

    ROOT = "root"
    KNOWN = "known"
    UNKNOWN = "unknown"


_REPLICATION_ROLE = {
    Replication.COMPETENT: ROLE_REPLICATION_COMPETENT,
    Replication.SINGLE_CYCLE: ROLE_REPLICATION_SINGLE_CYCLE,
}
_ROLE_REPLICATION = {role: replication for replication, role in _REPLICATION_ROLE.items()}
_TERMINAL = frozenset({ModuleKind.LEADER, ModuleKind.TRAILER})


def _require_display_id(identifier: str, what: str) -> None:
    if not identifier or display_id(identifier) != identifier:
        raise InputError(f"Invalid {what} {identifier!r}.")


@dataclass(frozen=True, slots=True)
class Module:
    """One element of an architecture, in genome order.

    ``locus`` is the zero-based transcription index when it is known. An unknown locus is a
    transcribed module whose paper does not state a position. The leader and the trailer are not
    transcribed. ``antigen_class`` is a slug such as ``filovirus-glycoprotein``, present only on
    antigen modules.
    """

    id: str
    kind: ModuleKind
    name: str
    locus: int | None = None
    antigen_class: str | None = None
    feature_identity: str | None = None
    locus_status: LocusStatus | None = None

    def __post_init__(self) -> None:
        _require_display_id(self.id, "module id")
        if not self.name:
            raise InputError(f"Module {self.id} needs a name.")
        if self.kind is ModuleKind.ANTIGEN:
            if self.antigen_class is None or _SLUG.fullmatch(self.antigen_class) is None:
                raise InputError(f"Module {self.id} needs an antigen class slug.")
        elif self.antigen_class is not None:
            raise InputError(f"Module {self.id} is not an antigen and has no antigen class.")
        if self.kind in _TERMINAL:
            terminal = self.locus_status in {None, LocusStatus.NOT_TRANSCRIBED}
            if self.locus is not None or not terminal:
                raise InputError(f"Module {self.id} is not transcribed and has no locus.")
            status = LocusStatus.NOT_TRANSCRIBED
        elif self.locus is None or self.locus_status is LocusStatus.UNKNOWN:
            if self.locus is not None:
                raise InputError(f"Module {self.id}: an unknown locus has no index.")
            status = LocusStatus.UNKNOWN
        elif self.locus < 0:
            raise InputError(f"Module {self.id} has a negative locus.")
        else:
            status = LocusStatus.KNOWN
        object.__setattr__(self, "locus_status", status)


@dataclass(frozen=True, slots=True)
class Derivation:
    """The transformation that produced an architecture from its parent.

    ``target`` is the module that was replaced or deleted, or the module an insertion follows.
    ``module_id`` is the module inserted or written in its place. ``replication`` is set when the
    transformation also changes replication class.
    """

    activity: str
    target: str
    module_id: str | None = None
    replication: Replication | None = None

    def __post_init__(self) -> None:
        if self.activity not in _ACTIVITIES:
            raise InputError(f"Unknown derivation {self.activity!r}.")
        _require_display_id(self.target, "derivation target")
        if self.module_id is not None:
            _require_display_id(self.module_id, "derivation module")
        if self.activity == "delete-module" and self.module_id is not None:
            raise InputError("A deletion derivation has no module.")
        if self.activity != "delete-module" and self.module_id is None:
            raise InputError(f"A {self.activity} derivation names the resulting module.")


@dataclass(frozen=True, slots=True)
class Architecture:
    """A VSV vector design: ordered modules, replication class, and lineage."""

    id: str
    name: str
    backbone: str
    replication: Replication
    modules: tuple[Module, ...]
    parent_id: str | None = None
    parent_status: ParentStatus = ParentStatus.ROOT
    citations: tuple[str, ...] = ()
    derivation: Derivation | None = None

    def __post_init__(self) -> None:
        _require_display_id(self.id, "architecture id")
        if not self.name:
            raise InputError(f"Architecture {self.id} needs a name.")
        if _BACKBONE.fullmatch(self.backbone) is None:
            raise InputError(f"Invalid backbone {self.backbone!r}.")
        if self.parent_id is not None:
            _require_display_id(self.parent_id, "parent id")
            status = ParentStatus.KNOWN
        elif self.parent_status is ParentStatus.UNKNOWN:
            status = ParentStatus.UNKNOWN
        else:
            status = ParentStatus.ROOT
        object.__setattr__(self, "parent_status", status)
        object.__setattr__(self, "modules", tuple(self.modules))
        object.__setattr__(self, "citations", tuple(self.citations))
        ids = [module.id for module in self.modules]
        if len(ids) != len(set(ids)):
            raise InputError(f"{self.id}: module ids must be unique.")


def canonical_chassis(
    *,
    identifier: str = "vsv_indiana",
    name: str = "VSV Indiana",
    backbone: str = "VSV-Indiana",
) -> Architecture:
    """The canonical VSV chassis: leader, N, P, M, G, L, and trailer."""
    modules: list[Module] = [Module("leader", ModuleKind.LEADER, "leader")]
    modules.extend(
        Module(gene, ModuleKind.VIRAL_GENE, gene, locus=locus) for locus, gene in enumerate(Gene)
    )
    modules.append(Module("trailer", ModuleKind.TRAILER, "trailer"))
    return Architecture(identifier, name, backbone, Replication.COMPETENT, tuple(modules))


def _canonical_features(genome: Genome) -> list[sbol3.SequenceFeature]:
    """Leader, the five genes, and trailer, for a genome that has every coding sequence."""
    found = genome.genes
    missing = [gene for gene in Gene if gene not in found or found[gene].cds_feature is None]
    if missing:
        names = ", ".join(missing)
        raise InputError(f"{genome.id}: canonical chassis requires coding sequences for {names}.")
    selected: list[sbol3.SequenceFeature] = []
    leaders = genome.features(ROLE_LEADER)
    trailers = genome.features(ROLE_TRAILER)
    if len(leaders) > 1 or len(trailers) > 1:
        raise InputError(f"{genome.id}: expected at most one leader and one trailer.")
    selected.extend(leaders)
    for gene in Gene:
        feature = found[gene].gene_feature or found[gene].cds_feature
        if feature is None:
            raise InputError(f"{genome.id}: gene {gene} has no feature.")
        selected.append(feature)
    selected.extend(trailers)
    return selected


def project(genome: Genome, *, backbone: str = "VSV") -> Architecture:
    """The canonical architecture realized by an identified five-gene genome.

    Genes must already be identified, as :func:`vsv.load` does. A genome that lacks one of
    N, P, M, G, or L, or a coding sequence for one of them, is not a canonical chassis.
    """
    _canonical_features(genome)
    found = genome.genes
    modules: list[Module] = []
    leaders = genome.features(ROLE_LEADER)
    if leaders:
        leader = leaders[0]
        modules.append(
            Module(
                "leader",
                ModuleKind.LEADER,
                str(leader.name or "leader"),
                feature_identity=str(leader.identity),
            )
        )
    for locus, gene in enumerate(Gene):
        feature = found[gene].gene_feature or found[gene].cds_feature
        if feature is None:
            raise InputError(f"{genome.id}: gene {gene} has no feature.")
        modules.append(
            Module(
                gene,
                ModuleKind.VIRAL_GENE,
                gene,
                locus=locus,
                feature_identity=str(feature.identity),
            )
        )
    trailers = genome.features(ROLE_TRAILER)
    if trailers:
        trailer = trailers[0]
        modules.append(
            Module(
                "trailer",
                ModuleKind.TRAILER,
                str(trailer.name or "trailer"),
                feature_identity=str(trailer.identity),
            )
        )
    return Architecture(
        display_id(genome.id),
        genome.id,
        backbone,
        Replication.COMPETENT,
        tuple(modules),
    )


def add_order_constraints(genome: Genome) -> None:
    """Record genome order as ``strictlyPrecedes`` constraints between canonical features.

    A second call adds nothing. The genome must already identify N, P, M, G, and L.
    """
    features = _canonical_features(genome)
    existing = {
        (str(constraint.subject), str(constraint.object), str(constraint.restriction))
        for constraint in genome.component.constraints
    }
    for earlier, later in pairwise(features):
        key = (str(earlier.identity), str(later.identity), sbol3.SBOL_STRICTLY_PRECEDES)
        if key in existing:
            continue
        genome.component.constraints.append(
            sbol3.Constraint(sbol3.SBOL_STRICTLY_PRECEDES, earlier, later)
        )
        existing.add(key)


def _module_roles(module: Module) -> list[str]:
    if module.kind is ModuleKind.LEADER:
        return [SO_REGION, ROLE_LEADER]
    if module.kind is ModuleKind.TRAILER:
        return [SO_REGION, ROLE_TRAILER]
    if module.kind is ModuleKind.ANTIGEN:
        return [SO_GENE, f"{_ANTIGEN_CLASS_NAMESPACE}/{module.antigen_class}"]
    if module.kind is ModuleKind.REPORTER:
        return [SO_GENE, ROLE_REPORTER]
    return [SO_GENE]


def _require_genome_order(architecture: Architecture) -> None:
    """Refuse to serialize a module list whose loci disagree with its order.

    Constraints record order only. Writing a list that is not already in locus order would
    assign new loci on the way back in.
    """
    body = [module for module in architecture.modules if module.kind not in _TERMINAL]
    if any(module.locus is None for module in body):
        raise InputError(f"{architecture.id}: every transcribed module needs a locus.")
    if [module.locus for module in body] != list(range(len(body))):
        raise InputError(f"{architecture.id}: modules are not in locus order.")
    leaders = [module for module in architecture.modules if module.kind is ModuleKind.LEADER]
    trailers = [module for module in architecture.modules if module.kind is ModuleKind.TRAILER]
    ordered = leaders + body + trailers
    if [module.id for module in architecture.modules] != [module.id for module in ordered]:
        raise InputError(f"{architecture.id}: modules are not in genome order.")


def _derivation_text(derivation: Derivation) -> str:
    parts = [f"target={derivation.target}"]
    if derivation.module_id is not None:
        parts.append(f"module={derivation.module_id}")
    if derivation.replication is not None:
        parts.append(f"replication={derivation.replication.value}")
    return ";".join(parts)


def _parse_derivation(activity: str, description: str) -> Derivation:
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
    return Derivation(activity, fields["target"], fields.get("module"), replication)


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
    if architecture.derivation is not None:
        if architecture.parent_id is None:
            raise InputError(f"{architecture.id}: a derivation requires a parent.")
        activity = sbol3.Activity(
            f"{architecture.id}_derivation",
            types=[f"{_ACTIVITY_NAMESPACE}/{architecture.derivation.activity}"],
            description=_derivation_text(architecture.derivation),
        )
        usage = sbol3.Usage(f"{NAMESPACE}/{architecture.parent_id}")
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
    if role_key(ROLE_REPORTER) in keys:
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


def _read_derivation(document: sbol3.Document, component: sbol3.Component) -> Derivation | None:
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
    parent_id = None if not derived else derived[0].rstrip("/").rsplit("/", 1)[-1]
    modules: list[Module] = []
    locus = 0
    for feature in _ordered_features(component):
        kind, antigen_class = _kind_of(feature)
        module_locus = None if kind in _TERMINAL else locus
        if module_locus is not None:
            locus += 1
        modules.append(
            Module(
                str(feature.display_id),
                kind,
                str(feature.name or feature.display_id),
                locus=module_locus,
                antigen_class=antigen_class,
                feature_identity=str(feature.identity),
            )
        )
    return Architecture(
        str(component.display_id),
        str(component.name or component.display_id),
        backbones[0],
        replication,
        tuple(modules),
        parent_id=parent_id,
        citations=_citation_sources(document, component),
        derivation=_read_derivation(document, component),
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
