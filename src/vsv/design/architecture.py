"""Sequence-free VSV architectures.

A genome is a nucleotide sequence and its features. An architecture is the design those features
realize: backbone, replication class, and modules in genome order. Projecting a canonical
five-gene genome produces the chassis other designs derive from.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import StrEnum
from itertools import pairwise

import sbol3

from vsv.design.edits import Delete, Edit, Insert, Replace
from vsv.design.modules import (
    Leader,
    Module,
    Trailer,
    Transcribed,
    ViralGene,
    at_locus,
    require_display_id,
)
from vsv.errors import InputError
from vsv.genes import Gene
from vsv.genome.sbol import (
    ROLE_LEADER,
    ROLE_TRAILER,
    Genome,
    Part,
    display_id,
    role_key,
)

_BACKBONE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")


class Replication(StrEnum):
    """Whether the vector completes a cycle on its own."""

    COMPETENT = "replication-competent"
    SINGLE_CYCLE = "single-cycle"
    UNKNOWN = "unknown"


class ParentStatus(StrEnum):
    """A parent that is not a named architecture.

    A string parent is the parent's id. Root and unknown are the only other states.
    """

    ROOT = "root"
    UNKNOWN = "unknown"


def _reindex(modules: list[Module]) -> tuple[Module, ...]:
    leaders = [module for module in modules if isinstance(module, Leader)]
    trailers = [module for module in modules if isinstance(module, Trailer)]
    body = [module for module in modules if isinstance(module, Transcribed)]
    return tuple(
        [at_locus(module, None) for module in leaders]
        + [at_locus(module, index) for index, module in enumerate(body)]
        + [at_locus(module, None) for module in trailers]
    )


@dataclass(frozen=True, slots=True)
class Architecture:
    """A VSV vector design: ordered modules, replication class, and lineage."""

    id: str
    name: str
    backbone: str
    replication: Replication
    modules: tuple[Module, ...]
    parent: str | ParentStatus = ParentStatus.ROOT
    citations: tuple[str, ...] = ()
    edit: Edit | None = None

    def __post_init__(self) -> None:
        require_display_id(self.id, "architecture id")
        if not self.name:
            raise InputError(f"Architecture {self.id} needs a name.")
        if _BACKBONE.fullmatch(self.backbone) is None:
            raise InputError(f"Invalid backbone {self.backbone!r}.")
        if isinstance(self.parent, ParentStatus):
            if self.edit is not None:
                raise InputError(f"{self.id}: a derivation requires a parent.")
        else:
            require_display_id(self.parent, "parent id")
        object.__setattr__(self, "modules", tuple(self.modules))
        object.__setattr__(self, "citations", tuple(self.citations))
        identifiers = [module.id for module in self.modules]
        if len(identifiers) != len(set(identifiers)):
            raise InputError(f"{self.id}: module ids must be unique.")

    @property
    def parent_id(self) -> str | None:
        """The parent architecture's id, when the parent is named."""
        if isinstance(self.parent, ParentStatus):
            return None
        return self.parent

    def module(self, module_id: str) -> Module:
        """The module with ``module_id``."""
        found = next((module for module in self.modules if module.id == module_id), None)
        if found is None:
            raise InputError(f"{self.id}: no module {module_id}.")
        return found

    @property
    def transcribed(self) -> tuple[Transcribed, ...]:
        """Modules that are transcribed, in genome order."""
        return tuple(module for module in self.modules if isinstance(module, Transcribed))

    @property
    def antigen_classes(self) -> tuple[str, ...]:
        """Antigen classes in module order, without duplicates."""
        classes: list[str] = []
        for module in self.modules:
            antigen_class = module.antigen_class
            if antigen_class is not None and antigen_class not in classes:
                classes.append(antigen_class)
        return tuple(classes)

    def replace(
        self,
        target: str,
        module: Transcribed,
        identifier: str,
        *,
        name: str | None = None,
        replication: Replication | None = None,
    ) -> Architecture:
        """Replace one transcribed module. The new module keeps the target's position."""
        current = self.module(target)
        if not isinstance(current, Transcribed) or not isinstance(module, Transcribed):
            raise InputError(f"{self.id}: leader and trailer cannot be replaced.")
        if any(item.id == module.id and item.id != current.id for item in self.modules):
            raise InputError(f"{self.id}: module {module.id} already exists.")
        modules = [module if item.id == current.id else item for item in self.modules]
        return self._child(
            identifier, name, replication, modules, Replace(target, module.id, replication)
        )

    def insert(
        self,
        after: str,
        module: Transcribed,
        identifier: str,
        *,
        name: str | None = None,
        replication: Replication | None = None,
    ) -> Architecture:
        """Insert a transcribed module after ``after`` and renumber later loci."""
        current = self.module(after)
        if isinstance(current, Trailer) or not isinstance(module, Transcribed):
            raise InputError(f"{self.id}: insertions go after the leader or a transcribed module.")
        if any(item.id == module.id for item in self.modules):
            raise InputError(f"{self.id}: module {module.id} already exists.")
        modules = list(self.modules)
        index = next(position for position, item in enumerate(modules) if item.id == current.id)
        modules.insert(index + 1, module)
        return self._child(
            identifier, name, replication, modules, Insert(after, module.id, replication)
        )

    def delete(
        self,
        target: str,
        identifier: str,
        *,
        name: str | None = None,
        replication: Replication | None = None,
    ) -> Architecture:
        """Delete one transcribed module and renumber the modules that follow it."""
        current = self.module(target)
        if not isinstance(current, Transcribed):
            raise InputError(f"{self.id}: leader and trailer cannot be deleted.")
        if len(self.transcribed) == 1:
            raise InputError(f"{self.id}: refusing to delete the last transcribed module.")
        modules = [item for item in self.modules if item.id != current.id]
        return self._child(identifier, name, replication, modules, Delete(target, replication))

    def _child(
        self,
        identifier: str,
        name: str | None,
        replication: Replication | None,
        modules: list[Module],
        edit: Edit,
    ) -> Architecture:
        require_display_id(identifier, "architecture id")
        if identifier == self.id:
            raise InputError(f"{self.id}: a derived architecture needs its own id.")
        return Architecture(
            identifier,
            self.name if name is None else name,
            self.backbone,
            self.replication if replication is None else replication,
            _reindex(modules),
            parent=self.id,
            citations=self.citations,
            edit=edit,
        )


def canonical_chassis(
    *,
    identifier: str = "vsv_indiana",
    name: str = "VSV Indiana",
    backbone: str = "VSV-Indiana",
) -> Architecture:
    """The canonical VSV chassis: leader, N, P, M, G, L, and trailer."""
    modules: list[Module] = [Leader()]
    modules.extend(ViralGene(gene, gene, locus=locus) for locus, gene in enumerate(Gene))
    modules.append(Trailer())
    return Architecture(identifier, name, backbone, Replication.COMPETENT, tuple(modules))


def _canonical_features(genome: Genome) -> list[Part]:
    """Leader, the five genes, and trailer, for a genome that has every coding sequence."""
    found = genome.genes
    missing = [gene for gene in Gene if gene not in found or found[gene].cds_feature is None]
    if missing:
        names = ", ".join(missing)
        raise InputError(f"{genome.id}: canonical chassis requires coding sequences for {names}.")
    selected: list[Part] = []
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
    features = _canonical_features(genome)
    modules: list[Module] = []
    leaders = [feature for feature in features if role_key(ROLE_LEADER) in _keys(feature)]
    trailers = [feature for feature in features if role_key(ROLE_TRAILER) in _keys(feature)]
    genes = [feature for feature in features if feature not in leaders and feature not in trailers]
    if leaders:
        leader = leaders[0]
        modules.append(
            Leader(name=str(leader.name or "leader"), feature_identity=str(leader.identity))
        )
    for locus, (gene, feature) in enumerate(zip(Gene, genes, strict=True)):
        modules.append(ViralGene(gene, gene, locus=locus, feature_identity=str(feature.identity)))
    if trailers:
        trailer = trailers[0]
        modules.append(
            Trailer(name=str(trailer.name or "trailer"), feature_identity=str(trailer.identity))
        )
    return Architecture(
        display_id(genome.id),
        genome.id,
        backbone,
        Replication.COMPETENT,
        tuple(modules),
    )


def _keys(feature: Part) -> set[str]:
    return {role_key(str(role)) for role in feature.roles}


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
