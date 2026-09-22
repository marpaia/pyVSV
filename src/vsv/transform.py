"""Transformations of a VSV architecture.

Each transformation returns a new architecture and leaves its input unchanged. The result records
its parent and the edit that produced it. Replication class changes only when the caller passes
one; deleting G does not by itself make a vector single-cycle.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from vsv.architecture import (
    Architecture,
    Derivation,
    LocusStatus,
    Module,
    ModuleKind,
    Replication,
)
from vsv.sbol import display_id
from vsv.sequence import InputError

_TERMINAL = frozenset({ModuleKind.LEADER, ModuleKind.TRAILER})


def _reindex(modules: list[Module]) -> tuple[Module, ...]:
    leaders = [module for module in modules if module.kind is ModuleKind.LEADER]
    trailers = [module for module in modules if module.kind is ModuleKind.TRAILER]
    body = [module for module in modules if module.kind not in _TERMINAL]
    return tuple(
        [
            replace(
                module, locus=None, locus_status=LocusStatus.NOT_TRANSCRIBED, feature_identity=None
            )
            for module in leaders
        ]
        + [
            replace(module, locus=index, locus_status=LocusStatus.KNOWN, feature_identity=None)
            for index, module in enumerate(body)
        ]
        + [
            replace(
                module, locus=None, locus_status=LocusStatus.NOT_TRANSCRIBED, feature_identity=None
            )
            for module in trailers
        ]
    )


def _require(architecture: Architecture, module_id: str) -> Module:
    match = next((module for module in architecture.modules if module.id == module_id), None)
    if match is None:
        raise InputError(f"{architecture.id}: no module {module_id}.")
    return match


def _derive(
    source: Architecture,
    identifier: str,
    name: str | None,
    replication: Replication | None,
    modules: list[Module],
    derivation: Derivation,
) -> Architecture:
    if display_id(identifier) != identifier:
        raise InputError(f"Invalid architecture id {identifier!r}.")
    if identifier == source.id:
        raise InputError(f"{source.id}: a derived architecture needs its own id.")
    return Architecture(
        identifier,
        source.name if name is None else name,
        source.backbone,
        source.replication if replication is None else replication,
        _reindex(modules),
        parent_id=source.id,
        citations=source.citations,
        derivation=derivation,
    )


@dataclass(frozen=True, slots=True)
class ReplaceModule:
    """Replace one transcribed module.

    The new module keeps the target's position. Its own locus is ignored and reassigned.
    """

    target_id: str
    module: Module
    identifier: str
    name: str | None = None
    replication: Replication | None = None

    def apply(self, source: Architecture) -> Architecture:
        target = _require(source, self.target_id)
        if target.kind in _TERMINAL or self.module.kind in _TERMINAL:
            raise InputError(f"{source.id}: leader and trailer cannot be replaced.")
        if any(module.id == self.module.id and module.id != target.id for module in source.modules):
            raise InputError(f"{source.id}: module {self.module.id} already exists.")
        modules = [self.module if module.id == target.id else module for module in source.modules]
        return _derive(
            source,
            self.identifier,
            self.name,
            self.replication,
            modules,
            Derivation("replace-module", self.target_id, self.module.id, self.replication),
        )


@dataclass(frozen=True, slots=True)
class InsertModule:
    """Insert a transcribed module after ``after_id`` and renumber later loci."""

    after_id: str
    module: Module
    identifier: str
    name: str | None = None
    replication: Replication | None = None

    def apply(self, source: Architecture) -> Architecture:
        after = _require(source, self.after_id)
        if after.kind is ModuleKind.TRAILER or self.module.kind in _TERMINAL:
            raise InputError(
                f"{source.id}: insertions go after the leader or a transcribed module."
            )
        if any(module.id == self.module.id for module in source.modules):
            raise InputError(f"{source.id}: module {self.module.id} already exists.")
        modules = list(source.modules)
        index = next(position for position, module in enumerate(modules) if module.id == after.id)
        modules.insert(index + 1, self.module)
        return _derive(
            source,
            self.identifier,
            self.name,
            self.replication,
            modules,
            Derivation("insert-module", self.after_id, self.module.id, self.replication),
        )


@dataclass(frozen=True, slots=True)
class DeleteModule:
    """Delete one transcribed module and renumber the modules that follow it."""

    target_id: str
    identifier: str
    name: str | None = None
    replication: Replication | None = None

    def apply(self, source: Architecture) -> Architecture:
        target = _require(source, self.target_id)
        if target.kind in _TERMINAL:
            raise InputError(f"{source.id}: leader and trailer cannot be deleted.")
        transcribed = [module for module in source.modules if module.kind not in _TERMINAL]
        if len(transcribed) == 1:
            raise InputError(f"{source.id}: refusing to delete the last transcribed module.")
        modules = [module for module in source.modules if module.id != target.id]
        return _derive(
            source,
            self.identifier,
            self.name,
            self.replication,
            modules,
            Derivation("delete-module", self.target_id, replication=self.replication),
        )
