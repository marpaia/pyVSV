"""The recorded edit that produced an architecture from its parent.

An edit names modules by display id. The module bodies live on the architecture. Applying an
edit is :meth:`vsv.design.architecture.Architecture.replace` and its siblings.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from vsv.design.modules import require_display_id

if TYPE_CHECKING:
    from vsv.design.architecture import Replication


@dataclass(frozen=True, slots=True)
class Replace:
    """``target`` was replaced by the module ``module``."""

    target: str
    module: str
    replication: Replication | None = None

    def __post_init__(self) -> None:
        require_display_id(self.target, "derivation target")
        require_display_id(self.module, "derivation module")


@dataclass(frozen=True, slots=True)
class Insert:
    """``module`` was inserted after ``after``."""

    after: str
    module: str
    replication: Replication | None = None

    def __post_init__(self) -> None:
        require_display_id(self.after, "derivation target")
        require_display_id(self.module, "derivation module")


@dataclass(frozen=True, slots=True)
class Delete:
    """``target`` was deleted."""

    target: str
    replication: Replication | None = None

    def __post_init__(self) -> None:
        require_display_id(self.target, "derivation target")


Edit = Replace | Insert | Delete
