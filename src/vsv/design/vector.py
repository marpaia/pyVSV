"""Replication of a VSV vector and the form of its glycoprotein."""

from __future__ import annotations

from enum import StrEnum


class Replication(StrEnum):
    """Whether the vector completes a cycle on its own."""

    COMPETENT = "replication-competent"
    SINGLE_CYCLE = "single-cycle"
    UNKNOWN = "unknown"


class GlycoproteinForm(StrEnum):
    """How much of a glycoprotein the insert encodes."""

    FULL_LENGTH = "full-length"
    SOLUBLE = "soluble"
    MUCIN_DELETED = "mucin-deleted"
    TAIL_CHIMERA = "tail-chimera"
    UNSPECIFIED = "unspecified"

    def supplies_entry(self) -> bool:
        """Whether this form can be the entry protein of a competent genome."""
        return self in {self.FULL_LENGTH, self.TAIL_CHIMERA, self.MUCIN_DELETED}

    def needs_anchor(self) -> bool:
        """Whether a competent replacement must show a hydrophobic anchor."""
        return self in {self.FULL_LENGTH, self.TAIL_CHIMERA}
