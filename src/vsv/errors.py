"""Errors raised when input cannot be read."""

from __future__ import annotations


class InputError(ValueError):
    """Raised when input cannot be read as a genome or a design."""
