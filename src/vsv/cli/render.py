"""Text layout shared by the VSV commands."""

from __future__ import annotations

import textwrap
from collections.abc import Sequence

from vsv.genes import Gene


def cli_name(text: str) -> str:
    """Show a VSV gene as its one-letter CLI name."""
    for gene in sorted(Gene, key=lambda item: len(item.value), reverse=True):
        if text == gene.value:
            return gene.short
        prefix = f"{gene.value} "
        if text.startswith(prefix):
            return f"{gene.short}{text[len(gene.value) :]}"
    return text


def cli_text(text: str) -> str:
    """Replace gene names in ``text`` with their one-letter CLI names."""
    for gene in sorted(Gene, key=lambda item: len(item.value), reverse=True):
        text = text.replace(gene.value, gene.short)
    return text


def row(cells: Sequence[str], widths: Sequence[int]) -> str:
    return "  ".join(cell.ljust(width) for cell, width in zip(cells, widths, strict=True)).rstrip()


def widths(header: Sequence[str], rows: Sequence[Sequence[str]]) -> list[int]:
    return [max(len(cell) for cell in column) for column in zip(header, *rows, strict=True)]


def table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    column_widths = widths(header, rows)
    return "\n".join(row(line, column_widths) for line in (header, *rows))


def indent(text: str, spaces: int = 2) -> str:
    pad = " " * spaces
    return "\n".join(pad + line if line else "" for line in text.splitlines())


def wrap(text: str, prefix: str, width: int = 78) -> list[str]:
    return textwrap.wrap(text, width=width, initial_indent=prefix, subsequent_indent=prefix)
