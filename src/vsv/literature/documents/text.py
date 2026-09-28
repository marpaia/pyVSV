"""Plain-text cleanup and plain-text block parsing."""

from __future__ import annotations

import html
import re

from vsv.ingest import redact
from vsv.literature.records import TextBlock

_HEADING = re.compile(r".{1,80}[^\.?!]\Z")
_COLUMNS = re.compile(r"\s{2,}|\t")


def clean_text(value: str) -> str:
    """Drop markup and long nucleotide runs from text that will be stored or quoted."""
    stripped = re.sub(r"<[^>]+>", " ", html.unescape(value))
    collapsed = re.sub(r"\s+", " ", stripped).strip()
    return redact(collapsed)


def section_text(text: str) -> tuple[TextBlock, ...]:
    """Blocks from plain article text.

    A short line that does not end a sentence becomes the section name. A run of
    lines with two or more columns becomes one table row per line. Other paragraphs
    stay paragraphs.
    """
    found: list[tuple[str, str, str]] = []
    section = "text"
    for chunk in re.split(r"\n\s*\n", text):
        lines = [line.strip() for line in chunk.splitlines() if line.strip()]
        if not lines:
            continue
        if len(lines) == 1 and _HEADING.fullmatch(lines[0]):
            section = lines[0]
            continue
        if len(lines) > 1 and all(_COLUMNS.search(line) for line in lines):
            for line in lines:
                cells = [cell.strip() for cell in _COLUMNS.split(line) if cell.strip()]
                found.append((section, "table-row", " | ".join(cells)))
            continue
        found.append((section, "paragraph", " ".join(lines)))
    return text_blocks(found)


def text_blocks(found: list[tuple[str, str, str]]) -> tuple[TextBlock, ...]:
    """Quotable blocks from section, kind, and text. Nucleotide runs are removed."""
    blocks: list[TextBlock] = []
    for section, kind, text in found:
        cleaned = redact(" ".join(text.split()))
        if not cleaned:
            continue
        blocks.append(TextBlock(f"b{len(blocks) + 1:04d}", section or "article", kind, cleaned))
    return tuple(blocks)
