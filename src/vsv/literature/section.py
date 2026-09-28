"""Split an article into the blocks a reading is allowed to quote."""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from vsv.errors import InputError
from vsv.ingest import redact
from vsv.literature.records import Block
from vsv.literature.xmlbody import local_name

_SKIP = frozenset({"ref-list", "ack", "funding-group", "author-notes", "fn-group"})
_ACK = frozenset({"acknowledgements", "acknowledgments", "funding"})
_HEADING = re.compile(r".{1,80}[^\.?!]\Z")
_COLUMNS = re.compile(r"\s{2,}|\t")


def section_xml(xml: str) -> tuple[Block, ...]:
    """Blocks from a JATS article, in document order.

    Reference lists, acknowledgements, funding statements, and author notes are omitted.
    A nucleotide run of 40 or more is replaced before the block is kept.
    """
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as error:
        raise InputError(f"Article XML could not be parsed: {error}") from error
    found: list[tuple[str, str, str]] = []
    _walk(root, (), found)
    return _blocks(found)


def section_text(text: str) -> tuple[Block, ...]:
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
    return _blocks(found)


def write_sections(
    path: Path,
    pmid: str,
    blocks: tuple[Block, ...],
    sha256: str = "",
) -> None:
    """Write the blocks for ``pmid``. The file is the article text.

    ``sha256`` names the document version that produced the blocks.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, object] = {
        "pmid": pmid,
        "blocks": [
            {"id": block.id, "section": block.section, "kind": block.kind, "text": block.text}
            for block in blocks
        ],
    }
    if sha256:
        payload["sha256"] = sha256
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _blocks(found: list[tuple[str, str, str]]) -> tuple[Block, ...]:
    blocks: list[Block] = []
    for section, kind, text in found:
        cleaned = redact(" ".join(text.split()))
        if not cleaned:
            continue
        blocks.append(Block(f"b{len(blocks) + 1:04d}", section or "article", kind, cleaned))
    return tuple(blocks)


def _walk(element: ET.Element, section: tuple[str, ...], found: list[tuple[str, str, str]]) -> None:
    if _skipped(element):
        return
    name = local_name(element.tag)
    if name == "sec":
        title = _direct_title(element)
        section = (*section, title) if title else section
    elif name == "abstract":
        title = _direct_title(element) or "Abstract"
        section = (title,)
    elif name == "p":
        found.append((_label(section), "paragraph", _text(element)))
        return
    elif name == "fig":
        _caption(element, section, "caption", found)
        return
    elif name == "table-wrap":
        _table(element, section, found)
        return
    elif name == "supplementary-material":
        _caption(element, section, "caption", found)
        return
    for child in list(element):
        _walk(child, section, found)


def _skipped(element: ET.Element) -> bool:
    name = local_name(element.tag)
    if name in _SKIP:
        return True
    if name != "sec":
        return False
    sec_type = element.attrib.get("sec-type", "").casefold()
    title = _direct_title(element).casefold()
    return "ack" in sec_type or "fund" in sec_type or title in _ACK


def _table(
    element: ET.Element, section: tuple[str, ...], found: list[tuple[str, str, str]]
) -> None:
    label = _label(section)
    _caption(element, section, "caption", found)
    header_rows = [row for node in _named(element, "thead") for row in _named(node, "tr")]
    header = _cells(header_rows[-1]) if header_rows else []
    body = [row for row in _named(element, "tr") if row not in header_rows]
    for row in body:
        cells = _cells(row)
        if not any(cells):
            continue
        if header and len(header) == len(cells):
            paired = [
                f"{name}: {cell}" if name else cell
                for name, cell in zip(header, cells, strict=True)
            ]
        else:
            paired = cells
        found.append((label, "table-row", " | ".join(paired)))


def _named(element: ET.Element, name: str) -> list[ET.Element]:
    return [node for node in element.iter() if local_name(node.tag) == name]


def _cells(row: ET.Element) -> list[str]:
    return [_text(cell) for cell in list(row) if local_name(cell.tag) in {"td", "th"}]


def _caption(
    element: ET.Element,
    section: tuple[str, ...],
    kind: str,
    found: list[tuple[str, str, str]],
) -> None:
    for child in list(element):
        if local_name(child.tag) == "caption":
            found.append((_label(section), kind, _text(child)))
            return


def _direct_title(element: ET.Element) -> str:
    for child in list(element):
        if local_name(child.tag) == "title":
            return _text(child)
    return ""


def _label(section: tuple[str, ...]) -> str:
    return " / ".join(part for part in section if part)


def _text(element: ET.Element) -> str:
    return " ".join("".join(element.itertext()).split())
