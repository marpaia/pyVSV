"""Article XML: body checks, license text, section blocks, and methods text."""

from __future__ import annotations

import xml.etree.ElementTree as ET

from vsv.errors import InputError
from vsv.literature.documents.text import clean_text, text_blocks
from vsv.literature.records import TextBlock

BODY_MIN = 2000

_SKIP = frozenset({"ref-list", "ack", "funding-group", "author-notes", "fn-group"})
_ACK = frozenset({"acknowledgements", "acknowledgments", "funding"})


def local_name(tag: str) -> str:
    """The element name without an XML namespace."""
    return tag.rsplit("}", 1)[-1]


def body_characters(xml: str) -> int | None:
    """Characters of article body text, or ``None`` when ``xml`` is not an article.

    Zero means the document parsed and has no body text. A count of ``BODY_MIN`` or
    more is a usable article.
    """
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return None
    if not _is_article(root):
        return None
    found = False
    total = 0
    for element in root.iter():
        if local_name(element.tag) != "body":
            continue
        found = True
        total += len(" ".join("".join(element.itertext()).split()))
    if not found:
        return 0
    return total


def license_sentence(xml: str) -> str:
    """The first license statement in ``xml``, or an empty string."""
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return ""
    for element in root.iter():
        if local_name(element.tag) != "license":
            continue
        text = " ".join("".join(element.itertext()).split())
        if text:
            return text[:240]
    return ""


def section_xml(xml: str) -> tuple[TextBlock, ...]:
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
    return text_blocks(found)


def methods_text(xml: str) -> str:
    """The methods section of a JATS full-text document, with nucleotide runs removed.

    The longest matching section is kept. Nested method sections are shorter than their parent.
    """
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return ""
    texts = ["".join(element.itertext()) for element in root.iter() if _is_methods(element)]
    if not texts:
        return ""
    return clean_text(max(texts, key=len))


def _is_article(root: ET.Element) -> bool:
    if local_name(root.tag) in {"article", "pmc-articleset"}:
        return True
    return any(local_name(element.tag) == "article" for element in root.iter())


def _is_methods(element: ET.Element) -> bool:
    if not str(element.tag).endswith("sec"):
        return False
    title = ""
    for child in list(element):
        if str(child.tag).endswith("title"):
            title = "".join(child.itertext())
            break
    label = f"{element.attrib.get('sec-type', '')} {title}".casefold()
    return "method" in label or "material" in label


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
