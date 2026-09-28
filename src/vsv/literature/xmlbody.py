"""Decide whether a PMC XML response contains an article body."""

from __future__ import annotations

import xml.etree.ElementTree as ET

BODY_MIN = 2000


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


def _is_article(root: ET.Element) -> bool:
    if local_name(root.tag) in {"article", "pmc-articleset"}:
        return True
    return any(local_name(element.tag) == "article" for element in root.iter())
