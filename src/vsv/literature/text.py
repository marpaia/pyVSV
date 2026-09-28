"""Parse local documents into readings."""

from __future__ import annotations

from pathlib import Path

from vsv.errors import InputError
from vsv.literature.records import PUBMED, REVIEW_REFERENCES, Block, Document
from vsv.literature.section import section_text, section_xml, write_sections
from vsv.literature.store import acquired_dir, load_catalog, sections_dir

_FORMAT_RANK = {"xml": 0, "html": 1, "pdf": 2}


def parse_collection(root: Path, collection: str) -> int:
    """Write a reading for each available document in ``collection``.

    XML is parsed from the stored article. HTML and PDF readings use the plain
    text saved beside the file. One reading is kept per paper, preferring XML.
    """
    if collection not in {PUBMED, REVIEW_REFERENCES}:
        raise InputError(f"{collection}: unknown collection.")
    catalog = load_catalog(root)
    selected = catalog.named(collection)
    if selected is None:
        raise InputError(f"{collection}: collection is not in this workspace.")
    members = set(selected.pmids)
    chosen: dict[str, Document] = {}
    for document in catalog.documents:
        if document.pmid not in members or document.availability != "available":
            continue
        current = chosen.get(document.pmid)
        if current is None or _rank(document) < _rank(current):
            chosen[document.pmid] = document
    written = 0
    directory = sections_dir(root)
    for pmid in sorted(chosen, key=int):
        document = chosen[pmid]
        blocks = _blocks(root, document)
        if not blocks:
            continue
        write_sections(directory / f"{pmid}.json", pmid, blocks, document.sha256)
        written += 1
    return written


def _blocks(root: Path, document: Document) -> tuple[Block, ...]:
    if document.format == "xml" and document.location:
        path = root / document.location
        if not path.is_file():
            return ()
        return section_xml(path.read_text(encoding="utf-8"))
    if document.format in {"html", "pdf"}:
        path = acquired_dir(root) / f"{document.pmid}.txt"
        if not path.is_file():
            return ()
        return section_text(path.read_text(encoding="utf-8"))
    return ()


def _rank(document: Document) -> int:
    return _FORMAT_RANK.get(document.format, 9)
