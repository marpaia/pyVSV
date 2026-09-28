"""Parsed document text and the checksum of the file it came from."""

from __future__ import annotations

import json
from pathlib import Path

from vsv.errors import InputError
from vsv.literature.records import DocumentText, TextBlock
from vsv.literature.storage.paths import sections_dir


def write_sections(
    path: Path,
    pmid: str,
    blocks: tuple[TextBlock, ...],
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


def load_texts(root: Path) -> tuple[DocumentText, ...]:
    """Parsed text stored under the workspace sections directory."""
    directory = sections_dir(root)
    if not directory.is_dir():
        return ()
    readings: list[DocumentText] = []
    for path in sorted(directory.glob("*.json")):
        readings.append(_reading(path))
    readings.sort(key=lambda item: int(item.pmid) if item.pmid.isdigit() else 0)
    return tuple(readings)


def _reading(path: Path) -> DocumentText:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise InputError(f"{path}: {error.msg}.") from error
    if not isinstance(raw, dict):
        raise InputError(f"{path}: reading must be an object.")
    pmid = raw.get("pmid")
    blocks = raw.get("blocks")
    if not isinstance(pmid, str) or not isinstance(blocks, list):
        raise InputError(f"{path}: reading needs pmid and blocks.")
    version = raw.get("sha256")
    parsed: list[TextBlock] = []
    for item in blocks:
        if not isinstance(item, dict):
            raise InputError(f"{path}: each block must be an object.")
        parsed.append(
            TextBlock(
                _string(item, "id", path),
                _string(item, "section", path),
                _string(item, "kind", path),
                _string(item, "text", path),
            )
        )
    digest = version if isinstance(version, str) else ""
    return DocumentText(pmid, digest, tuple(parsed))


def _string(raw: dict[str, object], key: str, path: Path) -> str:
    value = raw.get(key)
    if isinstance(value, str):
        return value
    raise InputError(f"{path}: {key} must be a string.")
