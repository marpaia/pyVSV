"""Local queries over the catalog."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from vsv.literature.cluster import clusters_from
from vsv.literature.query import VECTOR_KINDS, mentions
from vsv.literature.records import Extraction, Paper


def matching_papers(
    papers: Sequence[Paper],
    extractions: Sequence[Extraction] | None = None,
    *,
    members: set[str] | None = None,
    year: str | None = None,
    mesh: str | None = None,
    text: str | None = None,
    kind: str | None = None,
    antigen: str | None = None,
    unresolved: bool = False,
    unclustered: bool = False,
    aliases: Sequence[str] = (),
    use_aliases: bool = False,
    show_kind: bool = False,
) -> list[tuple[Paper, str]]:
    """Papers that match the filters, newest first.

    ``members`` restricts the result to one collection. ``None`` searches every paper
    that was passed in. The label is the extraction class when ``show_kind`` is set.
    """
    by_pmid: Mapping[str, Extraction] = {}
    if extractions is not None:
        by_pmid = {item.pmid: item for item in extractions}
    placed: set[str] = set()
    if unclustered and extractions is not None:
        placed = {pmid for cluster in clusters_from(list(extractions)) for pmid in cluster.pmids}
    matched: list[tuple[Paper, str]] = []
    for paper in papers:
        if members is not None and paper.pmid not in members:
            continue
        item = by_pmid.get(paper.pmid)
        if kind is not None and (item is None or item.kind != kind):
            continue
        if unresolved and (item is None or item.status != "unresolved"):
            continue
        if unclustered and (item is None or item.kind not in VECTOR_KINDS or paper.pmid in placed):
            continue
        if year is not None and paper.year != year:
            continue
        if mesh is not None and not any(mesh.casefold() in term.casefold() for term in paper.mesh):
            continue
        if (
            text is not None
            and text.casefold() not in f"{paper.title}\n{paper.abstract}".casefold()
        ):
            continue
        if use_aliases and not mentions(f"{paper.title}\n{paper.abstract}", tuple(aliases)):
            continue
        if antigen is not None:
            if item is None:
                continue
            blob = f"{item.antigen_name.value}\n{item.antigen_slug.value}"
            if antigen.casefold() not in blob.casefold():
                continue
        label = item.kind if show_kind and item is not None else "-"
        matched.append((paper, label))
    matched.sort(key=lambda item: (item[0].year or "-", int(item[0].pmid)), reverse=True)
    return matched
