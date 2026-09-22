"""The frozen PubMed query that defines the VSV corpus."""

from __future__ import annotations

from pathlib import Path

CORPUS_ROOT = Path(__file__).resolve().parents[3] / "data" / "corpus"
PAPERS_NAME = "papers.jsonl"
EXTRACTIONS_NAME = "extractions.jsonl"
LEDGER_NAME = "ledger.json"
AUDIT_NAME = "audit.json"
FULLTEXT_NAME = "fulltext"

# Title/abstract phrase plus the Indiana and New Jersey MeSH headings.
# A later ingest diffs this term against PubMed and records the count in the ledger.
PUBMED_TERM = (
    '"vesicular stomatitis virus"[tiab] OR '
    '"Vesicular stomatitis Indiana virus"[mh] OR '
    '"Vesicular stomatitis New Jersey virus"[mh]'
)

# Papers that say only "VSV" and are not already in PUBMED_TERM.
# The count is recorded. The hits are not part of the corpus.
RESIDUAL_TERM = f'"VSV"[tiab] NOT ({PUBMED_TERM})'

AUDIT_SEED = 20260921
AUDIT_SIZE = 100

PUBMED_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EUROPEPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUROPEPMC_FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest"

VECTOR_KINDS = frozenset({"vaccine-vector", "oncolytic", "reporter"})


def mentions(text: str, aliases: tuple[str, ...] | list[str]) -> bool:
    """Whether ``text`` contains one of ``aliases``, ignoring case."""
    folded = text.casefold()
    return any(alias.casefold() in folded for alias in aliases if alias)
