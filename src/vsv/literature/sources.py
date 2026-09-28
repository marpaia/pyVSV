"""NCBI and Europe PMC calls that discover papers and fetch article files."""

from vsv.literature.bibliography import (
    FetchJson,
    FetchText,
    efetch_pmc_url,
    europepmc_fulltext_url,
    idconv_url,
    oai_url,
    references_url,
)
from vsv.literature.bibliography import (
    pubmed_url as pubmed_record_url,
)
from vsv.literature.fetch import (
    Fetch,
    clean_text,
    europepmc_records,
    fulltext_url,
    live_fetch,
    live_fetch_text,
    methods_text,
    pubmed_count,
    pubmed_ids,
)

__all__ = [
    "Fetch",
    "FetchJson",
    "FetchText",
    "clean_text",
    "efetch_pmc_url",
    "europepmc_fulltext_url",
    "europepmc_records",
    "fulltext_url",
    "idconv_url",
    "live_fetch",
    "live_fetch_text",
    "methods_text",
    "oai_url",
    "pubmed_count",
    "pubmed_ids",
    "pubmed_record_url",
    "references_url",
]
