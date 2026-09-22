# Data conventions

Keep downloaded sequence files in `data/raw/` and derived intermediate files in `data/processed/`. Treat raw inputs as immutable and regenerate processed data through documented commands.

`data/designs/` is tracked. Each directory is one VSV architecture, its implementations, qualitative observations, and a provenance note. Those files describe antigen class, host, assay, and citation. A locus, parent, or replication class may be unknown. They do not contain nucleotide or amino-acid sequences.

`data/ingest/` holds the Europe PMC review queue written by `vsv ingest`. `papers/`, `queue/`, and `cursor.json` are local and gitignored. A pending proposal is not a curated design. See `data/ingest/README.md`.

`data/corpus/` is the PubMed VSV corpus: `papers.jsonl`, `extractions.jsonl`, and `ledger.json` are tracked. `fulltext/` is a local cache and is gitignored. See `data/corpus/README.md`.

For each dataset, add a small, tracked provenance document beside this README. Record its source URL or accession with version, retrieval date, checksum, local filename, and any usage or attribution requirements. Record the commands and parameters used for derived data, together with the input identifiers and code revision.

Small synthetic test fixtures belong in `tests/fixtures/` when needed. Large datasets and generated outputs should use external storage appropriate to the project; this scaffold does not configure storage or download data.
