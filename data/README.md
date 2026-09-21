# Data conventions

Keep downloaded source files in `data/raw/` and derived intermediate files in `data/processed/`. Both directories are ignored by Git. Treat raw inputs as immutable and regenerate processed data through documented commands.

For each dataset, add a small, tracked provenance document beside this README. Record its source URL or accession with version, retrieval date, checksum, local filename, and any usage or attribution requirements. Record the commands and parameters used for derived data, together with the input identifiers and code revision.

Small synthetic test fixtures belong in `tests/fixtures/` when needed. Large datasets and generated outputs should use external storage appropriate to the project; this scaffold does not configure storage or download data.
