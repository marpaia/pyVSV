# Data conventions

Keep downloaded sequence files in `data/raw/` and derived intermediate files in `data/processed/`. Treat raw inputs as immutable and regenerate processed data through documented commands.

`data/designs/` is tracked. The directories are the Indiana chassis and the filovirus vaccines. Each has `design.ttl`, an SBOL3 design file, and `PROVENANCE.md`, the reading note. The Turtle describes the ordered genome and the papers: antigen class, host, assay, and citation. A parent or replication class may be unknown. A nucleotide sequence is present only when the antigen feature has an `sbol:Sequence`. The shelf cards have none.

For each dataset, add a small, tracked provenance document beside this README. Record its source URL or accession with version, retrieval date, checksum, local filename, and any usage or attribution requirements. Record the commands and parameters used for derived data, together with the input identifiers and code revision.

Small synthetic test fixtures belong in `tests/fixtures/` when needed. Large datasets and generated outputs should use external storage appropriate to the project; this scaffold does not configure storage or download data.
