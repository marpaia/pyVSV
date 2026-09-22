# Literature review queue

`vsv ingest` searches Europe PMC and writes two local directories:

- `papers/` is one JSON record per PubMed identifier: title, abstract, DOI, and retrieval date.
- `queue/` is a pending design proposal when the title or abstract names a tracked antigen class.
- `cursor.json` is the Europe PMC cursor for the next run. `--since` starts a fresh date window and ignores the saved cursor.

A proposal may leave locus, replication, and parent unknown. That is the paper's silence, not a guessed coordinate. Nothing here is copied into `data/designs/` automatically.

The complete PubMed corpus lives in `data/corpus/`. `vsv ingest --all` writes it. This queue is the narrower Europe PMC search.

Runs of 40 or more nucleotide letters in a title or abstract are stored as `[sequence omitted]`. Proposals do not carry insert sequences. A person reviews a pending file before any design is added to the curated library.
