# VSV literature corpus

`papers.jsonl` is one PubMed record per line: PMID, DOI, PMCID, year, publication types, MeSH headings, title, abstract, open-access flag, source URL, and retrieval date. `extractions.jsonl` is the class and the architecture fields for each paper. `ledger.json` is the query, the PubMed count, the PMID list, the Europe PMC misses, and the residual count.

The query is:

```
"vesicular stomatitis virus"[tiab]
OR "Vesicular stomatitis Indiana virus"[mh]
OR "Vesicular stomatitis New Jersey virus"[mh]
```

`vsv ingest --all` refreshes the files from PubMed and Europe PMC. `vsv corpus` checks that every ledger PMID is stored. `vsv corpus extract` rewrites the extractions from the local papers. A value is `unknown` when the text does not state it. The sentence that supports a known value is stored with it.

`fulltext/` is the open-access XML cache used by `vsv corpus extract --fulltext`. It is gitignored. Titles, abstracts, and quotes have runs of 40 or more nucleotides removed before they are saved. The corpus does not hold insert sequences.

`audit.json`, when present, is a hand-labeled sample. `vsv corpus` prints the classifier's misses on that sample.

A cluster is papers that share a backbone, replication class, edit, antigen, and stated locus. `vsv promote` copies one cluster into `data/designs/` as `design.ttl` and does not run on its own.
