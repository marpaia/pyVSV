# Construction-paper census

`index.json` is the bibliography. It is the union of PubMed identifiers cited by seven reviews: Lichty 2004 (PMID 15121047), Clarke 2006 (PMC7079905), Geisbert 2011 (PMC3218670), Monath 2016 (PMC5220644), Suder 2018 (PMC6183239), Fathi 2019 (PMC6816421), and Liu 2021 (PMC8470541). Reference lists come from Europe PMC. A PMCID on a row comes from the PMC id converter, not from the frozen corpus flag.

`manifest.jsonl` is one download attempt per cited PMID. `status` is `xml-body`, `front-matter-only`, `no-pmc`, or `error`. An `xml-body` file is `data/corpus/fulltext/{PMCID}.xml`. A body is usable when its text is at least 2,000 characters. Front matter is not saved.

```sh
VSV_NCBI_EMAIL=name@example.edu vsv census index
VSV_NCBI_EMAIL=name@example.edu vsv census fetch
vsv census section
```

`vsv census section` writes `sections/{pmid}.json` for each `xml-body` row. A block is one paragraph, caption, or table row. Those files are gitignored because they contain the article.

`VSV_NCBI_EMAIL` is the contact address NCBI asks for on E-utilities requests. The commands sleep between requests. Re-running `vsv census fetch` leaves an existing XML file in place when its body still passes the length test.

`acquired/`, `sections/`, and `drafts/` are local caches and are gitignored. They are not part of this bibliography yet.
