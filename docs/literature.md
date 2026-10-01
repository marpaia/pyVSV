# Literature

Command order for the whole program is in [cli.md](cli.md).

`vsv.literature` keeps a workspace of VSV papers. For a PMID it stores any combination of a bibliographic identity, a collection membership, a saved article or the reason none was saved, text parsed from that article, and a quote-backed extraction. `load_catalog()` reads the workspace into one `Catalog`.

A curated genome plan is a design. `vsv promote` writes one from an extraction cluster, and does not change the literature files.

The workspace defaults to `data/` in the repository. `vsv literature --root` selects another. `corpus/` is the PubMed collection. `census/` is the bibliography assembled from review reference lists, plus the article files saved for those papers.

```text
data/
  corpus/
    papers.jsonl          one paper per line, in PMID order
    extractions.jsonl     one extraction per line, in PMID order
    ledger.json           PubMed snapshot: query, counts, PMID list, digest
    audit.json            hand-labeled sample, when one has been written
    fulltext/             {PMCID}.xml, gitignored
  census/
    index.json            reference snapshot: retrieval date, reviews, cited papers
    manifest.jsonl        one download attempt per line, in PMID order
    acquired/             {PMID}.html, .pdf, and .txt, gitignored
    sections/             {PMID}.json parsed blocks, gitignored
    drafts/               local notes, gitignored
  designs/                curated plans written by vsv promote
  ingest/                 Europe PMC review queue written by vsv ingest
```

`data/ingest/` is the Europe PMC review queue, separate from `corpus/` and `census/`. `papers/`, `queue/`, and `cursor.json` under it are gitignored. Tracked literature files are the bibliographic record. Article bodies are a local cache.

## One PMID, in order

Each step writes its own record and leaves the earlier ones in place.

```text
identity  ->  membership  ->  document  ->  parsed text  ->  extraction
```

### Identity

`Paper` is the bibliographic identity of one PMID: title, abstract, year, DOI, PMCID, publication types, MeSH headings, open-access flag, source URL, and retrieval date. The catalog holds at most one paper for a PMID.

`vsv literature sync pubmed` and `vsv ingest --all` ask PubMed for every PMID matching the frozen query, then ask Europe PMC for the bibliographic record. The query is `"vesicular stomatitis virus"[tiab]`, or the Indiana or New Jersey MeSH heading. When Europe PMC returns a PMID, the paper is stored with `in_europepmc` set. A PMID that PubMed lists and Europe PMC does not return is a stub, with an empty title and abstract. PMIDs that leave the PubMed result are dropped. An existing Europe PMC record is kept when its PMID is still in the result.

`corpus/ledger.json` is the snapshot of that result: `query`, `residual_query`, `retrieved`, `pubmed_count`, `residual_count`, `pmids`, `europepmc_misses`, and `digest`. The digest is the SHA-256 of the PMID list, one identifier per line. `residual_query` counts papers that say only `"VSV"` and fall outside the frozen query. That count is recorded. Those papers are not stored.

`vsv literature sync review-references` and `vsv census index` read the reference lists of seven seed reviews, confirm PMC identifiers, and write one cited-paper row per PMID. Both commands need `VSV_NCBI_EMAIL`, the contact address NCBI asks for on E-utilities requests. Requests are paced. The reviews are Lichty 2004 (15121047), Clarke 2006 (16977404), Geisbert 2011 (21987744), Monath 2016 (27395563), Suder 2018 (29757706), Fathi 2019 (31368826), and Liu 2021 (34578125).

If the PMID is already in `papers.jsonl`, that record's title, year, DOI, and open-access flag stay. Otherwise title, year, and DOI come from PubMed front matter, when NCBI returns them. `census/index.json` holds `retrieved`, the review list, and `papers`. Each cited paper has `pmid`, `pmcid`, `doi`, `title`, `year`, `open_access`, `in_corpus`, and `cited_by`. `in_corpus` records that `papers.jsonl` already held the PMID when the snapshot was written.

The catalog also adds a paper for a PMID that has no bibliographic row yet. It builds that paper from a cited paper, a manifest row, or a saved file whose name is a PMID.

### Membership

`Collection` is a named set of PMIDs admitted by one rule, with the provenance of that rule. Each member is a `Membership`: the PMID and the date of the snapshot. The two snapshots are the two collections.

| Name | Who belongs | Snapshot |
|---|---|---|
| `pubmed` | Every PMID in the frozen PubMed query | `corpus/ledger.json` |
| `review-references` | Every PMID cited by one of the seed reviews | `census/index.json` |

The PubMed sync replaces the `pubmed` snapshot. The review-reference sync replaces the `review-references` snapshot. A paper's other fields are independent of either list, so a PMID may appear in one snapshot, in both, or in neither. `vsv literature status` counts the overlap.

`Citation` records that one PMID cites another. In this workspace the source is a seed review and the target is a paper in `cited_by`. The target joins `review-references` because the snapshot lists it.

### Document

The status string on `census/manifest.jsonl` says whether an article file was kept, and why it was not. The catalog splits that string into a format (`xml`, `html`, `pdf`, or empty), a method (`download`, `cache`, `browser`, or empty), and an availability. A `Document` points at a location relative to the workspace, carries a SHA-256, and keeps each attempt as an `AcquisitionAttempt`.

A manifest row has `pmid`, `pmcid`, `doi`, `title`, `year`, `source`, `retrieved`, `sha256`, `bytes`, `body_characters`, `license`, `status`, `detail`, and `cited_by`. `sha256` is empty when no file was kept. `source` is the URL that supplied the file, `cache` when an existing XML file was kept, or empty when nothing was fetched.

A usable body is article text of at least 2,000 characters. A parsed article with no body text has a count of zero. Markup that is not an article has no count.

`vsv literature fetch --collection review-references` and `vsv census fetch` run this attempt for each cited paper:

1. No PMCID ends the attempt as `no-pmc`. Nothing is downloaded.
2. `corpus/fulltext/{PMCID}.xml` is kept when it already exists and its body is usable. The status is `xml-body` and the method is `cache`.
3. Otherwise Europe PMC full text is tried, then NCBI efetch. A usable response is written to `corpus/fulltext/{PMCID}.xml` with status `xml-body`.
4. An efetch response that parses but is shorter than a usable body ends as `front-matter-only`. The file on disk is left alone.
5. OAI-PMH is tried when efetch does not return article XML. A usable OAI document is saved as `xml-body`.
6. Three failed requests end as `error`. Any other miss ends as `front-matter-only`.

A row is acquired when its status is `xml-body`, `browser-html`, or `browser-pdf`, or when it already has a checksum. A later attempt that stores no file leaves an acquired row in place. Rows the run does not visit stay in the manifest.

`vsv literature fetch --collection pubmed` downloads Europe PMC XML for open-access corpus papers that have a PMCID. It writes `corpus/fulltext/` and leaves the manifest unchanged. An existing usable file is kept.

`scripts/census_browse.py` saves publisher HTML or PDF through a signed-in browser, at `census/acquired/{PMID}.html` or `.pdf`. The plain text beside that file is `{PMID}.txt`. Those attempts use the browser statuses below.

| Status | What happened | Format | Availability |
|---|---|---|---|
| `xml-body` | A usable article body is saved | `xml` | `available` |
| `browser-html` | Publisher HTML is saved | `html` | `available` |
| `browser-pdf` | A publisher PDF is saved | `pdf` | `available` |
| `front-matter-only` | The response parsed and the body is too short to keep | `xml` | `partial` |
| `needs-login` | The publisher page requires a sign-in | empty | `needs-login` |
| `unavailable` | No usable file was kept | `html` when a checksum is stored, otherwise empty | `unavailable` |
| `no-pmc` | The cited paper has no PMC identifier | empty | `unavailable` |
| `error` | The attempt failed | empty | `error` |

An XML file in `corpus/fulltext/` with no manifest row is still a document. Its location is `corpus/fulltext/{name}.xml`, its method is `download`, and its availability is `available`. The PMID comes from the corpus paper with that PMCID when one exists.

### Parsed text

An `available` document can be parsed into quotable blocks. `vsv literature parse --collection` writes one text file per available document in that collection. XML is preferred to HTML, and HTML to PDF. XML is read from the article file. HTML and PDF are read from `census/acquired/{PMID}.txt`. A paper with no readable file is skipped.

`vsv census section` writes the same kind of file for each manifest row whose status is `xml-body`, from `corpus/fulltext/{PMCID}.xml`.

The file is `census/sections/{PMID}.json`. It has `pmid`, `blocks`, and `sha256` when the parser was given a checksum. The checksum is the SHA-256 of the file that was parsed, so the text stays tied to that version of the document. A block has `id`, `section`, `kind`, and `text`, and is a paragraph, a caption, or a table row. Reference lists, acknowledgements, funding statements, and author notes are omitted from XML. A run of 40 or more nucleotides is removed before a block or a stored quote is kept. The catalog holds these records as `readings`, each a `DocumentText`.

### Extraction

`vsv corpus extract` classifies every paper in `papers.jsonl` and writes `extractions.jsonl`. The class comes from the title, abstract, publication types, and MeSH terms. The classes are `vaccine-vector`, `oncolytic`, `reporter`, `basic-virology`, `structure`, `epidemiology`, `pseudotype-reagent`, `review`, and `other`.

A paper that matches no rule is `other` and `unresolved`. Any other class is `classified`. A publication type of review sets the review flag. The class stays the subject class when one matches. A title that names vesicular stomatitis and matches no other rule is `basic-virology`. A review that matches neither is the class `review`.

A saved article at `census/acquired/{PMID}.txt` is read for papers classed `vaccine-vector`, `oncolytic`, or `reporter` when its text is at least 2,000 characters. With `--fulltext`, open-access methods XML is read from the local cache, or fetched once, for those papers that have no saved article. That text can fill architecture fields. It leaves the class unchanged.

Every corpus paper gets an extraction. Architecture fields are read for those three vector classes. On every other class the fields are `unknown`.

An architecture field is a `QuotedValue`: a value and the sentence it came from. When no sentence states it, the value is `unknown` and the quote is empty. Two sentences that state different values are a conflict. The field stays `unknown`, and the field name is listed in `conflicts`.

Claims are result sentences with a phenotype, host, assay, and outcome, each copied from the text that states them.

Vector papers that state an antigen and an unambiguous edit form a `PaperCluster`. The key is the edit, the antigen, and the stated locus. A paper that states a locus stays apart from one that does not.

An unstated backbone joins a single stated backbone. Two stated backbones stay apart. Unknown backbone, replication, locus, and parent are gaps on the cluster. A paper with no antigen, or with a conflicting edit or antigen, stays unclustered.

`vsv.design.literature` compares each cluster with the curated design library. `match` is the curated architecture id when the edit, antigen, locus, backbone, and replication agree. `vsv promote` writes one accepted cluster to `data/designs/<id>/` as `design.ttl` and `PROVENANCE.md`. The directory must be new. Promotion requires a known edit and an antigen class.

## Status, search, and checks

`vsv literature status` counts papers, the overlap of the two collections, citations, documents by availability, and parsed texts. Pending work is a review-reference member with no document, a document whose availability is `error` or `needs-login`, or an available document with no parsed text.

`vsv literature collections` lists each snapshot's name, size, date, rule, and provenance. `vsv literature show PMID` prints one paper with its collections, citations, and documents.

`vsv literature papers` searches the catalog. `--collection` limits the search to one snapshot. Other filters select by year, MeSH, title or abstract, class, antigen, an unresolved class, a vector paper that stayed unclustered, or an alias of a curated design.

`vsv literature check` exits 1 when the workspace disagrees with itself. A collection member has a paper, and a citation target has a paper. When `ledger.json` is present, its digest matches the PMID list, the list length matches the PubMed count, and `papers.jsonl` contains those PMIDs and no others. An available document names a file inside the workspace, and the file's checksum matches the document. A parsed text names a paper. A PMID listed in only one snapshot is consistent. An unknown architecture field is consistent.

`vsv corpus` prints the ledger and exits 1 when the stored papers disagree with it. When extractions exist, a corpus paper with no extraction is also a failure. The report lists class counts, cluster counts, and unclustered vector papers. When `audit.json` is present, the report prints class mismatches against that hand-labeled sample. Those mismatches leave the exit status unchanged.

## Package

The package root exports `Paper`, `Collection`, `Membership`, `Citation`, `Document`, `DocumentText`, `Catalog`, `load_catalog`, and `default_root`.

| Location | What it does |
|---|---|
| `records.py`, `catalog.py` | The shared records, and the catalog assembled from a workspace |
| `discovery/` | The PubMed query, the seed reviews, and both snapshots |
| `documents/` | Download, source fallback, and XML or plain-text parsing |
| `analysis/` | Classification, quoted extraction, clustering, and the audit |
| `storage/` | Workspace paths, and the reader and writer for each file |
| `search.py` | Filter, sort, and alias matching over papers and extractions |
| `checks.py` | The workspace check and the corpus coverage report |
| `http.py` | JSON and text fetch callables, retries, TLS, and request pacing |

`storage/paths.py` defines the default workspace. Callers pass a directory into the loaders and writers. Tests use a temporary directory and their own fetch callables. Importing the package does not contact the network.
