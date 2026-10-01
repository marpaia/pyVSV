# Command line

`vsv` with no subcommand prints help and exits 0. `vsv --version` prints the version. A usage error exits 2. A rejected input prints `vsv: error: …` and exits 1. A gene-naming warning prints `vsv: warning: …` and exits 0.

Paper records are in [literature.md](literature.md). Precedent scope is in `data/census/README.md`. On the command line the genes are `n`, `p`, `m`, `g`, and `l`.

## Read a genome

One or more FASTA, GenBank, or SBOL3 records.

```sh
# Length, base counts, ambiguous bases, and GC content.
vsv summary data/raw/NC_001560.fasta

# The n, p, m, g, and l genes and coding sequences. A naming problem is a warning.
vsv genes data/raw/NC_001560.fasta

# Leader, genes, signals, coding sequences, and trailer. Existing annotations are kept.
vsv annotate data/raw/NC_001560.gb

# GFF3. Other formats: turtle, ntriples, rdfxml, jsonld. The default format is a table.
vsv annotate --format gff3 data/raw/NC_001560.gb

# FASTA for the named genes. Repeat -g; the default is all five. --kind is gene, cds, or protein (default cds).
vsv extract -g g --kind protein data/raw/NC_001560.gb

# k-mer similarity and per-gene protein identity for at least two records. -k defaults to 12.
vsv compare -k 12 data/raw/NC_001560.fasta other.fasta
```

## Inspect the design library

Each curated plan is `data/designs/<id>/` with `design.ttl` and `PROVENANCE.md`. The plan stores no insert sequence.

```sh
# Curated plans, plus extracted clusters when extractions.jsonl exists. No match prints "no designs match".
vsv designs

# Exact antigen-class slug stored on the plan. Subclasses stay under their own slug.
vsv designs --antigen-class bandavirus-glycoprotein

# Plans with a matching observation. Phenotype: expression, growth, immunogenicity, safety, protection. Host: cell-line, mouse, hamster, nonhuman-primate, human, other.
vsv designs --phenotype protection --host nonhuman-primate

# Minimum evidence: in-vitro, rodent, nonhuman-primate, or human. preclinical means rodent.
vsv designs --evidence nonhuman-primate

# One design by id: backbone, modules, implementations, observations, unmodeled fields, and gaps. Route and challenge dose are omitted.
vsv show rvsv_dg_filovirus_gp

# The same view from a design directory. An unknown id is printed as an extraction cluster when one exists.
vsv show data/designs/rvsv_dg_filovirus_gp

# Structural violations exit 1. A gap alone (locus-unknown, parent-unknown, replication-unknown) exits 0.
vsv check data/designs/vsv_indiana
```

## Ask which published plans match an arrangement

```sh
# Same arrangement: G replaced, one gene after G, or one gene after the leader. Omits the query plan. --antigen-class includes subclasses.
vsv precedents data/designs/rvsv_dg_filovirus_gp \
  --antigen-class filovirus-glycoprotein \
  --phenotype protection \
  --host nonhuman-primate

# Plans with one gene after the leader. Columns: id, phenotype, host, assay, citation.
vsv precedents data/designs/rvsv_luc

# Protection rows for plans with one gene after G. A plan with no observation has no row.
vsv precedents data/designs/rvsv_gfp --phenotype protection

# "no precedents": outside those arrangements, an unknown locus, or a filter with no stored observation.
vsv precedents data/designs/rvsv_nj_lcrv
```

## Read a plasmid file

```sh
# Plan, edit, slot, insert, backbone, and product stored in a plasmid SBOL file.
vsv plasmid rvsv_dg_outbreak_gp_plasmid.ttl
```

## Papers

```sh
# Europe PMC review queue under data/ingest. Prints papers, proposals, and cursor.
vsv ingest --since 2024-01-01

# Frozen PubMed corpus under data/corpus. --since or a non-default --query is an error with --all.
vsv ingest --all

# Coverage ledger. Exits 1 when the digest or the stored PMIDs disagree. --dest selects the corpus directory.
vsv corpus

# Classify every stored paper into extractions.jsonl. Prints papers, vectors, unclustered, and fulltext.
vsv corpus extract

# Also read open-access methods for vector papers, from the local cache or one fetch.
vsv corpus extract --fulltext

# Papers whose title or abstract mentions an alias of this design. --limit defaults to 40, and 0 lists every match.
vsv papers --plan rvsv_dg_filovirus_gp

# Exact class name from literature.md, and vector papers that stayed unclustered. No match prints "no papers match".
vsv papers --class vaccine-vector --unclustered --limit 0

# Copy one extracted cluster to data/designs/<id>/ as design.ttl and PROVENANCE.md. Needs a new directory, a known edit, and an antigen class.
vsv promote <cluster-id>

# Review bibliography at data/census/index.json. Needs VSV_NCBI_EMAIL.
vsv census index

# Article XML into data/corpus/fulltext. --limit 0 fetches every cited paper.
vsv census fetch --limit 0

# Quotable blocks at data/census/sections/{pmid}.json for each xml-body row in the manifest.
vsv census section

# Paper count, collection overlap, citations, documents, readings, and work still open.
vsv literature status

# Each snapshot's name, size, date, rule, and provenance. An empty workspace prints "no collections".
vsv literature collections

# Refresh the frozen PubMed snapshot into the pubmed collection.
vsv literature sync pubmed

# Refresh the review bibliography. Needs VSV_NCBI_EMAIL.
vsv literature sync review-references

# papers filters limited to one collection. Collections are pubmed and review-references.
vsv literature papers --collection pubmed --plan rvsv_dg_filovirus_gp

# One PMID with its collections, the papers it cites, the papers that cite it, and its documents.
vsv literature show 15121047

# Documents for that collection. --limit 0 fetches every member. review-references needs VSV_NCBI_EMAIL.
vsv literature fetch --collection review-references --limit 0

# Quotable blocks for available local documents in that collection. Prints readings.
vsv literature parse --collection review-references

# Exits 1 when the workspace disagrees with itself. Prints ok, or one message per disagreement.
vsv literature check
```
