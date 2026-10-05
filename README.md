# pyVSV

Build a recombinant VSV vaccine plan from a pathogen genome and the Indiana chassis.

Mark the glycoprotein, place it by a method the literature states, and read which published
vaccines used that same arrangement.

<p align="center">
  <img src="https://www.arpaia.io/images/vsv/vsv-structure-genome.jpg" alt="VSV Structure and Genome">
</p>
<p align="center">
  <em>Liu G, Cao W, Salawudeen A, et al. Vesicular Stomatitis Virus: From Agricultural Pathogen to Vaccine Vector. Pathogens. 2021;10(9):1092. Published 2021 Aug 27. doi:<a href="https://doi.org/10.3390/pathogens10091092">10.3390/pathogens10091092</a></em>
</p>

## Install

```sh
pip install pyVSV
```

The distribution is `pyVSV`; the import package is `vsv`. Python 3.12 or newer. Compiling a plan to
an OT-2 bundle needs the `opentrons` extra of [`lab-compiler`](https://pypi.org/project/lab-compiler/).

## Quick start

```python
import vsv
from lab import compile
from lab.experiments.cloning import Assembly, BSAI, assembly_deck, build_assembly
from lab.targets import LiquidHandler
from vsv import AntigenClass, Gene, Replication

virus = vsv.read_record("outbreak.gb")  # FASTA, GenBank, or SBOL3
glycoprotein = virus.mark(virus.cds("GP"), AntigenClass.filovirus_glycoprotein)

chassis = vsv.indiana()
vaccine = chassis.replace(
    chassis.feature(Gene.GLYCOPROTEIN),
    glycoprotein,
    replication=Replication.COMPETENT,
    identifier="rvsv_dg_outbreak_gp",
    name="rVSV outbreak glycoprotein",
)
plasmid = vaccine.plasmid()

compile(
    build_assembly(
        [
            Assembly(
                product=plasmid.product,
                backbone=plasmid.backbone.component,
                parts=[plasmid.insert],
                restriction_enzyme=BSAI,
            )
        ],
        name=vaccine.id,
    ),
    deck=assembly_deck(),
    liquid_handler=LiquidHandler.OT2,
)
```

## Marking the antigen

`mark` records the antigen class and the coding sequence.

```python
glycoprotein = virus.mark(virus.cds("GP"), AntigenClass.filovirus_glycoprotein)
```

For a class the library has not named yet, construct one from its slug:

```python
AntigenClass("paramyxovirus-glycoprotein")
```

A FASTA record carries no features, so declare the span first. Spans are zero-based and half-open.

```python
glycoprotein = virus.mark(virus.add_cds("GP", start, end), AntigenClass.filovirus_glycoprotein)
```

## Placing the gene

`replace` removes VSV G and puts the marked gene in its place. The new document cites the 2004
glycoprotein-replacement method.

```python
vaccine = chassis.replace(
    chassis.feature(Gene.GLYCOPROTEIN),
    glycoprotein,
    replication=Replication.COMPETENT,
    identifier="rvsv_dg_outbreak_gp",
    name="rVSV outbreak glycoprotein",
)
```

`replication` is `Replication.COMPETENT` when the payload supplies entry, and
`Replication.SINGLE_CYCLE` when it does not.

To keep VSV G, or to place the gene first in the genome, use `insert`. It keeps the plan's existing
replication class, so it takes no `replication` argument. The published insertion sites are the
leader and the glycoprotein; anything else is rejected.

```python
# Keep VSV G, add the gene after it. Cites the 1997 method.
vaccine = chassis.insert(
    glycoprotein,
    after=chassis.feature(Gene.GLYCOPROTEIN),
    identifier="rvsv_gl_outbreak_gp",
    name="rVSV outbreak glycoprotein",
)

# Place the gene first, after the leader. Cites the 2014 method.
vaccine = chassis.insert(
    glycoprotein,
    after=chassis.leader,
    identifier="rvsv_leader_outbreak_gp",
    name="rVSV outbreak glycoprotein",
)
```

`vsv.indiana()` is the curated chassis. To start from a sequenced Indiana genome instead:

```python
chassis = vsv.plan(vsv.load("data/raw/NC_001560.fasta"))
```

## Checking a plan

`check()` is true when the plan passes, and failures are on `violations`:

```python
result = vaccine.check()
if not result:
    for violation in result.violations:
        print(violation)
```

`vaccine.observations` is empty; a new plan has no measured phenotypes. `vaccine.write(path)`
stores the plan without the coding sequence.

## Finding precedents

`precedents` compares a plan against a library of published ones and returns those that arrange a
gene the same way. `vsv.plans()` reads the curated library.

```python
for precedent in vaccine.precedents(vsv.plans()):
    print(precedent.id)
```

Filters narrow the result. `antigen_class` keeps plans whose antigen is that class or a subclass of
it; `phenotype` and `host` keep plans with a matching observation.

```python
from vsv import Host, Phenotype

vaccine.precedents(
    vsv.plans(),
    antigen_class=AntigenClass.filovirus_glycoprotein,
    phenotype=Phenotype.PROTECTION,
    host=Host.NONHUMAN_PRIMATE,
)
```

A plan with an unknown transcribed locus, or with any number of non-viral genes other than one, is
not compared and returns no precedents.

`data/designs/` holds one directory per published architecture, each with `design.ttl` and
`PROVENANCE.md`. The stored plan has no insert sequence.

The curated library lives in the repository, not in the wheel. `vsv.plans()` finds it when you work
from a checkout; from an installed package, pass the directory:

```python
from pathlib import Path

library = vsv.plans(Path("path/to/pyVSV/data/designs"))
```

## Plasmid and assembly

`plasmid.sequence` is the transcription unit: the Indiana G gene-start, the marked coding sequence,
and the gene-end. The backbone is the Indiana antigenome NC_001560.1, and the slot is the span that
unit replaces.

`plasmid.check()` fails on:

- an unclosed reading frame
- a gene-end inside the unit
- a replication-competent glycoprotein replacement with no hydrophobic anchor in its last 60 residues

`plasmid.write(path)` stores the unit and that antigenome.

`build_assembly` records the product, the backbone component, and the insert, transferring 2 µL.
`compile` writes the OT-2 bundle to `~/.lab/<plan-id>/OT-2/`. `assembly_deck()` places the reagent
block on the temperature module and the product plate on the thermocycler.

## Command line

```bash
vsv genes data/raw/NC_001560.fasta
vsv extract -g g --kind protein data/raw/NC_001560.gb
vsv designs
vsv show rvsv_dg_filovirus_gp
vsv check data/designs/vsv_indiana
vsv precedents rvsv_dg_outbreak_gp.ttl --antigen-class filovirus-glycoprotein --phenotype protection
vsv plasmid rvsv_dg_outbreak_gp_plasmid.ttl
```

On the command line the genes are `n`, `p`, `m`, `g`, and `l`. The full command reference (every
subcommand, its filters, and what an empty `precedents` answer means) is in
[docs/cli.md](https://github.com/marpaia/pyVSV/blob/master/docs/cli.md).

Corpus commands manage the literature workspace:

```bash
vsv ingest --all                            # store the PubMed corpus under data/corpus/
vsv papers --plan rvsv_dg_filovirus_gp      # papers mentioning an alias of that design
vsv promote <cluster-id>                    # write one accepted cluster into data/designs/
```

## Literature

### Methods the library cites

`replace` and `insert` record the citation for the method they implement, so a written plan carries
the paper its arrangement came from.

| Arrangement | Call | Method paper |
|---|---|---|
| Glycoprotein replacement | `replace(...)` | Garbutt et al., *Journal of Virology*, May 2004, [10.1128/jvi.78.10.5458-5465.2004](https://doi.org/10.1128/jvi.78.10.5458-5465.2004) |
| Gene after G | `insert(..., after=chassis.feature(Gene.GLYCOPROTEIN))` | [10.1128/jvi.71.8.5982-5989.1997](https://doi.org/10.1128/jvi.71.8.5982-5989.1997) |
| Gene after the leader | `insert(..., after=chassis.leader)` | [10.1128/jvi.00019-14](https://doi.org/10.1128/jvi.00019-14) |

Each curated design carries its own citations and quoted evidence in its `PROVENANCE.md`. The
literature workspace (collections, documents, and extractions) is described in
[docs/literature.md](https://github.com/marpaia/pyVSV/blob/master/docs/literature.md).

### Seed reviews

These seven reviews define the scope of the curated evidence. `vsv literature sync review-references`
walks their reference lists, and every paper they cite becomes a member of the `review-references`
collection, so the provenance of any precedent in this library traces back to one of them. The
snapshot, with retrieval date and per-paper citing review, is `data/census/index.json`.

| Review | PMID |
|---|---|
| Lichty et al. 2004, [Vesicular stomatitis virus: re-inventing the bullet](https://pubmed.ncbi.nlm.nih.gov/15121047/) | 15121047 |
| Clarke et al. 2006 | [16977404](https://pubmed.ncbi.nlm.nih.gov/16977404/) |
| Geisbert & Feldmann 2011 | [21987744](https://pubmed.ncbi.nlm.nih.gov/21987744/) |
| Monath et al. 2016 | [27395563](https://pubmed.ncbi.nlm.nih.gov/27395563/) |
| Suder et al. 2018 | [29757706](https://pubmed.ncbi.nlm.nih.gov/29757706/) |
| Fathi et al. 2019 | [31368826](https://pubmed.ncbi.nlm.nih.gov/31368826/) |
| Liu et al. 2021, [Vesicular Stomatitis Virus: From Agricultural Pathogen to Vaccine Vector](https://pubmed.ncbi.nlm.nih.gov/34578125/) | 34578125 |

The frozen PubMed query behind the other collection is `"vesicular stomatitis virus"[tiab]`, or the
Indiana or New Jersey MeSH heading; `data/corpus/ledger.json` records its counts and SHA-256 digest.

### Vector design and recovery

- [Vesicular stomatitis virus (VSV) as a vaccine vector for immunization against viral infections](https://pubmed.ncbi.nlm.nih.gov/24379275/) - Jan 2013
- [An efficient plasmid-based system for the recovery of recombinant vesicular stomatitis virus encoding foreign glycoproteins](https://www.nature.com/articles/s41598-024-65384-8) - Jun 2024
- [Leveraging Synthetic Virology for the Rapid Engineering of Vesicular Stomatitis Virus (VSV)](https://pubmed.ncbi.nlm.nih.gov/39459973/) - Oct 2024
- [Low-dose VSV-EBOV vaccination provides rapid protection from lethal Ebola virus challenge](https://pubmed.ncbi.nlm.nih.gov/41727071/) - Feb 2026

### Reference genomes

- [Complete Genome Sequences of Two Vesicular Stomatitis Virus Isolates Collected in Mexico](https://pmc.ncbi.nlm.nih.gov/articles/PMC5597772/) - Sep 2017 ([VSV Isolate 1](https://www.ncbi.nlm.nih.gov/nuccore/MF196237), [VSV Isolate 2](https://www.ncbi.nlm.nih.gov/nuccore/MF196238))

## Data

- [VSV Genome Sequence](https://www.ncbi.nlm.nih.gov/nuccore/NC_001560) - Complete genome sequence of Vesicular Stomatitis Virus (VSV)
- [VSV Protein Structures](https://www.rcsb.org/structure/6VSV) - 3D structures of VSV proteins available in the Protein Data Bank (PDB)
- [Vesiculovirus on ViralZone](https://viralzone.expasy.org/21) - Information on the Vesiculovirus genus, including VSV, on ViralZone

## Resources

- [VSV Gene Expression Vector on VectorBuilder](https://en.vectorbuilder.com/resources/vector-system/pVSV_Exp.html)

## License

Apache-2.0. See [LICENSE](https://github.com/marpaia/pyVSV/blob/master/LICENSE).
