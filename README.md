# pyVSV

Build a recombinant VSV vaccine plan from a pathogen genome and the Indiana chassis.

Take the glycoprotein coding sequence, place it by a method the literature states, and read
which published vaccines used that same arrangement.

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

The distribution is `pyVSV`; the import package is `vsv`. Python 3.12 or newer.

## Quick start

```python
import vsv
from vsv import Gene, GlycoproteinForm, Host, Replication, Timing
from vsv.genome.sbol import write_genome

genome = vsv.read_record("outbreak.gb")  # FASTA, GenBank, or an SBOL genome document
coding = genome.cds(Gene.GLYCOPROTEIN)
vaccine = vsv.indiana().replace(
    coding,
    identifier="rvsv_dg_outbreak_gp",
    name="rVSV outbreak glycoprotein",
    antigen_class="zaire-ebolavirus-glycoprotein",
    strain="Kikwit-1995",
    form=GlycoproteinForm.FULL_LENGTH,
    replication=Replication.COMPETENT,
)
vaccine.check()
vaccine.precedents(vsv.library(), host=Host.RHESUS, timing=Timing.POST_EXPOSURE)
print(write_genome(vaccine.antigenome()))
vaccine.write("kikwit.ttl")
```

`read_record` reads FASTA, GenBank, or an SBOL genome document. `write_genome` writes that
antigenome as SBOL Turtle. `write` stores the plan as an SBOL3 design file. `vsv.read` opens a
plan and returns a `Vaccine`. Calling code does not handle SBOL objects.

## The coding sequence

`cds` returns the DNA string for one of the five VSV genes and raises when that gene has no
coding sequence. Spans inside a `Genome` are zero-based and half-open. In an SBOL file those
spans are one-based inclusive ranges.

## Placing the gene

`replace` removes the VSV glycoprotein and puts the coding sequence in its place. The new vaccine cites the 2004
glycoprotein-replacement method. The antigen id defaults to `Glycoprotein`.

`replication` is `Replication.COMPETENT` when the antigen supplies entry, and
`Replication.SINGLE_CYCLE` when it does not.

To keep the VSV glycoprotein, or to place the gene first, use `insert`. It keeps the parent's replication class,
so it takes no `replication` argument. The published sites are `"g"` and `"leader"`.

```python
# Keep the VSV glycoprotein and add the gene after it. Cites the 1997 method.
vaccine = vsv.indiana().insert(
    coding,
    after="g",
    identifier="rvsv_gl_outbreak_gp",
    name="rVSV outbreak glycoprotein",
    antigen_class="zaire-ebolavirus-glycoprotein",
    module="outbreak_gp",
    form=GlycoproteinForm.FULL_LENGTH,
    strain="Kikwit-1995",
)

# Place the gene first, after the leader. Cites the 2014 method.
vaccine = vsv.indiana().insert(
    coding,
    after="leader",
    identifier="rvsv_leader_outbreak_gp",
    name="rVSV outbreak glycoprotein",
    antigen_class="zaire-ebolavirus-glycoprotein",
    module="outbreak_gp",
    form=GlycoproteinForm.FULL_LENGTH,
    strain="Kikwit-1995",
)
```

`vsv.indiana()` is nucleoprotein, phosphoprotein, matrix protein, glycoprotein, and large
polymerase, with realization `NC_001560`.

## Checking a vaccine

`check()` is true only when `violations` is empty. A knowledge gap is a violation: unknown
replication, or an unknown parent. `vsv check` exits 0 when every violation is a gap and
exits 1 on any other rule.

When the edited antigen has a coding sequence, `check()` also rejects:

- an unclosed reading frame
- a gene-end inside the transcription unit
- a gene-end joined directly to a gene-start
- a soluble glycoprotein used as a replication-competent glycoprotein replacement
- a full-length or tail-chimera competent glycoprotein replacement with no hydrophobic anchor in its last 60 residues
- a Zaire glycoprotein whose editing site is not eight A's when full-length, or seven A's when soluble

```python
result = vaccine.check()
if not result:
    for violation in result.violations:
        print(violation.rule, violation.message)
```

`vaccine.observations` is empty until a result is recorded. `write` stores the coding sequence
when one is set. The curated shelf cards have none.

## Finding precedents

`precedents` returns published vaccines with the same arrangement: the glycoprotein replaced,
one gene after the glycoprotein, or one gene after the leader. `vsv.library()` reads the curated shelf.

```python
for precedent in vaccine.precedents(vsv.library()):
    print(precedent.id)
```

`antigen_class` keeps vaccines whose antigen is that class or a subclass of it. `phenotype`,
`host`, and `timing` keep vaccines with a matching observation. A filovirus query includes the
species glycoproteins. A Zaire query does not include Sudan or Marburg. A nonhuman-primate query
includes cynomolgus and rhesus.

```python
from vsv import Phenotype

vaccine.precedents(
    vsv.library(),
    antigen_class="filovirus-glycoprotein",
    phenotype=Phenotype.PROTECTION,
    host=Host.NONHUMAN_PRIMATE,
)
```

A vaccine with any number of non-viral genes other than one is not compared.

`data/designs/` holds the Indiana chassis and one directory per published filovirus vaccine.
Each directory has `design.ttl` and `PROVENANCE.md`. The Turtle describes the ordered genome and
the papers. A nucleotide sequence is present only when the antigen feature has an `sbol:Sequence`.
The shelf cards have none.

The curated library lives in the repository, not in the wheel. `vsv.library()` finds it from a
checkout. From an installed package, pass the directory:

```python
from pathlib import Path

shelf = vsv.library(Path("path/to/pyVSV/data/designs"))
```

## Antigenome

`vaccine.antigenome()` splices the coding sequence into NC_001560.1 and returns a `Genome`.
The chassis sequence is the SBOL document `NC_001560.ttl`.
A glycoprotein replacement has five named genes, so `cds` returns the inserted glycoprotein.
The backbone is `https://vsv.bio/backbone/VSV-Indiana`. The reference accession is NC_001560.1.
The papers' infectious clones stay on `Implementation.plasmid`.

A glycoprotein replacement keeps both flanking intergenic dinucleotides on the backbone. An
insert between the glycoprotein and the large polymerase prefixes `CT`. An insert after the
leader suffixes `CT`, so a gene-end is not joined directly to the next gene-start.

A shelf card has no coding sequence, so `antigenome()` raises for those files.

## Command line

```bash
vsv genes data/raw/NC_001560.fasta
vsv genes data/raw/NC_001560.ttl
vsv extract -g g --kind protein data/raw/NC_001560.gb
vsv designs
vsv show rvsv_dg_zebov_kikwit_gp
vsv check data/designs/vsv_indiana
vsv precedents rvsv_dg_zebov_kikwit_gp --locus 3
vsv precedents data/designs/rvsv_gl_zebov_sgp
vsv antigenome kikwit.ttl
```

On the command line the genes are `n` (nucleoprotein), `p` (phosphoprotein), `m` (matrix protein),
`g` (glycoprotein), and `l` (large polymerase). The full command reference (every
subcommand, its filters, and what an empty `precedents` answer means) is in
[docs/cli.md](https://github.com/marpaia/pyVSV/blob/master/docs/cli.md).

## Literature

### Methods the library cites

`replace` and `insert` record the citation for the method they implement, so a written plan carries
the paper its arrangement came from.

| Arrangement | Call | Method paper |
|---|---|---|
| Glycoprotein replacement | `replace(...)` | Garbutt et al., *Journal of Virology*, May 2004, [10.1128/jvi.78.10.5458-5465.2004](https://doi.org/10.1128/jvi.78.10.5458-5465.2004) |
| Gene after the glycoprotein | `insert(..., after="g")` | [10.1128/jvi.71.8.5982-5989.1997](https://doi.org/10.1128/jvi.71.8.5982-5989.1997) |
| Gene after the leader | `insert(..., after="leader")` | [10.1128/jvi.00019-14](https://doi.org/10.1128/jvi.00019-14) |

Each curated design carries its own citations and quoted evidence in its `PROVENANCE.md`.

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
