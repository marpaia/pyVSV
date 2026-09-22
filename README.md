# Vesicular Stomatitis Virus (VSV)

Vesicular Stomatitis Virus (VSV) is a negative-strand RNA virus that belongs to the Rhabdoviridae family. It primarily affects livestock, but VSV is also of interest in virology research due to its potential use as a viral vector for vaccine development and gene therapy. This repository contains information and resources related to VSV and a software ecosystem for the architecture, design, and simulation of VSV-based viral vectors.

<p align="center">
  <img src="./docs/images/vsv-structure-genome.jpg" alt="VSV Structure and Genome">
</p>
<p align="center">
  <em>Liu G, Cao W, Salawudeen A, et al. Vesicular Stomatitis Virus: From Agricultural Pathogen to Vaccine Vector. Pathogens. 2021;10(9):1092. Published 2021 Aug 27. doi:<a href="https://doi.org/10.3390/pathogens10091092">10.3390/pathogens10091092</a></em>
</p>

## Usage

VSV transcribes five genes, always in the order Nucleoprotein, Phosphoprotein, MatrixProtein, Glycoprotein, LargePolymerase. `vsv` finds them in a genome and gives typed access to each gene's coordinates and sequences. FASTA genomes are annotated from their sequence. GenBank and SBOL3 records keep their curated annotations, and their genes are named by position.

### Python

```python
import vsv
from vsv import Gene, Product

genome = vsv.load("data/raw/NC_001560.fasta")  # FASTA, GenBank, or SBOL3

for gene, annotation in genome.genes.items():  # genome order
    print(gene, annotation.gene_span, annotation.cds_span, annotation.protein)

n = genome.genes[Gene.NUCLEOPROTEIN]
n.gene_span  # (50, 1376): zero-based, half-open, so genome.bases[50:1376]
n.cds_dna  # coding sequence with its stop codon
n.sequence(Product.PROTEIN)  # the same as n.protein
```

`vsv.load_all(paths)` reads several files, and a genome whose five genes cannot all be identified raises a `vsv.GeneWarning`. Spans and sequences are `None` for anything a genome does not annotate.

### Command line

```bash
vsv genes data/raw/NC_001560.fasta                        # table of the five genes
vsv extract -g p --kind protein data/raw/NC_001560.gb  # FASTA; --kind is gene, cds, or protein
vsv annotate data/raw/NC_001560.fasta --format gff3       # every feature, also SBOL3 formats
vsv summary data/raw/*.fasta                              # length and base composition
vsv compare data/raw/MF196237.fasta data/raw/NC_001560.fasta
```

`extract` writes all five genes unless `-g` is given, once per gene wanted. On the command line the genes are `n`, `p`, `m`, `g`, and `l`. Genomes that lack a gene get a warning on standard error.

### Designs

`data/designs/` is a curated library of VSV architectures, the implementations of those architectures, and the qualitative evidence attached to each implementation. Antigen modules are named by class. A locus, a parent, or a replication class is unknown when the source does not state it. The library does not store insert sequences.

```bash
vsv designs
vsv designs --antigen-class henipavirus-glycoprotein
vsv show rvsv_niv_g
vsv show data/designs/rvsv_lasv_gpc
vsv check data/designs/vsv_indiana
```

`vsv check` prints structural violations and knowledge gaps. An unknown locus, replication class, or parent is a gap. The command exits with status 1 when a structural rule fails.

`vsv ingest` searches Europe PMC for recombinant VSV vector papers and writes a local review queue under `data/ingest/`. Proposals stay pending until a person promotes one into `data/designs/`. Long nucleotide runs are omitted before anything is saved. The queue is not part of the curated library.

## Literature

- [Vesicular stomatitis virus: re-inventing the bullet](https://pubmed.ncbi.nlm.nih.gov/15121047/) - May 2004
- [Vesicular stomatitis virus (VSV) as a vaccine vector for immunization against viral infections](https://pubmed.ncbi.nlm.nih.gov/24379275/) - Jan 2013
- [Complete Genome Sequences of Two Vesicular Stomatitis Virus Isolates Collected in Mexico](https://pmc.ncbi.nlm.nih.gov/articles/PMC5597772/) - Sep 2017 ([VSV Isolate 1](https://www.ncbi.nlm.nih.gov/nuccore/MF196237), [VSV Isolate 2](https://www.ncbi.nlm.nih.gov/nuccore/MF196238))
- [Vesicular Stomatitis Virus: From Agricultural Pathogen to Vaccine Vector](https://pubmed.ncbi.nlm.nih.gov/34578125/) - Aug 2021
- [An efficient plasmid-based system for the recovery of recombinant vesicular stomatitis virus encoding foreign glycoproteins](https://www.nature.com/articles/s41598-024-65384-8) - Jun 2024
- [Leveraging Synthetic Virology for the Rapid Engineering of Vesicular Stomatitis Virus (VSV)](https://pubmed.ncbi.nlm.nih.gov/39459973/) - Oct 2024
- [Low-dose VSV-EBOV vaccination provides rapid protection from lethal Ebola virus challenge](https://pubmed.ncbi.nlm.nih.gov/41727071/) - Feb 2026

## Data

- [VSV Genome Sequence](https://www.ncbi.nlm.nih.gov/nuccore/NC_001560) - Complete genome sequence of Vesicular Stomatitis Virus (VSV)
- [VSV Protein Structures](https://www.rcsb.org/structure/6VSV) - 3D structures of VSV proteins available in the Protein Data Bank (PDB)
- [Vesicuovirus on ViralZone](https://viralzone.expasy.org/21) - Information on the Vesiculovirus genus, including VSV, on ViralZone

## Resources

- [VSV Gene Expression Vector on VectorBuilder](https://en.vectorbuilder.com/resources/vector-system/pVSV_Exp.html)
