# Command line

`vsv` with no subcommand prints help and exits 0. `vsv --version` prints the version. A usage error exits 2. A rejected input prints `vsv: error: …` and exits 1. A gene-naming warning prints `vsv: warning: …` and exits 0.

On the command line the genes are `n`, `p`, `m`, `g`, and `l`.

## Read a genome

One or more FASTA, GenBank, or SBOL genome records. `vsv genes` and `vsv annotate` print tables.
`vsv extract` prints SBOL Turtle.

```sh
# The n, p, m, g, and l genes and coding sequences. A naming problem is a warning.
vsv genes data/raw/NC_001560.fasta
vsv genes data/raw/NC_001560.ttl

# Leader, genes, signals, coding sequences, and trailer.
vsv annotate data/raw/NC_001560.gb

# SBOL for the named genes. Repeat -g; the default is all five. --kind is gene, cds, or protein (default cds).
vsv extract -g g --kind protein data/raw/NC_001560.gb
```

## Inspect the design library

Each curated vaccine is the Indiana chassis or a filovirus vaccine, at `data/designs/<id>/` with `design.ttl` and `PROVENANCE.md`. The Turtle describes the ordered genome and the papers. A nucleotide sequence is present only when the antigen feature has an `sbol:Sequence`. The shelf cards have none.

```sh
# Curated vaccines. Columns: id, edit, antigens, replication, evidence, gaps.
# No match prints "no designs match".
vsv designs

# Antigen class, including subclasses. A Zaire query does not return Sudan or Marburg.
# A filovirus query returns every antigen card and not the Indiana chassis.
vsv designs --antigen-class filovirus-glycoprotein

# Vaccines with a matching observation. Phenotype: expression, growth, immunogenicity, safety, protection.
# Host: cell-line, mouse, hamster, guinea-pig, ferret, nonhuman-primate, cynomolgus, rhesus, human, other.
# A nonhuman-primate query includes cynomolgus and rhesus. --strain matches an antigen strain.
# --form is full-length, soluble, mucin-deleted, or tail-chimera. --timing is prophylactic or post-exposure.
vsv designs --phenotype protection --host nonhuman-primate

# Minimum evidence: in-vitro, rodent, nonhuman-primate, or human. preclinical means rodent.
vsv designs --evidence human

# One vaccine by id, directory, or design.ttl: backbone, replication, parent, derivation, realization,
# modules, strain, form, implementations, plasmid name, observations, timing, day, vaccine dose,
# challenge, and gaps.
# Show states that a gene closer to the 3' end is transcribed more.
vsv show rvsv_dg_zebov_kikwit_gp

# The same view from a design directory. An unknown id exits 1.
vsv show data/designs/rvsv_dg_zebov_kikwit_gp

# Structural violations exit 1. A gap alone (parent-unknown, replication-unknown) exits 0.
vsv check data/designs/vsv_indiana
vsv check rvsv_dg_bombali_gp
```

## Ask which published vaccines match an arrangement

```sh
# Same arrangement: G replaced, one gene after G, or one gene after the leader.
# Omits the query vaccine. --antigen-class includes subclasses.
vsv precedents data/designs/rvsv_dg_zebov_mayinga_gp \
  --antigen-class filovirus-glycoprotein \
  --phenotype protection \
  --host nonhuman-primate

# Columns: id, phenotype, host, assay, citation.
# A vaccine with no observations is omitted from the table. Bombali is one of those.
vsv precedents data/designs/rvsv_dg_zebov_mayinga_gp

# The soluble glycoprotein is the only curated insert after G, so this query returns no precedents.
vsv precedents data/designs/rvsv_gl_zebov_sgp

# "no precedents": the Indiana chassis has no added gene, the transcription index matches nothing,
# or a filter has no stored observation.
vsv precedents data/designs/vsv_indiana

# --strain, --form, and --timing filter the same way as vsv designs.
# --locus is the transcription index of the non-viral gene. Position 1 is 0.
# A G replacement sits at 3. Locus 0 matches nothing on Kikwit.
vsv precedents data/designs/rvsv_dg_zebov_kikwit_gp --locus 3
vsv precedents data/designs/rvsv_dg_zebov_kikwit_gp --locus 0
vsv precedents data/designs/rvsv_gl_zebov_sgp --locus 4
```

## Read an antigenome

```sh
# SBOL Turtle of the spliced NC_001560.1 antigenome for a design id, directory, or design.ttl
# that has a coding sequence. A paper's infectious clone is recorded on the implementation.
# A shelf card has no coding sequence, so the command exits 1.
vsv antigenome kikwit.ttl
```
