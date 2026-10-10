# Synthetic virology ontology

`synthetic-virology.ttl` is an OWL ontology for recombinant RNA-virus vector genome plans. Its IRI is `https://vsv.bio/ontology` and the current version IRI is `https://vsv.bio/ontology/2026-10-07`. Term IRIs are the `https://vsv.bio/...` URIs that `design.ttl` already writes.

## What a term is

A genome plan in `data/designs/` is an SBOL component. Gene order is `sbol:strictlyPrecedes`. A feature that is a gene also carries a Sequence Ontology role. Those facts are SBOL and the Sequence Ontology. This file does not redefine them.

The Sequence Ontology has no term for an Indiana backbone, for replication competence, or for replacing VSV G versus inserting a gene downstream of VSV G. Those are the terms defined here. Where a term is a sequence feature, it is a subclass of the Sequence Ontology class with the same meaning:

| This ontology | Sequence Ontology |
|---|---|
| vector genome plan, Indiana backbone, replication-competent, single-cycle | `SO:0001041` viral sequence |
| antigen class, VSV gene | `SO:0000704` gene |
| leader, trailer | `SO:0000001` region |
| gene start | `SO:0005836` regulatory region |
| gene end | `SO:0000141` terminator |

The ontology imports `http://purl.obolibrary.org/obo/so.owl` and uses the OBO IRIs (`http://purl.obolibrary.org/obo/SO_0000704`). Design files write the same terms as `https://identifiers.org/SO:0000704`. The library treats those spellings as one term.

Species glycoprotein classes are subclasses of `filovirus-glycoprotein`, including Zaire and Marburg. A query for the genus class includes those species. Strain and glycoprotein form are properties of an antigen, not classes.

## What stays out

The file does not contain a copy of the Sequence Ontology. Host, route, assay, phenotype, outcome, stage, dose, and timing are strings on a design file. The library reads those strings as its own enums. Glycoprotein form is a property of an antigen. The Zaire class records an editing-run length for the full-length form and the soluble form.

`alias` is a product name for one strain-resolved genome. `strain` and `isolate-accession` name that antigen. `plasmid` and `cloning-sites` name the infectious clone a paper used. `index` is the authored order of an implementation or observation. `index` is not a biological relation.

`day` is days relative to challenge, and a negative day is before challenge. `vaccine-amount` and `challenge-amount` are dose amounts. `challenge-dose` keeps a dose the paper printed as prose.

SBOL and PROV relations are not redefined. A built virus is an `sbol:Implementation` whose `sbol:built` property is the plan. Parentage is `prov:wasDerivedFrom`. The edit is a `prov:Activity` whose `sbol:type` is `replace-module` or `insert-module`.

## Reading it

The Turtle file is the ontology. Protege and any RDF library open it. Imports are followed only when a reasoner is asked to load them. To inspect the local terms without downloading the Sequence Ontology, parse the file and do not resolve `owl:imports`.

The library reads subclass edges and editing-run annotations from this file into a `Vaccine` check. It does not run a reasoner and it does not resolve `owl:imports`. A design file is still read into a `Vaccine` on its own.

```python
import rdflib

graph = rdflib.Graph()
graph.parse("ontology/synthetic-virology.ttl")
```
