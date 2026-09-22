# Synthetic virology ontology

`synthetic-virology.ttl` is an OWL ontology for recombinant RNA-virus vector genome plans. Its IRI is `https://vsv.bio/ontology` and the current version IRI is `https://vsv.bio/ontology/2026-09-22`. Term IRIs are the `https://vsv.bio/...` URIs that `design.ttl` already writes.

## What a term is

A genome plan in `data/designs/` is an SBOL component. Gene order is `sbol:strictlyPrecedes`. A feature that is a gene also carries a Sequence Ontology role. Those facts are SBOL and the Sequence Ontology. This file does not redefine them.

The Sequence Ontology has no term for an Indiana backbone, for replication competence, or for replacing VSV G versus inserting a gene downstream of VSV G. Those are the terms defined here. Where a term is a sequence feature, it is a subclass of the Sequence Ontology class with the same meaning:

| This ontology | Sequence Ontology |
|---|---|
| vector genome plan, Indiana backbone, New Jersey backbone, replication-competent, single-cycle | `SO:0001041` viral sequence |
| antigen class, VSV gene, reporter | `SO:0000704` gene |
| leader, trailer | `SO:0000001` region |
| gene start | `SO:0005836` regulatory region |
| gene end | `SO:0000141` terminator |

The ontology imports `http://purl.obolibrary.org/obo/so.owl` and uses the OBO IRIs (`http://purl.obolibrary.org/obo/SO_0000704`). Design files write the same terms as `https://identifiers.org/SO:0000704`. The library treats those spellings as one term.

Species glycoprotein classes are subclasses of `filovirus-glycoprotein`. A query for the genus class therefore includes Sudan, Bundibugyo, Reston, and Taï Forest. The Zaire plans in the library use the genus class itself.

## What stays out

The file does not contain a copy of the Sequence Ontology. It does not define assays, hosts, or vaccine products as a replacement for the Ontology for Biomedical Investigations or the Vaccine Ontology. Phenotype, host, assay, outcome, and development stage are included because a design file stores those labels, and a reader of the file needs a definition. They are the library's labels.

`alias` is the literature name of a plan. `index` is the authored order of an implementation or observation. `index` is not a biological relation.

SBOL and PROV relations are not redefined. A built virus is an `sbol:Implementation` whose `sbol:built` property is the plan. Parentage is `prov:wasDerivedFrom`. The edit is a `prov:Activity` whose `sbol:type` is `replace-module`, `insert-module`, or `delete-module`.

## Reading it

The Turtle file is the ontology. Protege and any RDF library open it. Imports are followed only when a reasoner is asked to load them. To inspect the local terms without downloading the Sequence Ontology, parse the file and do not resolve `owl:imports`.

```python
import rdflib

graph = rdflib.Graph()
graph.parse("ontology/synthetic-virology.ttl")
```
