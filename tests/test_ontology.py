"""The synthetic virology ontology names the IRIs that design files write."""

from pathlib import Path

import rdflib
from rdflib import OWL, RDF, RDFS, URIRef

from vsv.corpus.extract import _REPORTER_SLUGS, _SLUGS
from vsv.design.evidence import Assay, Host, Outcome, Phenotype, Stage
from vsv.design.query import antigen_matches
from vsv.genes import Gene

_ROOT = Path(__file__).resolve().parents[1]
_ONTOLOGY = _ROOT / "ontology" / "synthetic-virology.ttl"
_NS = "https://vsv.bio/"
_SO_GENE = URIRef("http://purl.obolibrary.org/obo/SO_0000704")
_ONTOLOGY_IRI = URIRef("https://vsv.bio/ontology")


def _graph() -> rdflib.Graph:
    graph = rdflib.Graph()
    graph.parse(_ONTOLOGY)
    return graph


def _is_class(graph: rdflib.Graph, iri: str) -> bool:
    return (URIRef(iri), RDF.type, OWL.Class) in graph


def test_ontology_declares_itself_and_imports_the_sequence_ontology() -> None:
    graph = _graph()
    assert (_ONTOLOGY_IRI, RDF.type, OWL.Ontology) in graph
    imports = set(graph.objects(_ONTOLOGY_IRI, OWL.imports))
    assert URIRef("http://purl.obolibrary.org/obo/so.owl") in imports


def test_antigen_classes_and_evidence_labels_are_classes() -> None:
    graph = _graph()
    slugs = {slug for _pattern, slug, _name in _SLUGS if slug not in _REPORTER_SLUGS}
    for slug in slugs:
        assert _is_class(graph, f"{_NS}antigen-class/{slug}"), slug
    for value in (*Host, *Phenotype, *Assay, *Outcome, *Stage):
        prefix = {
            Host: "host",
            Phenotype: "phenotype",
            Assay: "assay",
            Outcome: "outcome",
            Stage: "stage",
        }[type(value)]
        assert _is_class(graph, f"{_NS}{prefix}/{value.value}")
    for gene in Gene:
        assert _is_class(graph, f"{_NS}gene/{gene.value}")
    for name in ("replication-competent", "single-cycle", "leader", "trailer", "reporter"):
        assert _is_class(graph, f"{_NS}role/{name}")
    for name in ("replace-module", "insert-module", "delete-module"):
        assert _is_class(graph, f"{_NS}activity/{name}")


def test_an_antigen_class_is_a_sequence_ontology_gene() -> None:
    graph = _graph()
    antigen = URIRef(f"{_NS}AntigenClass")
    filovirus = URIRef(f"{_NS}antigen-class/filovirus-glycoprotein")
    sudan = URIRef(f"{_NS}antigen-class/sudan-ebolavirus-glycoprotein")
    assert (antigen, RDFS.subClassOf, _SO_GENE) in graph
    assert (filovirus, RDFS.subClassOf, antigen) in graph
    assert (sudan, RDFS.subClassOf, filovirus) in graph


def test_antigen_queries_follow_the_ontology() -> None:
    graph = _graph()
    prefix = f"{_NS}antigen-class/"
    edges = [
        (str(child).removeprefix(prefix), str(parent).removeprefix(prefix))
        for child, parent in graph.subject_objects(RDFS.subClassOf)
        if str(child).startswith(prefix) and str(parent).startswith(prefix)
    ]
    assert ("sudan-ebolavirus-glycoprotein", "filovirus-glycoprotein") in edges
    for child, parent in edges:
        assert antigen_matches(child, parent)
    assert not antigen_matches("marburgvirus-glycoprotein", "filovirus-glycoprotein")


def test_alias_is_an_annotation_and_the_sequence_ontology_is_not_redefined() -> None:
    graph = _graph()
    alias = URIRef(f"{_NS}alias")
    assert (alias, RDF.type, OWL.AnnotationProperty) in graph
    local_gene = URIRef("https://vsv.bio/SO_0000704")
    assert (local_gene, RDF.type, OWL.Class) not in graph
