"""The synthetic virology ontology names the IRIs that design files write."""

from pathlib import Path

import rdflib
from rdflib import OWL, RDF, RDFS, URIRef

from vsv.design.classes import antigen_matches, editing_length
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


def test_antigen_classes_are_classes() -> None:
    graph = _graph()
    prefix = f"{_NS}antigen-class/"
    filovirus = URIRef(prefix + "filovirus-glycoprotein")
    children = {
        str(child).removeprefix(prefix) for child in graph.subjects(RDFS.subClassOf, filovirus)
    }
    assert children
    for slug in {*children, "filovirus-glycoprotein"}:
        assert _is_class(graph, f"{prefix}{slug}"), slug
    for gene in Gene:
        assert _is_class(graph, f"{_NS}gene/{gene.value}")
    for name in ("replication-competent", "single-cycle", "leader", "trailer"):
        assert _is_class(graph, f"{_NS}role/{name}")
    for name in ("replace-module", "insert-module"):
        assert _is_class(graph, f"{_NS}activity/{name}")
    absent = (
        f"{_NS}host/mouse",
        f"{_NS}role/reporter",
        f"{_NS}activity/delete-module",
        f"{_NS}backbone/VSV-NewJersey",
        f"{_NS}tail/native",
    )
    for iri in absent:
        assert not _is_class(graph, iri), iri


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
    for child, parent in edges:
        assert antigen_matches(child, parent)
        assert antigen_matches(child, child)
    assert antigen_matches("zaire-ebolavirus-glycoprotein", "filovirus-glycoprotein")
    assert antigen_matches("marburgvirus-glycoprotein", "filovirus-glycoprotein")
    assert not antigen_matches("arenavirus-gpc", "filovirus-glycoprotein")
    assert not antigen_matches("zaire-ebolavirus-glycoprotein", "sudan-ebolavirus-glycoprotein")
    assert not antigen_matches("influenza-hemagglutinin", "filovirus-glycoprotein")


def test_the_zaire_class_declares_an_editing_run() -> None:
    assert editing_length("zaire-ebolavirus-glycoprotein", "full-length") == 8
    assert editing_length("zaire-ebolavirus-glycoprotein", "soluble") == 7
    assert editing_length("zaire-ebolavirus-glycoprotein", "mucin-deleted") is None
    assert editing_length("sudan-ebolavirus-glycoprotein", "full-length") is None


def test_alias_is_an_annotation_and_the_sequence_ontology_is_not_redefined() -> None:
    graph = _graph()
    alias = URIRef(f"{_NS}alias")
    assert (alias, RDF.type, OWL.AnnotationProperty) in graph
    local_gene = URIRef("https://vsv.bio/SO_0000704")
    assert (local_gene, RDF.type, OWL.Class) not in graph
