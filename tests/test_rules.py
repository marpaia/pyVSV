"""Design rules for architectures and the observations attached to them."""

from dataclasses import dataclass, replace

from vsv.design import (
    Antigen,
    Architecture,
    Leader,
    Module,
    ParentStatus,
    Replace,
    Replication,
    Trailer,
    canonical_chassis,
)
from vsv.design.modules import at_locus
from vsv.design.rules import Property, check, satisfies

_DOI = "https://doi.org/10.1126/science.aar7204"
_PMID = "PMID:36508878"


@dataclass(frozen=True, slots=True)
class Observed:
    citation: str
    host: str


def _without_g(architecture: Architecture) -> Architecture:
    locus = 0
    modules: list[Module] = []
    for module in architecture.modules:
        if module.id == "Glycoprotein":
            continue
        if isinstance(module, Leader | Trailer):
            modules.append(module)
        else:
            modules.append(at_locus(module, locus))
            locus += 1
    return replace(architecture, modules=tuple(modules))


def _rules(architecture: Architecture, *observations: Observed) -> set[str]:
    return {violation.rule for violation in check(architecture, observations)}


def test_the_canonical_chassis_meets_every_structural_rule() -> None:
    assert check(canonical_chassis()) == []
    assert satisfies(canonical_chassis(), (), Property.CANONICAL_ORDER)
    assert satisfies(canonical_chassis(), (), Property.CITED)
    assert not satisfies(canonical_chassis(), (), Property.PRECLINICAL_OR_ABOVE)


def test_viral_genes_out_of_order() -> None:
    architecture = canonical_chassis()
    modules = []
    for module in architecture.modules:
        if module.id == "Nucleoprotein":
            modules.append(replace(module, id="Phosphoprotein", name="Phosphoprotein"))
        elif module.id == "Phosphoprotein":
            modules.append(replace(module, id="Nucleoprotein", name="Nucleoprotein"))
        else:
            modules.append(module)
    swapped = replace(architecture, modules=tuple(modules))
    assert _rules(swapped) == {"viral-gene-order"}
    assert not satisfies(swapped, (), Property.CANONICAL_ORDER)


def test_transcribed_loci_must_be_contiguous() -> None:
    architecture = canonical_chassis()
    modules = [
        at_locus(module, 5) if module.id == "LargePolymerase" else module
        for module in architecture.modules
    ]
    assert "loci" in _rules(replace(architecture, modules=tuple(modules)))


def test_the_leader_is_first_and_the_trailer_is_last() -> None:
    architecture = canonical_chassis()
    leader, *rest = architecture.modules
    moved = replace(architecture, modules=(*rest, leader))
    assert _rules(moved) == {"termini"}


def test_an_antigen_without_a_stated_locus_is_a_knowledge_gap() -> None:
    architecture = canonical_chassis()
    antigen = Antigen("antigen_1", "filovirus glycoprotein", "filovirus-glycoprotein")
    modules = (*architecture.modules[:-1], antigen, architecture.modules[-1])
    rules = _rules(replace(architecture, modules=modules))
    assert "locus-unknown" in rules
    assert "loci" not in rules


def test_a_derived_architecture_needs_a_resolvable_citation() -> None:
    child = replace(canonical_chassis(), id="child", parent="vsv_indiana")
    assert _rules(child) == {"cited"}
    assert not satisfies(child, (), Property.CITED)
    unreadable = replace(child, citations=("not-a-paper",))
    assert "cited" in _rules(unreadable)
    cited = replace(child, citations=(_DOI,))
    assert "cited" not in _rules(cited)
    assert satisfies(cited, (), Property.CITED)
    assert "cited" not in _rules(replace(child, citations=(_PMID,)))


def test_observations_need_a_citation_and_a_known_host() -> None:
    architecture = canonical_chassis()
    assert "observation" in _rules(architecture, Observed("not-a-paper", "mouse"))
    assert "observation" in _rules(architecture, Observed(_PMID, "mars"))
    assert check(architecture, (Observed(_PMID, "mouse"),)) == []


def test_preclinical_evidence_is_a_rodent_or_higher_host() -> None:
    architecture = canonical_chassis()
    rodent = (Observed(_PMID, "hamster"),)
    culture = (Observed(_PMID, "cell-line"),)
    assert satisfies(architecture, rodent, Property.PRECLINICAL_OR_ABOVE)
    assert satisfies(architecture, (Observed(_DOI, "human"),), Property.PRECLINICAL_OR_ABOVE)
    assert not satisfies(architecture, culture, Property.PRECLINICAL_OR_ABOVE)


def test_single_cycle_vectors_do_not_keep_viral_gene_g() -> None:
    single_cycle = replace(canonical_chassis(), replication=Replication.SINGLE_CYCLE)
    assert "replication" in _rules(single_cycle)
    antigen = Antigen("antigen_1", "filovirus glycoprotein", "filovirus-glycoprotein", locus=3)
    modules = tuple(
        antigen if module.id == "Glycoprotein" else module for module in canonical_chassis().modules
    )
    replaced = Architecture(
        "rvsv_dg",
        "rVSV",
        "VSV-Indiana",
        Replication.SINGLE_CYCLE,
        modules,
        parent="vsv_indiana",
        citations=(_DOI,),
        edit=Replace("Glycoprotein", "antigen_1", Replication.SINGLE_CYCLE),
    )
    assert "replication" not in _rules(replaced)


def test_unknown_replication_and_parent_are_gaps() -> None:
    architecture = replace(
        canonical_chassis(),
        id="partial",
        replication=Replication.UNKNOWN,
        parent=ParentStatus.UNKNOWN,
        citations=(_DOI,),
    )
    assert _rules(architecture) == {"replication-unknown", "parent-unknown"}
    named = replace(architecture, parent="vsv_indiana")
    assert named.parent_id == "vsv_indiana"
    assert "parent-unknown" not in _rules(named)


def test_a_competent_vector_keeps_a_g_locus() -> None:
    missing = _without_g(canonical_chassis())
    assert "replication" in _rules(missing)
    antigen = Antigen("antigen_1", "filovirus glycoprotein", "filovirus-glycoprotein", locus=3)
    modules = tuple(
        antigen if module.id == "Glycoprotein" else module for module in canonical_chassis().modules
    )
    replaced = Architecture(
        "competent_swap",
        "swap",
        "VSV-Indiana",
        Replication.COMPETENT,
        modules,
        parent="vsv_indiana",
        citations=(_DOI,),
        edit=Replace("Glycoprotein", "antigen_1"),
    )
    assert "replication" not in _rules(replaced)
