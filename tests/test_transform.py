"""Architecture transformations record a parent and leave the source alone."""

from collections.abc import Callable

import pytest

from vsv.design import (
    Antigen,
    Architecture,
    Delete,
    Insert,
    Leader,
    ModuleKind,
    Replace,
    Replication,
    Reporter,
    Trailer,
    ViralGene,
    canonical_chassis,
)
from vsv.design.rules import Property, satisfies
from vsv.errors import InputError
from vsv.genes import Gene

_ANTIGEN = Antigen("antigen_1", "filovirus glycoprotein", "filovirus-glycoprotein", locus=99)
_REPORTER = Reporter("reporter_1", "GFP", locus=99)


def test_replace_insert_and_delete_preserve_canonical_order() -> None:
    source = canonical_chassis()
    replaced = source.replace(Gene.GLYCOPROTEIN, _ANTIGEN, "replaced")
    inserted = source.insert(Gene.NUCLEOPROTEIN, _REPORTER, "inserted")
    deleted = source.delete(Gene.LARGE_POLYMERASE, "deleted")
    for design in (replaced, inserted, deleted):
        assert design.parent_id == source.id
        assert satisfies(design, (), Property.CANONICAL_ORDER)
    assert [(module.id, module.locus) for module in source.modules] == [
        ("leader", None),
        ("Nucleoprotein", 0),
        ("Phosphoprotein", 1),
        ("MatrixProtein", 2),
        ("Glycoprotein", 3),
        ("LargePolymerase", 4),
        ("trailer", None),
    ]


def test_replace_records_the_edit_and_ignores_the_supplied_locus() -> None:
    source = canonical_chassis()
    result = source.replace(
        Gene.GLYCOPROTEIN, _ANTIGEN, "rvsv_dg", replication=Replication.SINGLE_CYCLE
    )
    antigen = next(module for module in result.modules if module.id == "antigen_1")
    assert antigen.locus == 3
    assert antigen.kind is ModuleKind.ANTIGEN
    assert antigen.feature_identity is None
    assert Gene.GLYCOPROTEIN not in {module.id for module in result.modules}
    assert result.replication is Replication.SINGLE_CYCLE
    assert result.edit == Replace("Glycoprotein", "antigen_1", Replication.SINGLE_CYCLE)
    assert source.replication is Replication.COMPETENT
    assert source.edit is None


def test_insert_renumbers_later_loci() -> None:
    result = canonical_chassis().insert(Gene.NUCLEOPROTEIN, _REPORTER, "with_reporter")
    assert [(module.id, module.locus) for module in result.modules] == [
        ("leader", None),
        ("Nucleoprotein", 0),
        ("reporter_1", 1),
        ("Phosphoprotein", 2),
        ("MatrixProtein", 3),
        ("Glycoprotein", 4),
        ("LargePolymerase", 5),
        ("trailer", None),
    ]
    assert result.edit == Insert("Nucleoprotein", "reporter_1")


def test_delete_renumbers_and_can_record_replication() -> None:
    result = canonical_chassis().delete(
        Gene.GLYCOPROTEIN, "delta_g", replication=Replication.SINGLE_CYCLE
    )
    assert [module.id for module in result.modules] == [
        "leader",
        "Nucleoprotein",
        "Phosphoprotein",
        "MatrixProtein",
        "LargePolymerase",
        "trailer",
    ]
    assert [module.locus for module in result.modules if module.locus is not None] == [0, 1, 2, 3]
    assert result.edit == Delete("Glycoprotein", Replication.SINGLE_CYCLE)


@pytest.mark.parametrize(
    ("message", "apply"),
    [
        ("no module", lambda source: source.replace("X", _ANTIGEN, "child")),
        ("cannot be replaced", lambda source: source.replace("leader", _ANTIGEN, "child")),
        (
            "already exists",
            lambda source: source.insert(
                "Phosphoprotein", ViralGene("Nucleoprotein", "Nucleoprotein", locus=0), "child"
            ),
        ),
        ("insertions", lambda source: source.insert("trailer", _REPORTER, "child")),
        ("cannot be deleted", lambda source: source.delete("leader", "child")),
        ("own id", lambda source: source.delete("LargePolymerase", "vsv_indiana")),
    ],
)
def test_edits_that_are_rejected(
    message: str, apply: Callable[[Architecture], Architecture]
) -> None:
    with pytest.raises(InputError, match=message):
        apply(canonical_chassis())


def test_the_last_transcribed_module_cannot_be_deleted() -> None:
    solo = Architecture(
        "solo",
        "solo",
        "VSV",
        Replication.COMPETENT,
        (
            Leader(),
            ViralGene("Nucleoprotein", "Nucleoprotein", locus=0),
            Trailer(),
        ),
    )
    with pytest.raises(InputError, match="last transcribed"):
        solo.delete("Nucleoprotein", "empty")
