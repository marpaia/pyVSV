"""Architecture transformations record a parent and leave the source alone."""

import pytest

from vsv.architecture import (
    Architecture,
    Derivation,
    Module,
    ModuleKind,
    Replication,
    canonical_chassis,
)
from vsv.genes import Gene
from vsv.rules import Property, satisfies
from vsv.sequence import InputError
from vsv.transform import DeleteModule, InsertModule, ReplaceModule

_ANTIGEN = Module(
    "antigen_1",
    ModuleKind.ANTIGEN,
    "filovirus glycoprotein",
    locus=99,
    antigen_class="filovirus-glycoprotein",
)
_REPORTER = Module("reporter_1", ModuleKind.REPORTER, "GFP", locus=99)


def test_replace_insert_and_delete_preserve_canonical_order() -> None:
    source = canonical_chassis()
    replaced = ReplaceModule(Gene.GLYCOPROTEIN, _ANTIGEN, "replaced").apply(source)
    inserted = InsertModule(Gene.NUCLEOPROTEIN, _REPORTER, "inserted").apply(source)
    deleted = DeleteModule(Gene.LARGE_POLYMERASE, "deleted").apply(source)
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
    result = ReplaceModule(
        Gene.GLYCOPROTEIN, _ANTIGEN, "rvsv_dg", replication=Replication.SINGLE_CYCLE
    ).apply(source)
    antigen = next(module for module in result.modules if module.id == "antigen_1")
    assert antigen.locus == 3
    assert antigen.kind is ModuleKind.ANTIGEN
    assert antigen.feature_identity is None
    assert Gene.GLYCOPROTEIN not in {module.id for module in result.modules}
    assert result.replication is Replication.SINGLE_CYCLE
    assert result.derivation == Derivation(
        "replace-module", "Glycoprotein", "antigen_1", Replication.SINGLE_CYCLE
    )
    assert source.replication is Replication.COMPETENT
    assert source.derivation is None


def test_insert_renumbers_later_loci() -> None:
    result = InsertModule(Gene.NUCLEOPROTEIN, _REPORTER, "with_reporter").apply(canonical_chassis())
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
    assert result.derivation == Derivation("insert-module", "Nucleoprotein", "reporter_1")


def test_delete_renumbers_and_can_record_replication() -> None:
    result = DeleteModule(Gene.GLYCOPROTEIN, "delta_g", replication=Replication.SINGLE_CYCLE).apply(
        canonical_chassis()
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
    assert result.derivation == Derivation(
        "delete-module", "Glycoprotein", replication=Replication.SINGLE_CYCLE
    )


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        (ReplaceModule("X", _ANTIGEN, "child"), "no module"),
        (ReplaceModule("leader", _ANTIGEN, "child"), "cannot be replaced"),
        (
            InsertModule(
                "Phosphoprotein",
                Module("Nucleoprotein", ModuleKind.VIRAL_GENE, "Nucleoprotein", locus=0),
                "child",
            ),
            "already exists",
        ),
        (InsertModule("trailer", _REPORTER, "child"), "insertions"),
        (DeleteModule("leader", "child"), "cannot be deleted"),
        (DeleteModule("LargePolymerase", "vsv_indiana"), "own id"),
    ],
)
def test_edits_that_are_rejected(
    edit: ReplaceModule | InsertModule | DeleteModule, message: str
) -> None:
    with pytest.raises(InputError, match=message):
        edit.apply(canonical_chassis())


def test_the_last_transcribed_module_cannot_be_deleted() -> None:
    solo = Architecture(
        "solo",
        "solo",
        "VSV",
        Replication.COMPETENT,
        (
            Module("leader", ModuleKind.LEADER, "leader"),
            Module("Nucleoprotein", ModuleKind.VIRAL_GENE, "Nucleoprotein", locus=0),
            Module("trailer", ModuleKind.TRAILER, "trailer"),
        ),
    )
    with pytest.raises(InputError, match="last transcribed"):
        DeleteModule("Nucleoprotein", "empty").apply(solo)
