"""The SBOL3 design file is the shareable record of one vaccine."""

from pathlib import Path

import sbol3

import vsv
from vsv import GlycoproteinForm, Replication, indiana, library
from vsv.genome.record import Genome, Span
from vsv.genome.sbol import write_genome
from vsv.sbol import NAMESPACE, ROLE_LEADER, SO_REGION

_CODING = "ATG" + "GCT" * 20 + "TAA"
_CLASS = "https://vsv.bio/antigen-class/zaire-ebolavirus-glycoprotein"
_GENE = "https://identifiers.org/SO:0000704"


def _ordered(component: sbol3.Component) -> list[str]:
    """Display ids in the strictly-precedes chain, including leader and trailer."""
    features = [
        feature for feature in component.features if isinstance(feature, sbol3.LocalSubComponent)
    ]
    by_identity = {str(feature.identity): feature for feature in features}
    outgoing: dict[str, str] = {}
    incoming: set[str] = set()
    for constraint in component.constraints:
        if str(constraint.restriction) != sbol3.SBOL_STRICTLY_PRECEDES:
            continue
        subject, target = str(constraint.subject), str(constraint.object)
        outgoing[subject] = target
        incoming.add(target)
    starts = [feature for feature in features if str(feature.identity) not in incoming]
    assert len(starts) == 1
    ordered: list[str] = []
    current: sbol3.Feature | None = starts[0]
    seen: set[str] = set()
    while current is not None:
        identity = str(current.identity)
        assert identity not in seen
        seen.add(identity)
        ordered.append(str(current.display_id))
        nxt = outgoing.get(identity)
        current = by_identity[nxt] if nxt is not None else None
    assert len(ordered) == len(features)
    return ordered


def test_a_design_file_orders_the_genome_and_records_its_parent(tmp_path: Path) -> None:
    vaccine = indiana().replace(
        _CODING,
        identifier="rvsv_dg_outbreak_gp",
        name="rVSV outbreak glycoprotein",
        antigen_class="zaire-ebolavirus-glycoprotein",
        replication=Replication.COMPETENT,
        strain="Kikwit-1995",
        form=GlycoproteinForm.FULL_LENGTH,
    )
    path = tmp_path / "design.ttl"
    vaccine.write(path)
    document = sbol3.Document()
    document.read(str(path))
    (component,) = (obj for obj in document.objects if isinstance(obj, sbol3.Component))
    assert str(component.identity) == "https://vsv.bio/rvsv_dg_outbreak_gp"
    assert "https://vsv.bio/vsv_indiana" in [str(item) for item in component.derived_from]
    assert _ordered(component) == [
        "leader",
        "Nucleoprotein",
        "Phosphoprotein",
        "MatrixProtein",
        "Glycoprotein",
        "LargePolymerase",
        "trailer",
    ]
    antigen = next(
        feature for feature in component.features if str(feature.display_id) == "Glycoprotein"
    )
    roles = {str(role) for role in antigen.roles}
    assert _CLASS in roles
    assert _GENE in roles
    assert any(isinstance(obj, sbol3.Sequence) for obj in document.objects)


def test_a_second_shelf_round_trip_writes_no_sequence(tmp_path: Path) -> None:
    kikwit = next(item for item in library() if item.id == "rvsv_dg_zebov_kikwit_gp")
    path = tmp_path / "design.ttl"
    kikwit.write(path)
    vsv.read(path).write(path)
    document = sbol3.Document()
    document.read(str(path))
    assert "hasSequence" not in path.read_text()
    assert not any(isinstance(obj, sbol3.Sequence) for obj in document.objects)


def test_a_plan_and_a_genome_share_the_leader_role(tmp_path: Path) -> None:
    design = tmp_path / "design.ttl"
    indiana().write(design)
    genome = tmp_path / "genome.ttl"
    genome.write_text(write_genome(Genome("tiny", "tiny genome", "ACGTACGT", leader=Span(0, 4))))
    assert _leader_roles(design) == {SO_REGION, ROLE_LEADER}
    assert _leader_roles(genome) == {SO_REGION, ROLE_LEADER}


def _leader_roles(path: Path) -> set[str]:
    sbol3.set_namespace(NAMESPACE)
    document = sbol3.Document()
    document.read(str(path))
    found: set[str] = set()
    for component in document.objects:
        if not isinstance(component, sbol3.Component):
            continue
        assert str(component.identity).startswith(f"{NAMESPACE}/")
        for feature in component.features:
            roles = {str(role) for role in feature.roles}
            if ROLE_LEADER in roles:
                found = roles
    assert found
    return found
