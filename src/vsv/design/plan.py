"""A genome plan held as the SBOL component that describes it.

``indiana`` is the canonical chassis. ``replace`` and ``insert`` write a new document: ordered
features, a precedence chain, a derivation activity, and the citation for that method. The plan
file carries no sequence. A plan built from a marked feature keeps that coding sequence in
memory, and ``plasmid`` reads it. Precedents are curated plans with the same precedence chain.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace as data_replace
from pathlib import Path

import sbol3
from lab.part import Part

from vsv import NAMESPACE
from vsv.design.architecture import Architecture, Replication, canonical_chassis, project
from vsv.design.document import architecture_of, design_of, to_document
from vsv.design.evidence import Host, Observation, Phenotype
from vsv.design.library import LIBRARY_ROOT
from vsv.design.modules import Antigen, ViralGene, require_display_id
from vsv.design.plasmid import Plasmid, TranscriptionUnit, backbone
from vsv.design.query import antigen_matches
from vsv.design.rules import CheckResult, check
from vsv.errors import InputError
from vsv.genes import Gene
from vsv.genome.payload import AntigenClass, AntigenFeature
from vsv.genome.sbol import ROLE_LEADER, Genome, display_id, role_key
from vsv.genome.sequence import Nucleotides

_BACKBONE = f"{NAMESPACE}/backbone/"
_ANTIGEN_CLASS = f"{NAMESPACE}/antigen-class/"
_REPORTER = f"{NAMESPACE}/role/reporter"
_REPLACE_GLYCOPROTEIN = "https://doi.org/10.1128/jvi.78.10.5458-5465.2004"
_INSERT_AFTER_GLYCOPROTEIN = "https://doi.org/10.1128/jvi.71.8.5982-5989.1997"
_INSERT_AFTER_LEADER = "https://doi.org/10.1128/jvi.00019-14"


class Plan:
    """One genome plan and the SBOL document that holds it.

    ``sequences`` maps an antigen module id to the coding sequence recorded when that module was
    marked. ``write`` does not store those sequences. A plan read from a file has none.
    """

    def __init__(
        self,
        document: sbol3.Document,
        component: sbol3.Component,
        *,
        sequences: Mapping[str, Nucleotides] | None = None,
    ) -> None:
        self.document = document
        self.component = component
        self._sequences = dict(sequences or {})

    @property
    def id(self) -> str:
        """The plan's display id."""
        return str(self.component.display_id)

    @property
    def leader(self) -> sbol3.LocalSubComponent:
        """The leader feature."""
        found = [feature for feature in self._features() if _has_role(feature, ROLE_LEADER)]
        if len(found) != 1:
            raise InputError(f"{self.id}: expected one leader.")
        return found[0]

    @property
    def observations(self) -> tuple[Observation, ...]:
        """Observations stored on this document."""
        return design_of(self.document).observations

    def feature(self, gene: Gene) -> sbol3.LocalSubComponent:
        """The feature whose display id is ``gene``."""
        found = [feature for feature in self._features() if feature.display_id == gene.value]
        if len(found) != 1:
            raise InputError(f"{self.id}: expected one {gene.value} feature.")
        return found[0]

    def antigen(self, antigen_class: AntigenClass) -> sbol3.LocalSubComponent:
        """The feature that carries ``antigen_class``."""
        found = [
            feature for feature in self._features() if _antigen_slug(feature) == antigen_class.slug
        ]
        if len(found) != 1:
            raise InputError(
                f"{self.id}: expected one {antigen_class.slug} feature, found {len(found)}."
            )
        return found[0]

    def write(self, path: str | Path) -> None:
        """Write this plan as an SBOL3 file. The suffix selects the format."""
        self.document.write(path)

    def replace(
        self,
        target: sbol3.LocalSubComponent,
        payload: AntigenFeature,
        *,
        replication: Replication,
        identifier: str,
        name: str,
    ) -> Plan:
        """Replace the viral glycoprotein with ``payload``.

        ``replication`` is competent when the payload supplies entry, and single-cycle when it
        does not. The new document cites the 2004 glycoprotein-replacement method. The new plan
        keeps the payload's coding sequence.
        """
        self._require_owned(target)
        if not _is_viral_glycoprotein(target):
            raise InputError(f"{self.id}: {Gene.GLYCOPROTEIN} is not a viral gene on this plan.")
        antigen = _require_sequence(self.id, payload)
        built = architecture_of(self.document).replace(
            str(target.display_id),
            antigen,
            identifier,
            name=name,
            replication=replication,
        )
        sequences = dict(self._sequences)
        sequences[antigen.id] = payload.sequence
        return _publish(built, _REPLACE_GLYCOPROTEIN, sequences)

    def insert(
        self,
        payload: AntigenFeature,
        *,
        after: sbol3.LocalSubComponent,
        identifier: str,
        name: str,
    ) -> Plan:
        """Insert ``payload`` after the leader or after the viral glycoprotein.

        An insert after the glycoprotein keeps that gene and cites the 1997 method. An insert
        after the leader cites the 2014 method. The new plan keeps this plan's replication class
        and the payload's coding sequence.
        """
        self._require_owned(after)
        if _is_leader(after):
            citation = _INSERT_AFTER_LEADER
        elif _is_viral_glycoprotein(after):
            citation = _INSERT_AFTER_GLYCOPROTEIN
        else:
            raise InputError(
                f"{self.id}: a published insertion site is the leader or {Gene.GLYCOPROTEIN}."
            )
        antigen = _require_sequence(self.id, payload)
        built = architecture_of(self.document).insert(
            str(after.display_id), antigen, identifier, name=name
        )
        sequences = dict(self._sequences)
        sequences[antigen.id] = payload.sequence
        return _publish(built, citation, sequences)

    def precedents(
        self,
        plans: Sequence[Plan],
        *,
        antigen_class: AntigenClass | None = None,
        phenotype: Phenotype | None = None,
        host: Host | None = None,
    ) -> list[Plan]:
        """Published plans that arrange a gene the same way this one does.

        A glycoprotein replacement returns the curated replacements, an insert after G returns
        the curated inserts after G, and an insert after the leader returns the curated insert
        there. ``antigen_class`` keeps plans whose antigen is that class or a subclass of it.
        ``phenotype`` and ``host`` keep plans with a matching observation. Each result is that
        published plan, observations included. This candidate keeps its own empty observation
        list.
        """
        shape = _shape(self)
        if shape is None:
            return []
        found: list[Plan] = []
        for item in plans:
            if item.component.identity == self.component.identity or _shape(item) != shape:
                continue
            if antigen_class is not None and not _matches_antigen(item, antigen_class):
                continue
            if phenotype is not None or host is not None:
                observations = item.observations
                if not any(
                    (phenotype is None or observation.phenotype is phenotype)
                    and (host is None or observation.host is host)
                    for observation in observations
                ):
                    continue
            found.append(item)
        return found

    def plasmid(self) -> Plasmid:
        """The plasmid for this plan, using the coding sequence it already holds.

        A replacement uses the glycoprotein slot on ``pvsv-dg``. An insert after G uses the
        site between G and L. An insert after the leader uses the site there. The sequence is
        a transcription unit: the slot's gene-start, this coding sequence, and the gene-end.
        """
        shape = _shape(self)
        if shape is None:
            raise InputError(f"{self.id}: the plan has no published insertion site.")
        if not self._sequences:
            raise InputError(f"{self.id}: the plan has no coding sequence.")
        if len(self._sequences) != 1:
            raise InputError(f"{self.id}: the plan has more than one coding sequence.")
        insert_id, coding = next(iter(self._sequences.items()))
        edit, slot_name, backbone_id = _PLASMID_SITES[shape]
        cited = backbone(backbone_id)
        unit = TranscriptionUnit(cited.slot.gene_start, coding, cited.slot.gene_end)
        sequence = unit.sequence
        antigenome = Nucleotides(
            cited.sequence.bases[: cited.slot.start]
            + sequence.bases
            + cited.sequence.bases[cited.slot.end :]
        )
        return Plasmid(
            plan_id=self.id,
            edit=edit,
            slot=slot_name,
            insert_id=insert_id,
            insert=Part(f"{NAMESPACE}/{self.id}/{insert_id}"),
            unit=unit,
            sequence=sequence,
            replication=architecture_of(self.document).replication,
            backbone=cited,
            antigenome=antigenome,
            product=Part(f"{NAMESPACE}/{self.id}/plasmid"),
        )

    def check(self) -> CheckResult:
        """Whether this plan meets the structural rules.

        The result is true when it passes. ``violations`` lists every failure.
        """
        record = design_of(self.document)
        return CheckResult(tuple(check(record.architecture, record.observations)))

    def _features(self) -> list[sbol3.LocalSubComponent]:
        return [
            feature
            for feature in self.component.features
            if isinstance(feature, sbol3.LocalSubComponent)
        ]

    def _require_owned(self, feature: sbol3.LocalSubComponent) -> None:
        owned = any(str(item.identity) == str(feature.identity) for item in self._features())
        if not owned:
            raise InputError(f"{self.id}: the feature is not on this plan.")


def indiana() -> Plan:
    """The Indiana chassis: leader, N, P, M, G, L, and trailer."""
    document = to_document(canonical_chassis())
    return Plan(document, _plan_component(document))


def plan(genome: Genome, *, backbone: str = "VSV-Indiana") -> Plan:
    """The sequence-free plan realized by an identified five-gene genome.

    The genome's sequence stays on ``genome``. The returned document has no sequence.
    """
    document = to_document(project(genome, backbone=backbone))
    return Plan(document, _plan_component(document))


def read(path: str | Path) -> Plan:
    """Read the plan in an SBOL3 file."""
    document = sbol3.Document()
    try:
        document.read(path)
    except (OSError, ValueError) as error:
        raise InputError(f"{path}: {error}") from error
    return Plan(document, _plan_component(document))


def plans(root: Path | None = None) -> list[Plan]:
    """Every curated plan under ``root``, in directory-name order."""
    directory = LIBRARY_ROOT if root is None else root
    if not directory.is_dir():
        raise InputError(f"{directory}: design library not found.")
    found: list[Plan] = []
    for path in sorted(item for item in directory.iterdir() if item.is_dir()):
        turtle = path / "design.ttl"
        if not turtle.is_file():
            raise InputError(f"{path}: missing design.ttl.")
        found.append(read(turtle))
    identifiers = [item.id for item in found]
    if len(identifiers) != len(set(identifiers)):
        raise InputError(f"{directory}: duplicate design ids.")
    return found


def _publish(built: Architecture, citation: str, sequences: Mapping[str, Nucleotides]) -> Plan:
    document = to_document(data_replace(built, citations=(citation,)))
    return Plan(document, _plan_component(document), sequences=sequences)


def _require_sequence(plan_id: str, payload: AntigenFeature) -> Antigen:
    if not payload.sequence.bases:
        raise InputError(f"{plan_id}: {payload.feature.display_id} has no coding sequence.")
    return _antigen(payload)


def _matches_antigen(plan: Plan, antigen_class: AntigenClass) -> bool:
    slugs = architecture_of(plan.document).antigen_classes
    return any(antigen_matches(slug, antigen_class.slug) for slug in slugs)


_G = Gene.GLYCOPROTEIN.short
_PLASMID_SITES = {
    "replace-glycoprotein": (f"replace {_G}", _G, "pvsv-dg"),
    "insert-after-glycoprotein": (f"insert after {_G}", f"{_G}-l", "pvsv-g-l"),
    "insert-after-leader": ("insert after leader", "leader", "pvsv-leader"),
}


def _plan_component(document: sbol3.Document) -> sbol3.Component:
    found = [
        obj
        for obj in document.objects
        if isinstance(obj, sbol3.Component)
        and any(str(role).startswith(_BACKBONE) for role in obj.roles)
    ]
    if len(found) != 1:
        raise InputError(f"Expected one plan, found {len(found)}.")
    return found[0]


def _antigen(payload: AntigenFeature) -> Antigen:
    raw = str(payload.feature.display_id)
    if raw.endswith("_cds"):
        raw = raw[: -len("_cds")]
    module_id = display_id(raw)
    require_display_id(module_id, "antigen")
    name = str(payload.feature.name or module_id)
    return Antigen(module_id, name, payload.antigen_class.slug)


def _has_role(feature: sbol3.Feature, role: str) -> bool:
    return any(role_key(str(item)) == role_key(role) for item in feature.roles)


def _antigen_slug(feature: sbol3.Feature) -> str | None:
    matches = [str(role) for role in feature.roles if str(role).startswith(_ANTIGEN_CLASS)]
    if len(matches) > 1:
        raise InputError(f"{feature.display_id}: more than one antigen class.")
    if not matches:
        return None
    return matches[0].rsplit("/", 1)[-1]


def _is_leader(feature: sbol3.LocalSubComponent) -> bool:
    return _has_role(feature, ROLE_LEADER)


def _is_viral_glycoprotein(feature: sbol3.LocalSubComponent) -> bool:
    return (
        feature.display_id == Gene.GLYCOPROTEIN
        and _antigen_slug(feature) is None
        and not _has_role(feature, _REPORTER)
    )


def _shape(plan: Plan) -> str | None:
    """The published edit this precedence chain states, when it states exactly one."""
    architecture = architecture_of(plan.document)
    if any(module.locus is None for module in architecture.transcribed):
        return None
    body = list(architecture.transcribed)
    extra = [index for index, module in enumerate(body) if not isinstance(module, ViralGene)]
    if len(extra) != 1:
        return None
    index = extra[0]
    viral_g = any(
        isinstance(module, ViralGene) and module.id == Gene.GLYCOPROTEIN for module in body
    )
    if (
        not viral_g
        and 0 < index < len(body) - 1
        and isinstance(body[index - 1], ViralGene)
        and body[index - 1].id == Gene.MATRIX_PROTEIN
        and isinstance(body[index + 1], ViralGene)
        and body[index + 1].id == Gene.LARGE_POLYMERASE
    ):
        return "replace-glycoprotein"
    if (
        viral_g
        and index > 0
        and isinstance(body[index - 1], ViralGene)
        and body[index - 1].id == Gene.GLYCOPROTEIN
    ):
        return "insert-after-glycoprotein"
    if (
        index == 0
        and len(body) > 1
        and isinstance(body[1], ViralGene)
        and body[1].id == Gene.NUCLEOPROTEIN
    ):
        return "insert-after-leader"
    return None
