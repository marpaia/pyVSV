"""A VSV vaccine: ordered genes, one edit, and the papers that report it.

The Indiana chassis is nucleoprotein, phosphoprotein, matrix protein, glycoprotein, and
large polymerase. A replacement puts an antigen where the glycoprotein was. An insert adds an
antigen after the leader or after the glycoprotein. Leader and trailer are written into the
design file and are not modules. Position is the index in ``modules``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from vsv.design.antigenome import SLOT_CITATION, cassette_violations, splice
from vsv.design.classes import antigen_matches, is_leaf
from vsv.design.coding import coding_violations
from vsv.design.evidence import (
    CheckResult,
    Host,
    Implementation,
    Observation,
    Phenotype,
    Timing,
    Violation,
    host_matches,
    is_citation,
)
from vsv.design.modules import Antigen as Antigen
from vsv.design.modules import Edit, Module, _require_id
from vsv.design.modules import Insert as Insert
from vsv.design.modules import Replace as Replace
from vsv.design.modules import Viral as Viral
from vsv.design.sbol.write import write as _write_plan
from vsv.design.vector import GlycoproteinForm as GlycoproteinForm
from vsv.design.vector import Replication as Replication
from vsv.errors import InputError
from vsv.genes import Gene
from vsv.genome.record import Genome

_BACKBONE = "VSV-Indiana"

KNOWLEDGE_GAPS = frozenset({"replication-unknown", "parent-unknown"})


@dataclass(frozen=True, slots=True)
class Vaccine:
    """One genome plan and the implementations and observations recorded for it."""

    id: str
    name: str
    replication: Replication
    modules: tuple[Module, ...]
    parent: str | None = None
    edit: Edit | None = None
    citations: tuple[str, ...] = ()
    aliases: tuple[str, ...] = ()
    implementations: tuple[Implementation, ...] = ()
    observations: tuple[Observation, ...] = ()
    realization: str | None = None

    def __post_init__(self) -> None:
        _require_id(self.id, "vaccine id")
        if not self.name:
            raise InputError(f"Vaccine {self.id} needs a name.")
        if self.parent is not None:
            _require_id(self.parent, "parent id")
        object.__setattr__(self, "modules", tuple(self.modules))
        object.__setattr__(self, "citations", tuple(self.citations))
        object.__setattr__(self, "aliases", tuple(self.aliases))
        object.__setattr__(self, "implementations", tuple(self.implementations))
        object.__setattr__(self, "observations", tuple(self.observations))
        identifiers = [module.id for module in self.modules]
        if len(identifiers) != len(set(identifiers)):
            raise InputError(f"{self.id}: module ids must be unique.")
        if self.aliases and not _alias_ready(self.modules):
            message = f"{self.id}: an alias needs a leaf class, a strain, and a glycoprotein form."
            raise InputError(message)
        known = {item.id for item in self.implementations}
        for observation in self.observations:
            if observation.implementation_id not in known:
                message = (
                    f"{self.id}: observation {observation.id} names unknown "
                    f"implementation {observation.implementation_id}."
                )
                raise InputError(message)

    @property
    def backbone(self) -> str:
        """The chassis this library builds on."""
        return _BACKBONE

    @property
    def antigens(self) -> tuple[Antigen, ...]:
        """Antigen modules in genome order."""
        return tuple(module for module in self.modules if isinstance(module, Antigen))

    def replace(
        self,
        coding: str,
        *,
        identifier: str,
        name: str,
        antigen_class: str,
        replication: Replication,
        strain: str = "",
        accession: str = "",
        form: GlycoproteinForm = GlycoproteinForm.UNSPECIFIED,
        module: str = "Glycoprotein",
        antigen_name: str = "",
    ) -> Vaccine:
        """Replace the VSV glycoprotein with an antigen and cite the 2004 method.

        Citations on the result are that method paper. ``module`` defaults to ``Glycoprotein``.
        ``antigen_name`` defaults to ``name``.
        """
        if not _has_viral_g(self.modules):
            raise InputError(f"{self.id}: {Gene.GLYCOPROTEIN.value} is not on this genome.")
        antigen = _antigen(
            module,
            antigen_name or name,
            antigen_class,
            strain,
            accession,
            form,
            coding,
        )
        modules: list[Module] = [
            antigen if isinstance(item, Viral) and item.gene is Gene.GLYCOPROTEIN else item
            for item in self.modules
        ]
        if sum(item.id == antigen.id for item in modules) != 1:
            raise InputError(f"{self.id}: module {antigen.id} already exists.")
        return Vaccine(
            identifier,
            name,
            replication,
            tuple(modules),
            parent=self.id,
            edit=Replace(antigen.id),
            citations=(SLOT_CITATION["g"],),
        )

    def insert(
        self,
        coding: str,
        *,
        after: Literal["leader", "g"],
        identifier: str,
        name: str,
        antigen_class: str,
        module: str,
        strain: str = "",
        accession: str = "",
        form: GlycoproteinForm = GlycoproteinForm.UNSPECIFIED,
        antigen_name: str = "",
    ) -> Vaccine:
        """Insert an antigen after the leader or after the glycoprotein.

        An insert after the glycoprotein cites the 1997 method. An insert after the leader
        cites the 2014 method. Replication stays the replication of this vaccine.
        """
        if after not in {"leader", "g"}:
            raise InputError(f"{self.id}: a published insertion site is the leader or g.")
        if after == "g" and not _has_viral_g(self.modules):
            raise InputError(f"{self.id}: {Gene.GLYCOPROTEIN.value} is not on this genome.")
        antigen = _antigen(
            module,
            antigen_name or name,
            antigen_class,
            strain,
            accession,
            form,
            coding,
        )
        if any(item.id == antigen.id for item in self.modules):
            raise InputError(f"{self.id}: module {antigen.id} already exists.")
        modules = list(self.modules)
        modules.insert(_insert_index(modules, after), antigen)
        citation = SLOT_CITATION["g-l"] if after == "g" else SLOT_CITATION["leader"]
        return Vaccine(
            identifier,
            name,
            self.replication,
            tuple(modules),
            parent=self.id,
            edit=Insert(after, antigen.id),
            citations=(citation,),
        )

    def check(self) -> CheckResult:
        """Every structural rule and, when a coding sequence is set, every sequence rule.

        The result is true when ``violations`` is empty. A knowledge gap is a violation.
        ``vsv check`` exits 0 when every violation is a gap.
        """
        violations = [
            *_cited(self),
            *_observations(self),
            *_replication(self),
            *_knowledge(self),
            *_sequence(self),
        ]
        return CheckResult(tuple(violations))

    def precedents(
        self,
        vaccines: list[Vaccine] | tuple[Vaccine, ...],
        *,
        antigen_class: str | None = None,
        phenotype: Phenotype | None = None,
        host: Host | None = None,
        strain: str | None = None,
        form: GlycoproteinForm | None = None,
        timing: Timing | None = None,
        locus: int | None = None,
    ) -> list[Vaccine]:
        """Published vaccines with this vaccine's edit.

        A glycoprotein replacement matches other replacements. An insert matches inserts at
        the same site. A vaccine with any number of antigens other than one is not compared.
        ``locus`` is the transcription index of the antigen. Position 1 is 0.
        """
        shape = _shape(self)
        if shape is None:
            return []
        found: list[Vaccine] = []
        for item in vaccines:
            if item.id == self.id or _shape(item) != shape:
                continue
            classes = [antigen.antigen_class for antigen in item.antigens]
            if antigen_class is not None and not any(
                antigen_matches(slug, antigen_class) for slug in classes
            ):
                continue
            if strain is not None and strain not in {antigen.strain for antigen in item.antigens}:
                continue
            if form is not None and form not in {antigen.form for antigen in item.antigens}:
                continue
            if locus is not None and locus not in _loci(item):
                continue
            narrowed = phenotype is not None or host is not None or timing is not None
            matched = any(
                _observation_matches(obs, phenotype, host, timing) for obs in item.observations
            )
            if narrowed and not matched:
                continue
            found.append(item)
        return found

    def antigenome(self) -> Genome:
        """The edited NC_001560.1 antigenome, with its transcription units identified."""
        slot = _slot(self)
        antigen = _edited_antigen(self)
        if antigen is None or not antigen.coding:
            raise InputError(f"{self.id}: the plan has no coding sequence.")
        if slot is None:
            raise InputError(f"{self.id}: the plan has no published insertion site.")
        return splice(antigen.coding, slot, identifier=self.id, description=self.name)

    def write_plan(self, path: str | Path) -> None:
        """Write this vaccine as an SBOL3 design file. Callers do not handle SBOL objects."""
        _write_plan(self, path)

    def write_antigenome(self, path: str | Path) -> None:
        """Write the spliced NC_001560.1 antigenome as an SBOL3 document."""
        self.antigenome().write(path)


def indiana() -> Vaccine:
    """The Indiana chassis, with no edit.

    Nucleoprotein, phosphoprotein, matrix protein, glycoprotein, and large polymerase.
    """
    return Vaccine(
        "vsv_indiana",
        "VSV Indiana",
        Replication.COMPETENT,
        tuple(Viral(gene) for gene in Gene),
        realization="NC_001560",
    )


def _antigen(
    module: str,
    name: str,
    antigen_class: str,
    strain: str,
    accession: str,
    form: GlycoproteinForm,
    coding: str,
) -> Antigen:
    return Antigen(module, name, antigen_class, strain, accession, form, coding)


def _has_viral_g(modules: tuple[Module, ...] | list[Module]) -> bool:
    return any(isinstance(module, Viral) and module.gene is Gene.GLYCOPROTEIN for module in modules)


def _insert_index(modules: list[Module], after: Literal["leader", "g"]) -> int:
    if after == "leader":
        return 0
    for index, module in enumerate(modules):
        if isinstance(module, Viral) and module.gene is Gene.GLYCOPROTEIN:
            return index + 1
    raise InputError(f"{Gene.GLYCOPROTEIN.value} is not on this genome.")


def _alias_ready(modules: tuple[Module, ...]) -> bool:
    antigens = [module for module in modules if isinstance(module, Antigen)]
    if not antigens:
        return False
    return all(
        is_leaf(module.antigen_class)
        and bool(module.strain)
        and module.form is not GlycoproteinForm.UNSPECIFIED
        for module in antigens
    )


def _edited_antigen(vaccine: Vaccine) -> Antigen | None:
    if vaccine.edit is None:
        return None
    return next((module for module in vaccine.antigens if module.id == vaccine.edit.module), None)


def _slot(vaccine: Vaccine) -> str | None:
    edit = vaccine.edit
    if isinstance(edit, Replace):
        return "g"
    if isinstance(edit, Insert) and edit.after == "g":
        return "g-l"
    if isinstance(edit, Insert) and edit.after == "leader":
        return "leader"
    return None


def _shape(vaccine: Vaccine) -> tuple[str, str] | None:
    if len(vaccine.antigens) != 1 or vaccine.edit is None:
        return None
    if isinstance(vaccine.edit, Replace):
        return ("replace", "g")
    if vaccine.edit.after == "g":
        return ("insert", "g")
    return ("insert", "leader")


def _loci(vaccine: Vaccine) -> set[int]:
    return {index for index, module in enumerate(vaccine.modules) if isinstance(module, Antigen)}


def _observation_matches(
    observation: Observation,
    phenotype: Phenotype | None,
    host: Host | None,
    timing: Timing | None,
) -> bool:
    if phenotype is not None and observation.phenotype is not phenotype:
        return False
    if host is not None and not host_matches(observation.host, host):
        return False
    return timing is None or observation.timing is timing


def _cited(vaccine: Vaccine) -> list[Violation]:
    derived = vaccine.parent is not None or vaccine.edit is not None
    if not derived:
        return []
    if any(is_citation(citation) for citation in vaccine.citations):
        return []
    return [Violation("cited", f"{vaccine.id}: a derived vaccine needs a DOI or PMID citation.")]


def _observations(vaccine: Vaccine) -> list[Violation]:
    violations = []
    known_hosts = {host.value for host in Host}
    for index, observation in enumerate(vaccine.observations):
        if not is_citation(observation.citation):
            violations.append(
                Violation(
                    "observation",
                    f"{vaccine.id}: observation {index} needs a DOI or PMID citation.",
                )
            )
        if observation.host.value not in known_hosts:
            violations.append(
                Violation(
                    "observation",
                    f"{vaccine.id}: observation {index} names an unknown host.",
                )
            )
    return violations


def _replication(vaccine: Vaccine) -> list[Violation]:
    viral_g = _has_viral_g(vaccine.modules)
    if vaccine.replication is Replication.SINGLE_CYCLE and viral_g:
        return [
            Violation(
                "replication",
                f"{vaccine.id}: a single-cycle genome still carries {Gene.GLYCOPROTEIN.value}.",
            )
        ]
    if vaccine.replication is not Replication.COMPETENT or viral_g:
        return []
    antigen = _edited_antigen(vaccine)
    if antigen is not None and antigen.form.supplies_entry():
        return []
    return [
        Violation(
            "replication",
            f"{vaccine.id}: a replication-competent genome whose glycoprotein was removed "
            "does not name an entry protein.",
        )
    ]


def _knowledge(vaccine: Vaccine) -> list[Violation]:
    gaps: list[Violation] = []
    if vaccine.replication is Replication.UNKNOWN:
        gaps.append(Violation("replication-unknown", f"{vaccine.id}: replication is unknown."))
    if vaccine.parent is None and vaccine.edit is not None:
        gaps.append(Violation("parent-unknown", f"{vaccine.id}: parent design is unknown."))
    return gaps


def _sequence(vaccine: Vaccine) -> list[Violation]:
    antigen = _edited_antigen(vaccine)
    slot = _slot(vaccine)
    if antigen is None or not antigen.coding or slot is None or vaccine.edit is None:
        return []
    return [
        *cassette_violations(vaccine.id, antigen.coding, slot),
        *coding_violations(
            vaccine.id,
            antigen.coding,
            replication=vaccine.replication,
            form=antigen.form,
            antigen_class=antigen.antigen_class,
            replaced=isinstance(vaccine.edit, Replace),
        ),
    ]
