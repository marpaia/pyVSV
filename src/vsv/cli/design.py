"""Commands for the curated VSV vaccine library."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import click

from vsv.cli.render import cli_name, indent, row, table, widths, wrap
from vsv.design.evidence import (
    Challenge,
    DoseUnit,
    EvidenceLevel,
    Host,
    Observation,
    Phenotype,
    Route,
    Timing,
    evidence_level,
    host_matches,
)
from vsv.design.library import LIBRARY_ROOT, library, load, select
from vsv.design.vaccine import (
    KNOWLEDGE_GAPS,
    Antigen,
    GlycoproteinForm,
    Insert,
    Replace,
    Vaccine,
    Viral,
)
from vsv.errors import InputError
from vsv.genome.sbol import write_genome

_EVIDENCE = {
    "none": EvidenceLevel.NONE,
    "in-vitro": EvidenceLevel.IN_VITRO,
    "rodent": EvidenceLevel.RODENT,
    "preclinical": EvidenceLevel.RODENT,
    "nonhuman-primate": EvidenceLevel.NONHUMAN_PRIMATE,
    "human": EvidenceLevel.HUMAN,
}
_FORMS = (
    GlycoproteinForm.FULL_LENGTH,
    GlycoproteinForm.SOLUBLE,
    GlycoproteinForm.MUCIN_DELETED,
    GlycoproteinForm.TAIL_CHIMERA,
)
_TIMINGS = (Timing.PROPHYLACTIC, Timing.POST_EXPOSURE)


class EvidenceName(click.ParamType[EvidenceLevel]):
    """An evidence threshold, including the ``preclinical`` alias for rodent."""

    name = "evidence"

    def convert(
        self, value: object, param: click.Parameter | None, ctx: click.Context | None
    ) -> EvidenceLevel:
        if isinstance(value, EvidenceLevel):
            return value
        text = str(value)
        try:
            return _EVIDENCE[text]
        except KeyError:
            self.fail(f"unknown evidence {text!r}.", param, ctx)


def _phenotype(value: str | None) -> Phenotype | None:
    return None if value is None else Phenotype(value)


def _host(value: str | None) -> Host | None:
    return None if value is None else Host(value)


def _form(value: str | None) -> GlycoproteinForm | None:
    return None if value is None else GlycoproteinForm(value)


def _timing(value: str | None) -> Timing | None:
    return None if value is None else Timing(value)


def _resolve(identifier: str) -> Vaccine:
    path = Path(identifier)
    if path.is_dir() or path.is_file():
        return load(path)
    matches = [vaccine for vaccine in library() if vaccine.id == identifier]
    if len(matches) != 1:
        raise InputError(f"{identifier}: design not found.")
    return matches[0]


def _design_rows(records: Sequence[Vaccine]) -> str:
    rows: list[list[str]] = []
    for record in records:
        rows.append(
            [
                record.id,
                _derivation(record) if record.edit is not None else "-",
                ",".join(antigen.antigen_class for antigen in record.antigens) or "-",
                record.replication.value,
                evidence_level([item.host for item in record.observations]).value,
                _gap_text(record),
            ]
        )
    if not rows:
        return "no designs match"
    return table(["id", "edit", "antigens", "replication", "evidence", "gaps"], rows)


def _gap_text(vaccine: Vaccine) -> str:
    names = [item.rule for item in vaccine.check().violations if item.rule in KNOWLEDGE_GAPS]
    return ",".join(names) or "-"


def _parent_text(vaccine: Vaccine) -> str:
    if vaccine.parent is not None:
        return vaccine.parent
    if vaccine.edit is None:
        return "root"
    return "unknown"


def _derivation(vaccine: Vaccine) -> str:
    edit = vaccine.edit
    if isinstance(edit, Replace):
        text = "replace g"
        shown = cli_name(edit.module)
        if shown != "g":
            text += f" with {shown}"
    elif isinstance(edit, Insert):
        text = f"insert after {edit.after} with {cli_name(edit.module)}"
    else:
        text = "-"
    return f"{text}, {vaccine.replication.value}"


def _fact_lines(pairs: Sequence[tuple[str, str]]) -> list[str]:
    width = max(len(key) for key, _value in pairs)
    return [f"{key:<{width}}  {value}" for key, value in pairs]


def _module_table(vaccine: Vaccine) -> str:
    rows = []
    for index, module in enumerate(vaccine.modules):
        rows.append(
            [
                str(index),
                cli_name(module.id),
                module.antigen_class if isinstance(module, Antigen) else "",
                *_antigen_fields(module),
            ]
        )
    return table(["locus", "module", "class", "strain", "form"], rows)


def _antigen_fields(module: Viral | Antigen) -> tuple[str, str]:
    if not isinstance(module, Antigen):
        return ("-", "-")
    strain = module.strain or "-"
    form = "-" if module.form is GlycoproteinForm.UNSPECIFIED else module.form.value
    return (strain, form)


def _dose_text(amount: str, unit: DoseUnit) -> str:
    if unit is DoseUnit.UNSPECIFIED:
        return amount
    return f"{amount} {unit.value}".strip()


def _challenge_text(challenge: Challenge) -> str:
    parts = [challenge.virus or "-"]
    if challenge.strain:
        parts.append(challenge.strain)
    if challenge.amount or challenge.unit is not DoseUnit.UNSPECIFIED:
        parts.append(_dose_text(challenge.amount, challenge.unit))
    if challenge.route is not Route.UNKNOWN:
        parts.append(challenge.route.value)
    return "  ".join(parts)


def _observation_lines(observations: Sequence[Observation], citation: str) -> list[str]:
    if not observations:
        return ["  no observations"]
    header = ["phenotype", "host", "assay", "outcome"]
    rows = [
        [item.phenotype.value, item.host.value, item.assay.value, item.outcome.value]
        for item in observations
    ]
    column_widths = widths(header, rows)
    lines = ["  " + row(header, column_widths)]
    for item, cells in zip(observations, rows, strict=True):
        lines.append("  " + row(cells, column_widths))
        if item.citation != citation:
            lines.append(f"    {item.citation}")
        if item.note:
            lines.extend(wrap(item.note, "    "))
        lines.extend(_observation_details(item))
        lines.append("")
    return lines[:-1]


def _observation_details(item: Observation) -> list[str]:
    details: list[str] = []
    if item.timing is not Timing.UNKNOWN:
        details.append(f"timing  {item.timing.value}")
    if item.day is not None:
        details.append(f"day  {item.day}")
    if item.vaccine_amount or item.vaccine_unit is not DoseUnit.UNSPECIFIED:
        details.append(f"vaccine  {_dose_text(item.vaccine_amount, item.vaccine_unit)}")
    if item.challenge is not None:
        details.append(f"challenge  {_challenge_text(item.challenge)}")
    return [f"    {line}" for line in details]


def _show(vaccine: Vaccine) -> str:
    facts = [
        ("backbone", vaccine.backbone),
        ("replication", vaccine.replication.value),
        ("parent", _parent_text(vaccine)),
    ]
    if vaccine.edit is not None:
        facts.append(("derivation", _derivation(vaccine)))
    if vaccine.realization:
        facts.append(("realization", vaccine.realization))
    lines = [
        vaccine.id,
        vaccine.name,
        "",
        *_fact_lines(facts),
        "",
        "A gene closer to the 3' end is transcribed more, and position 1 attenuates the virus.",
        "",
        "modules",
        indent(_module_table(vaccine)),
    ]
    grouped: dict[str, list[Observation]] = {item.id: [] for item in vaccine.implementations}
    for observation in vaccine.observations:
        grouped[observation.implementation_id].append(observation)
    if not vaccine.implementations:
        lines.extend(["", "implementations", "  none"])
    for implementation in vaccine.implementations:
        lines.extend(
            [
                "",
                f"{implementation.stage.value}  {implementation.id}",
                f"  {implementation.citation}",
            ]
        )
        if implementation.plasmid:
            lines.append(f"  plasmid  {implementation.plasmid}")
        if implementation.note:
            lines.extend(wrap(implementation.note, "  "))
        lines.append("")
        lines.extend(_observation_lines(grouped[implementation.id], implementation.citation))
    gaps = [item.message for item in vaccine.check().violations if item.rule in KNOWLEDGE_GAPS]
    if gaps:
        lines.extend(["", "gaps", *[f"  {gap}" for gap in gaps]])
    return "\n".join(lines)


@click.command()
@click.option("--antigen-class", default=None, help="Antigen class slug.")
@click.option(
    "--phenotype",
    type=click.Choice([item.value for item in Phenotype]),
    default=None,
    help="Observation phenotype.",
)
@click.option(
    "--host",
    type=click.Choice([item.value for item in Host]),
    default=None,
    help="Observation host. nonhuman-primate includes cynomolgus and rhesus.",
)
@click.option("--strain", default=None, help="Antigen strain.")
@click.option(
    "--form",
    type=click.Choice([item.value for item in _FORMS]),
    default=None,
    help="Glycoprotein form.",
)
@click.option(
    "--timing",
    type=click.Choice([item.value for item in _TIMINGS]),
    default=None,
    help="Vaccine timing relative to challenge.",
)
@click.option(
    "--evidence",
    type=EvidenceName(),
    default=None,
    help="Minimum evidence: in-vitro, rodent, preclinical, nonhuman-primate, or human.",
)
@click.option(
    "--library",
    "library_root",
    type=click.Path(path_type=Path),
    default=LIBRARY_ROOT,
    show_default=True,
    help="Design library.",
)
def designs(
    antigen_class: str | None,
    phenotype: str | None,
    host: str | None,
    strain: str | None,
    form: str | None,
    timing: str | None,
    evidence: EvidenceLevel | None,
    library_root: Path,
) -> None:
    """List curated VSV vaccines."""
    records = select(
        library(library_root),
        antigen_class=antigen_class,
        phenotype=_phenotype(phenotype),
        host=_host(host),
        min_evidence=evidence,
        strain=strain,
        form=_form(form),
        timing=_timing(timing),
    )
    click.echo(_design_rows(records))


@click.command()
@click.argument("design")
@click.pass_context
def check(ctx: click.Context, design: str) -> None:
    """Report design-rule violations for one curated vaccine."""
    vaccine = _resolve(design)
    violations = vaccine.check().violations
    structural = [item for item in violations if item.rule not in KNOWLEDGE_GAPS]
    gaps = [item for item in violations if item.rule in KNOWLEDGE_GAPS]
    if not structural and not gaps:
        click.echo(f"{vaccine.id}: ok")
        return
    blocks = []
    if structural:
        blocks.append(
            table(["rule", "message"], [(item.rule, item.message) for item in structural])
        )
    if gaps:
        blocks.append(table(["gap", "message"], [(item.rule, item.message) for item in gaps]))
    click.echo("\n".join(blocks))
    if structural:
        ctx.exit(1)


@click.command()
@click.argument("plan")
@click.option("--antigen-class", default=None, help="Antigen class slug, including subclasses.")
@click.option(
    "--phenotype",
    type=click.Choice([item.value for item in Phenotype]),
    default=None,
    help="Observation phenotype.",
)
@click.option(
    "--host",
    type=click.Choice([item.value for item in Host]),
    default=None,
    help="Observation host. nonhuman-primate includes cynomolgus and rhesus.",
)
@click.option("--strain", default=None, help="Antigen strain.")
@click.option(
    "--form",
    type=click.Choice([item.value for item in _FORMS]),
    default=None,
    help="Glycoprotein form.",
)
@click.option(
    "--timing",
    type=click.Choice([item.value for item in _TIMINGS]),
    default=None,
    help="Vaccine timing relative to challenge.",
)
@click.option(
    "--locus",
    type=int,
    default=None,
    help="Transcription index of the non-viral gene. Position 1 is 0.",
)
def precedents(
    plan: str,
    antigen_class: str | None,
    phenotype: str | None,
    host: str | None,
    strain: str | None,
    form: str | None,
    timing: str | None,
    locus: int | None,
) -> None:
    """List published vaccines with this vaccine's arrangement."""
    candidate = _resolve(plan)
    rows: list[list[str]] = []
    for published in candidate.precedents(
        library(),
        antigen_class=antigen_class,
        phenotype=_phenotype(phenotype),
        host=_host(host),
        strain=strain,
        form=_form(form),
        timing=_timing(timing),
        locus=locus,
    ):
        for result in published.observations:
            if phenotype is not None and result.phenotype is not Phenotype(phenotype):
                continue
            if host is not None and not host_matches(result.host, Host(host)):
                continue
            if timing is not None and result.timing is not Timing(timing):
                continue
            rows.append(
                [
                    published.id,
                    result.phenotype.value,
                    result.host.value,
                    result.assay.value,
                    result.citation,
                ]
            )
    if not rows:
        click.echo("no precedents")
        return
    click.echo(table(["id", "phenotype", "host", "assay", "citation"], rows))


@click.command()
@click.argument("design")
def antigenome(design: str) -> None:
    """Write the spliced NC_001560.1 antigenome as SBOL."""
    click.echo(write_genome(_resolve(design).antigenome()).rstrip("\n"))


@click.command()
@click.argument("design")
def show(design: str) -> None:
    """Show one curated vaccine."""
    click.echo(_show(_resolve(design)))
