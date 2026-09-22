"""Load the curated VSV design library.

Each design is a directory of JSON. The files describe architecture, citations, and qualitative
observations. They do not carry nucleotide or amino-acid sequences.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from enum import StrEnum
from pathlib import Path

from vsv.design.architecture import Architecture, ParentStatus, Replication
from vsv.design.edits import Delete, Edit, Insert, Replace
from vsv.design.evidence import (
    Assay,
    Design,
    Host,
    Implementation,
    Observation,
    Outcome,
    Phenotype,
    Stage,
    Unmodeled,
)
from vsv.design.modules import (
    Antigen,
    Leader,
    LocusStatus,
    Module,
    ModuleKind,
    Reporter,
    Trailer,
    ViralGene,
)
from vsv.errors import InputError

LIBRARY_ROOT = Path(__file__).resolve().parents[3] / "data" / "designs"

_SEQUENCE_KEYS = frozenset(
    {
        "sequence",
        "sequences",
        "cds",
        "protein",
        "bases",
        "nucleotides",
        "amino_acids",
        "translation",
    }
)
_DESIGN_KEYS = frozenset(
    {
        "id",
        "name",
        "backbone",
        "replication",
        "modules",
        "parent_id",
        "parent_status",
        "citations",
        "derivation",
        "realization",
        "model_source",
        "implementations",
        "unmodeled",
    }
)
_MODULE_KEYS = frozenset({"id", "kind", "name", "locus", "locus_status", "antigen_class"})
_DERIVATION_KEYS = frozenset({"activity", "target", "module", "replication"})
_OBSERVATION_KEYS = frozenset(
    {"id", "implementation_id", "phenotype", "host", "assay", "outcome", "citation", "note"}
)
_IMPLEMENTATION_KEYS = frozenset({"id", "stage", "citation", "note"})
_UNMODELED_KEYS = frozenset({"key", "value", "source"})


def _reject_sequences(value: object, path: Path) -> None:
    if isinstance(value, dict):
        forbidden = _SEQUENCE_KEYS.intersection(value)
        if forbidden:
            names = ", ".join(sorted(forbidden))
            raise InputError(f"{path}: sequence payload is not allowed ({names}).")
        for child in value.values():
            _reject_sequences(child, path)
    elif isinstance(value, list):
        for child in value:
            _reject_sequences(child, path)


def _as_text(value: object) -> str:
    if isinstance(value, str):
        return value
    return json.dumps(value)


def _unmodeled_item(value: object, path: Path) -> Unmodeled:
    raw = value if isinstance(value, dict) else None
    if raw is None or not {"key", "value"}.issubset(raw) or set(raw).difference(_UNMODELED_KEYS):
        raise InputError(f"{path}: unmodeled entries need key, value, and an optional source.")
    if str(raw["key"]) in _SEQUENCE_KEYS:
        raise InputError(f"{path}: sequence payload is not allowed ({raw['key']}).")
    source = "" if raw.get("source") is None else str(raw["source"])
    return Unmodeled(str(raw["key"]), _as_text(raw["value"]), source)


def _split(
    value: object, path: Path, allowed: frozenset[str]
) -> tuple[Mapping[str, object], tuple[Unmodeled, ...]]:
    if not isinstance(value, dict):
        raise InputError(f"{path}: expected a JSON object.")
    extra: list[Unmodeled] = []
    declared = value.get("unmodeled", [])
    if not isinstance(declared, list):
        raise InputError(f"{path}: unmodeled must be a list.")
    extra.extend(_unmodeled_item(item, path) for item in declared)
    stray = set(value).difference(allowed, {"unmodeled"})
    for key in sorted(stray):
        extra.append(Unmodeled(key, _as_text(value[key]), ""))
    kept = {key: item for key, item in value.items() if key in allowed}
    return kept, tuple(extra)


def _load_json(path: Path) -> object:
    try:
        value = json.loads(path.read_text())
    except json.JSONDecodeError as error:
        raise InputError(f"{path}: {error}") from error
    _reject_sequences(value, path)
    return value


def _req_str(raw: Mapping[str, object], key: str, path: Path) -> str:
    value = raw.get(key)
    if not isinstance(value, str):
        raise InputError(f"{path}: {key} must be a string.")
    return value


def _opt_str(raw: Mapping[str, object], key: str, path: Path) -> str | None:
    if key not in raw or raw[key] is None:
        return None
    value = raw[key]
    if not isinstance(value, str):
        raise InputError(f"{path}: {key} must be a string.")
    return value


def _opt_int(raw: Mapping[str, object], key: str, path: Path) -> int | None:
    if key not in raw or raw[key] is None:
        return None
    value = raw[key]
    if isinstance(value, bool) or not isinstance(value, int):
        raise InputError(f"{path}: {key} must be an integer.")
    return value


def _enum[E: StrEnum](raw: Mapping[str, object], key: str, kind: type[E], path: Path) -> E:
    text = _req_str(raw, key, path)
    try:
        return kind(text)
    except ValueError as error:
        raise InputError(f"{path}: unknown {key} {text!r}.") from error


def _locus_status(raw: Mapping[str, object], path: Path) -> LocusStatus | None:
    if "locus_status" not in raw or raw["locus_status"] is None:
        return None
    text = raw["locus_status"]
    if not isinstance(text, str):
        raise InputError(f"{path}: locus_status must be a string.")
    try:
        return LocusStatus(text)
    except ValueError as error:
        raise InputError(f"{path}: unknown locus status.") from error


def _check_locus(
    module_id: str,
    kind: ModuleKind,
    locus: int | None,
    status: LocusStatus | None,
    antigen_class: str | None,
) -> None:
    terminal = kind in {ModuleKind.LEADER, ModuleKind.TRAILER}
    if terminal:
        if locus is not None or status not in {None, LocusStatus.NOT_TRANSCRIBED}:
            raise InputError(f"Module {module_id} is not transcribed and has no locus.")
        if antigen_class is not None:
            raise InputError(f"Module {module_id} is not an antigen and has no antigen class.")
        return
    if antigen_class is not None and kind is not ModuleKind.ANTIGEN:
        raise InputError(f"Module {module_id} is not an antigen and has no antigen class.")
    if kind is ModuleKind.ANTIGEN and antigen_class is None:
        raise InputError(f"Module {module_id} needs an antigen class slug.")
    if status is LocusStatus.NOT_TRANSCRIBED:
        raise InputError(f"Module {module_id} is not transcribed and has no locus.")
    if status is LocusStatus.UNKNOWN and locus is not None:
        raise InputError(f"Module {module_id}: an unknown locus has no index.")
    if status is LocusStatus.KNOWN and locus is None:
        raise InputError(f"Module {module_id}: a known locus needs an index.")


def _module(value: object, path: Path) -> tuple[Module, tuple[Unmodeled, ...]]:
    raw, extra = _split(value, path, _MODULE_KEYS)
    kind = _enum(raw, "kind", ModuleKind, path)
    module_id = _req_str(raw, "id", path)
    name = _req_str(raw, "name", path)
    locus = _opt_int(raw, "locus", path)
    status = _locus_status(raw, path)
    antigen_class = _opt_str(raw, "antigen_class", path)
    _check_locus(module_id, kind, locus, status, antigen_class)
    match kind:
        case ModuleKind.LEADER:
            built: Module = Leader(name, module_id)
        case ModuleKind.TRAILER:
            built = Trailer(name, module_id)
        case ModuleKind.VIRAL_GENE:
            built = ViralGene(module_id, name, locus)
        case ModuleKind.REPORTER:
            built = Reporter(module_id, name, locus)
        case ModuleKind.ANTIGEN:
            if antigen_class is None:
                raise InputError(f"Module {module_id} needs an antigen class slug.")
            built = Antigen(module_id, name, antigen_class, locus)
    return built, extra


def _edit(value: object, path: Path) -> Edit | None:
    if value is None:
        return None
    raw, _extra = _split(value, path, _DERIVATION_KEYS)
    activity = _req_str(raw, "activity", path)
    target = _req_str(raw, "target", path)
    module = _opt_str(raw, "module", path)
    replication_text = _opt_str(raw, "replication", path)
    try:
        replication = None if replication_text is None else Replication(replication_text)
    except ValueError as error:
        raise InputError(f"{path}: invalid derivation.") from error
    match activity:
        case "delete-module":
            if module is not None:
                raise InputError(f"{path}: a deletion derivation has no module.")
            return Delete(target, replication)
        case "replace-module":
            if module is None:
                raise InputError(f"{path}: a replace-module derivation names the resulting module.")
            return Replace(target, module, replication)
        case "insert-module":
            if module is None:
                raise InputError(f"{path}: an insert-module derivation names the resulting module.")
            return Insert(target, module, replication)
        case _:
            raise InputError(f"{path}: invalid derivation.")


def _parent(raw: Mapping[str, object], path: Path) -> str | ParentStatus:
    parent_id = _opt_str(raw, "parent_id", path)
    status = _opt_str(raw, "parent_status", path)
    if parent_id is None and status in {None, ParentStatus.ROOT}:
        return ParentStatus.ROOT
    if parent_id is None and status == ParentStatus.UNKNOWN:
        return ParentStatus.UNKNOWN
    if parent_id is not None and status in {None, "known"}:
        return parent_id
    if status not in {None, "known", ParentStatus.ROOT, ParentStatus.UNKNOWN}:
        raise InputError(f"{path}: unknown parent status.")
    raise InputError(f"{path}: parent id and parent status disagree.")


def _implementation(value: object, path: Path) -> tuple[Implementation, tuple[Unmodeled, ...]]:
    raw, extra = _split(value, path, _IMPLEMENTATION_KEYS)
    built = Implementation(
        _req_str(raw, "id", path),
        _enum(raw, "stage", Stage, path),
        _req_str(raw, "citation", path),
        _opt_str(raw, "note", path) or "",
    )
    return built, extra


def _architecture(
    path: Path,
) -> tuple[Architecture, str | None, str | None, tuple[Implementation, ...], tuple[Unmodeled, ...]]:
    raw, unmodeled = _split(_load_json(path), path, _DESIGN_KEYS)
    modules = raw.get("modules")
    if not isinstance(modules, list):
        raise InputError(f"{path}: modules must be a list.")
    citations = raw.get("citations", [])
    if not isinstance(citations, list) or not all(isinstance(item, str) for item in citations):
        raise InputError(f"{path}: citations must be a list of strings.")
    implementations = raw.get("implementations", [])
    if not isinstance(implementations, list):
        raise InputError(f"{path}: implementations must be a list.")
    parsed_modules = [_module(module, path) for module in modules]
    parsed_implementations = [_implementation(item, path) for item in implementations]
    unmodeled = (
        unmodeled
        + tuple(item for _module_value, extra in parsed_modules for item in extra)
        + tuple(item for _built, extra in parsed_implementations for item in extra)
    )
    architecture = Architecture(
        _req_str(raw, "id", path),
        _req_str(raw, "name", path),
        _req_str(raw, "backbone", path),
        _enum(raw, "replication", Replication, path),
        tuple(module for module, _extra in parsed_modules),
        parent=_parent(raw, path),
        citations=tuple(citations),
        edit=_edit(raw.get("derivation"), path),
    )
    return (
        architecture,
        _opt_str(raw, "realization", path),
        _opt_str(raw, "model_source", path),
        tuple(built for built, _extra in parsed_implementations),
        unmodeled,
    )


def _observation(value: object, path: Path) -> tuple[Observation, tuple[Unmodeled, ...]]:
    raw, extra = _split(value, path, _OBSERVATION_KEYS)
    built = Observation(
        _req_str(raw, "id", path),
        _req_str(raw, "implementation_id", path),
        _enum(raw, "phenotype", Phenotype, path),
        _enum(raw, "host", Host, path),
        _enum(raw, "assay", Assay, path),
        _enum(raw, "outcome", Outcome, path),
        _req_str(raw, "citation", path),
        _opt_str(raw, "note", path) or "",
    )
    return built, extra


def load_design(directory: Path) -> Design:
    """Load one design directory."""
    design_path = directory / "design.json"
    observations_path = directory / "observations.json"
    if not design_path.is_file():
        raise InputError(f"{directory}: missing design.json.")
    architecture, realization, model_source, implementations, unmodeled = _architecture(design_path)
    known = {item.id for item in implementations}
    if not observations_path.is_file():
        observations: tuple[Observation, ...] = ()
    else:
        loaded = _load_json(observations_path)
        if not isinstance(loaded, list):
            raise InputError(f"{observations_path}: expected a list.")
        parsed = [_observation(item, observations_path) for item in loaded]
        observations = tuple(observation for observation, _extra in parsed)
        unmodeled = unmodeled + tuple(
            item for _observation_value, extra in parsed for item in extra
        )
    for observation in observations:
        if observation.implementation_id not in known:
            raise InputError(
                f"{observations_path}: {observation.id} cites unknown implementation "
                f"{observation.implementation_id}."
            )
    return Design(architecture, observations, implementations, realization, model_source, unmodeled)


def load_library(root: Path | None = None) -> list[Design]:
    """Load every design directory under ``root``, in directory-name order."""
    directory = LIBRARY_ROOT if root is None else root
    if not directory.is_dir():
        raise InputError(f"{directory}: design library not found.")
    directories = sorted(path for path in directory.iterdir() if path.is_dir())
    records = [load_design(path) for path in directories]
    identifiers = [record.architecture.id for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise InputError(f"{directory}: duplicate design ids.")
    return records
