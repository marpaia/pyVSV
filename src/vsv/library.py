"""Load the curated VSV design library.

Each design is a directory of JSON. The files describe architecture, citations, and qualitative
observations. They do not carry nucleotide or amino-acid sequences.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path

from vsv.architecture import (
    Architecture,
    Derivation,
    LocusStatus,
    Module,
    ModuleKind,
    ParentStatus,
    Replication,
)
from vsv.evidence import Assay, Implementation, Observation, Outcome, Phenotype, Stage, Unmodeled
from vsv.query import DesignRecord
from vsv.rules import Host
from vsv.sequence import InputError

LIBRARY_ROOT = Path(__file__).resolve().parents[2] / "data" / "designs"

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


def _module(value: object, path: Path) -> tuple[Module, tuple[Unmodeled, ...]]:
    raw, extra = _split(value, path, _MODULE_KEYS)
    try:
        kind = ModuleKind(str(raw["kind"]))
    except (KeyError, ValueError) as error:
        raise InputError(f"{path}: a module needs a known kind.") from error
    locus = raw.get("locus")
    if locus is not None and not isinstance(locus, int):
        raise InputError(f"{path}: locus must be an integer.")
    status = raw.get("locus_status")
    antigen_class = raw.get("antigen_class")
    try:
        locus_status = None if status is None else LocusStatus(str(status))
    except ValueError as error:
        raise InputError(f"{path}: unknown locus status.") from error
    return Module(
        str(raw["id"]),
        kind,
        str(raw["name"]),
        None if locus is None else locus,
        None if antigen_class is None else str(antigen_class),
        locus_status=locus_status,
    ), extra


def _derivation(value: object, path: Path) -> Derivation | None:
    if value is None:
        return None
    raw, _extra = _split(value, path, _DERIVATION_KEYS)
    replication = raw.get("replication")
    module = raw.get("module")
    try:
        return Derivation(
            str(raw["activity"]),
            str(raw["target"]),
            None if module is None else str(module),
            None if replication is None else Replication(str(replication)),
        )
    except (KeyError, ValueError) as error:
        raise InputError(f"{path}: invalid derivation.") from error


def _implementation(value: object, path: Path) -> tuple[Implementation, tuple[Unmodeled, ...]]:
    raw, extra = _split(value, path, _IMPLEMENTATION_KEYS)
    try:
        built = Implementation(
            str(raw["id"]),
            Stage(str(raw["stage"])),
            str(raw["citation"]),
            "" if raw.get("note") is None else str(raw["note"]),
        )
    except (KeyError, ValueError) as error:
        raise InputError(f"{path}: invalid implementation.") from error
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
    parent = raw.get("parent_id")
    parent_status = raw.get("parent_status")
    try:
        status = ParentStatus.ROOT if parent_status is None else ParentStatus(str(parent_status))
    except ValueError as error:
        raise InputError(f"{path}: unknown parent status.") from error
    realization = raw.get("realization")
    model_source = raw.get("model_source")
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
    try:
        architecture = Architecture(
            str(raw["id"]),
            str(raw["name"]),
            str(raw["backbone"]),
            Replication(str(raw["replication"])),
            tuple(module for module, _extra in parsed_modules),
            parent_id=None if parent is None else str(parent),
            parent_status=status,
            citations=tuple(str(item) for item in citations),
            derivation=_derivation(raw.get("derivation"), path),
        )
    except (KeyError, ValueError) as error:
        raise InputError(f"{path}: invalid design.") from error
    return (
        architecture,
        None if realization is None else str(realization),
        None if model_source is None else str(model_source),
        tuple(built for built, _extra in parsed_implementations),
        unmodeled,
    )


def _observation(
    value: object, path: Path, architecture_id: str
) -> tuple[Observation, tuple[Unmodeled, ...]]:
    raw, extra = _split(value, path, _OBSERVATION_KEYS)
    try:
        built = Observation(
            str(raw["id"]),
            architecture_id,
            str(raw["implementation_id"]),
            Phenotype(str(raw["phenotype"])),
            Host(str(raw["host"])),
            Assay(str(raw["assay"])),
            Outcome(str(raw["outcome"])),
            str(raw["citation"]),
            "" if raw.get("note") is None else str(raw["note"]),
        )
    except (KeyError, ValueError) as error:
        raise InputError(f"{path}: invalid observation.") from error
    return built, extra


def load_design(directory: Path) -> DesignRecord:
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
        parsed = [_observation(item, observations_path, architecture.id) for item in loaded]
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
    return DesignRecord(
        architecture, observations, implementations, realization, model_source, unmodeled
    )


def load_library(root: Path | None = None) -> list[DesignRecord]:
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
