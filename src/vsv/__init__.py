"""Computational genetics tools for vesicular stomatitis virus research."""

from importlib.metadata import version

import sbol3

__version__ = version("vsv")

NAMESPACE = "https://vsv.bio"

# SBOL objects created without a namespace trigger a warning, and sbol_utilities creates some at
# import time, so the namespace is set before any vsv module is imported.
sbol3.set_namespace(NAMESPACE)

from vsv.api import GeneWarning, load, load_all  # noqa: E402
from vsv.architecture import (  # noqa: E402
    Architecture,
    Derivation,
    LocusStatus,
    Module,
    ModuleKind,
    ParentStatus,
    Replication,
    add_order_constraints,
    architecture_of,
    canonical_chassis,
    project,
    to_document,
)
from vsv.evidence import ChassisCard, Implementation, Observation, Stage, Unmodeled  # noqa: E402
from vsv.genes import Gene, Product  # noqa: E402
from vsv.library import load_library  # noqa: E402
from vsv.model import parameters_for  # noqa: E402
from vsv.query import select  # noqa: E402
from vsv.rules import Host, Property, Violation, check, satisfies  # noqa: E402
from vsv.sbol import GeneAnnotation, Genome  # noqa: E402
from vsv.sequence import InputError  # noqa: E402
from vsv.transform import DeleteModule, InsertModule, ReplaceModule  # noqa: E402

__all__ = [
    "NAMESPACE",
    "Architecture",
    "ChassisCard",
    "DeleteModule",
    "Derivation",
    "Gene",
    "GeneAnnotation",
    "GeneWarning",
    "Genome",
    "Host",
    "Implementation",
    "InputError",
    "InsertModule",
    "LocusStatus",
    "Module",
    "ModuleKind",
    "Observation",
    "ParentStatus",
    "Product",
    "Property",
    "ReplaceModule",
    "Replication",
    "Stage",
    "Unmodeled",
    "Violation",
    "__version__",
    "add_order_constraints",
    "architecture_of",
    "canonical_chassis",
    "check",
    "load",
    "load_all",
    "load_library",
    "parameters_for",
    "project",
    "satisfies",
    "select",
    "to_document",
]
