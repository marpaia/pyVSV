"""Computational genetics tools for vesicular stomatitis virus research."""

from importlib.metadata import version

import sbol3

__version__ = version("vsv")

NAMESPACE = "https://vsv.bio"

# SBOL objects created without a namespace trigger a warning, and sbol_utilities creates some at
# import time, so the namespace is set before any vsv module is imported.
sbol3.set_namespace(NAMESPACE)

from vsv.design import (  # noqa: E402
    Antigen,
    Architecture,
    ChassisCard,
    Delete,
    Design,
    Host,
    Implementation,
    Insert,
    Leader,
    LocusStatus,
    Module,
    ModuleKind,
    Observation,
    ParentStatus,
    Property,
    Replace,
    Replication,
    Reporter,
    Stage,
    Trailer,
    Unmodeled,
    Violation,
    ViralGene,
    add_order_constraints,
    architecture_of,
    canonical_chassis,
    check,
    design_of,
    design_to_document,
    load_library,
    parameters_for,
    project,
    satisfies,
    select,
    to_document,
)
from vsv.errors import InputError  # noqa: E402
from vsv.genes import Gene, Product  # noqa: E402
from vsv.genome import GeneAnnotation, GeneWarning, Genome, load, load_all  # noqa: E402

__all__ = [
    "NAMESPACE",
    "Antigen",
    "Architecture",
    "ChassisCard",
    "Delete",
    "Design",
    "Gene",
    "GeneAnnotation",
    "GeneWarning",
    "Genome",
    "Host",
    "Implementation",
    "InputError",
    "Insert",
    "Leader",
    "LocusStatus",
    "Module",
    "ModuleKind",
    "Observation",
    "ParentStatus",
    "Product",
    "Property",
    "Replace",
    "Replication",
    "Reporter",
    "Stage",
    "Trailer",
    "Unmodeled",
    "Violation",
    "ViralGene",
    "__version__",
    "add_order_constraints",
    "architecture_of",
    "canonical_chassis",
    "check",
    "design_of",
    "design_to_document",
    "load",
    "load_all",
    "load_library",
    "parameters_for",
    "project",
    "satisfies",
    "select",
    "to_document",
]
