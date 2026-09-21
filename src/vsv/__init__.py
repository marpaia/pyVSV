"""Computational genetics tools for vesicular stomatitis virus research."""

from importlib.metadata import version

import sbol3

__version__ = version("vsv")

NAMESPACE = "https://vsv.bio"

# SBOL objects created without a namespace trigger a warning, and sbol_utilities creates some at
# import time, so the namespace is set before any vsv module is imported.
sbol3.set_namespace(NAMESPACE)

from vsv.api import GeneWarning, load, load_all  # noqa: E402
from vsv.genes import Gene, Product  # noqa: E402
from vsv.sbol import GeneAnnotation, Genome  # noqa: E402
from vsv.sequence import InputError  # noqa: E402

__all__ = [
    "NAMESPACE",
    "Gene",
    "GeneAnnotation",
    "GeneWarning",
    "Genome",
    "InputError",
    "Product",
    "__version__",
    "load",
    "load_all",
]
