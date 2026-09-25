"""A plasmid for one vaccine plan.

The plan file has no coding sequence. The plasmid file carries the marked coding sequence on one
insert part and names the bought backbone for the plan's slot.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import sbol3
from lab.protocols import Part

from vsv import NAMESPACE
from vsv.errors import InputError
from vsv.genome.sbol import SO_CDS, new_document
from vsv.genome.sequence import Nucleotides

_BACKBONE = f"{NAMESPACE}/backbone/"
_KEYS = ("plan", "edit", "slot", "insert")


@dataclass(frozen=True, slots=True)
class Plasmid:
    """The insert, backbone, and product for one plan."""

    plan_id: str
    edit: str
    slot: str
    insert_id: str
    insert: Part
    sequence: Nucleotides
    backbone: Part
    product: Part

    def write(self, path: str | Path) -> None:
        """Write this plasmid as an SBOL3 file. The suffix selects the format."""
        document = _document(self)
        document.write(path)


def read_plasmid(path: str | Path) -> Plasmid:
    """Read the plasmid in an SBOL3 file."""
    document = sbol3.Document()
    try:
        document.read(str(path))
    except (OSError, ValueError) as error:
        raise InputError(f"{path}: {error}") from error
    products = [
        obj
        for obj in document.objects
        if isinstance(obj, sbol3.Component) and str(obj.identity).endswith("/plasmid")
    ]
    if len(products) != 1:
        raise InputError(f"{path}: expected one plasmid, found {len(products)}.")
    product = products[0]
    facts = _facts(str(product.description or ""))
    missing = [key for key in _KEYS if key not in facts]
    if missing:
        raise InputError(f"{path}: missing {', '.join(missing)}.")
    backbones = [str(role) for role in product.roles if str(role).startswith(_BACKBONE)]
    if len(backbones) != 1:
        raise InputError(f"{path}: expected one backbone, found {len(backbones)}.")
    located = [feature for feature in product.features if isinstance(feature, sbol3.SubComponent)]
    if len(located) != 1:
        raise InputError(f"{path}: expected one insert, found {len(located)}.")
    insert = document.find(str(located[0].instance_of))
    if not isinstance(insert, sbol3.Component) or len(insert.sequences) != 1:
        raise InputError(f"{path}: the insert has no coding sequence.")
    sequence = document.find(str(insert.sequences[0]))
    if not isinstance(sequence, sbol3.Sequence) or not sequence.elements:
        raise InputError(f"{path}: the insert has no coding sequence.")
    return Plasmid(
        facts["plan"],
        facts["edit"],
        facts["slot"],
        facts["insert"],
        Part(str(insert.identity)),
        Nucleotides(str(sequence.elements).upper()),
        Part(backbones[0]),
        Part(str(product.identity)),
    )


def _document(plasmid: Plasmid) -> sbol3.Document:
    document = new_document()
    sequence = sbol3.Sequence(
        f"{plasmid.insert.iri}_seq",
        elements=plasmid.sequence.bases.lower(),
        encoding=sbol3.IUPAC_DNA_ENCODING,
    )
    insert = sbol3.Component(
        plasmid.insert.iri,
        [sbol3.SBO_DNA],
        roles=[SO_CDS],
        sequences=[sequence],
        name=plasmid.insert_id,
    )
    product = sbol3.Component(
        plasmid.product.iri,
        [sbol3.SBO_DNA],
        roles=[plasmid.backbone.iri],
        name=plasmid.plan_id,
        description=(
            f"plan={plasmid.plan_id};edit={plasmid.edit};"
            f"slot={plasmid.slot};insert={plasmid.insert_id}"
        ),
        features=[sbol3.SubComponent(insert, name=plasmid.insert_id)],
    )
    document.add(sequence)
    document.add(insert)
    document.add(product)
    return document


def _facts(description: str) -> dict[str, str]:
    facts: dict[str, str] = {}
    for part in description.split(";"):
        key, separator, value = part.partition("=")
        if separator:
            facts[key] = value
    return facts
