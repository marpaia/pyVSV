"""An antigen class marked on a coding feature.

The class is a role on the feature in the genome's SBOL document. Known classes are the ones the
ontology already names. A class for a virus the ontology does not name yet is
``AntigenClass("paramyxovirus-glycoprotein")``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import ClassVar

import sbol3

from vsv import NAMESPACE
from vsv.errors import InputError
from vsv.genome.sequence import Nucleotides

_SLUG = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*\Z")
_ANTIGEN_CLASS = f"{NAMESPACE}/antigen-class"

# (attribute, slug). The attribute names are declared on AntigenClass below.
_KNOWN: tuple[tuple[str, str], ...] = (
    ("filovirus_glycoprotein", "filovirus-glycoprotein"),
    ("sudan_ebolavirus_glycoprotein", "sudan-ebolavirus-glycoprotein"),
    ("bundibugyo_ebolavirus_glycoprotein", "bundibugyo-ebolavirus-glycoprotein"),
    ("reston_ebolavirus_glycoprotein", "reston-ebolavirus-glycoprotein"),
    ("taiforest_ebolavirus_glycoprotein", "taiforest-ebolavirus-glycoprotein"),
    ("marburgvirus_glycoprotein", "marburgvirus-glycoprotein"),
    ("arenavirus_gpc", "arenavirus-gpc"),
    ("henipavirus_fusion", "henipavirus-fusion"),
    ("henipavirus_glycoprotein", "henipavirus-glycoprotein"),
    ("sarbecovirus_spike", "sarbecovirus-spike"),
    ("influenza_hemagglutinin", "influenza-hemagglutinin"),
    ("hiv_envelope", "hiv-envelope"),
    ("pneumovirus_fusion", "pneumovirus-fusion"),
    ("nairovirus_glycoprotein", "nairovirus-glycoprotein"),
    ("zika_envelope", "zika-envelope"),
    ("hantavirus_glycoprotein", "hantavirus-glycoprotein"),
    ("chikungunya_envelope", "chikungunya-envelope"),
    ("rabies_glycoprotein", "rabies-glycoprotein"),
    ("bandavirus_glycoprotein", "bandavirus-glycoprotein"),
    ("herpesvirus_glycoprotein", "herpesvirus-glycoprotein"),
    ("papillomavirus_capsid", "papillomavirus-capsid"),
    ("dengue_envelope", "dengue-envelope"),
    ("gm_csf", "gm-csf"),
    ("interferon_beta", "interferon-beta"),
    ("interleukin_4", "interleukin-4"),
    ("thymidine_kinase", "thymidine-kinase"),
    ("hepacivirus_glycoprotein", "hepacivirus-glycoprotein"),
    ("hepacivirus_structural", "hepacivirus-structural"),
    ("west_nile_envelope", "west-nile-envelope"),
    ("yersinia_lcrv", "yersinia-lcrv"),
    ("vesiculovirus_glycoprotein", "vesiculovirus-glycoprotein"),
    ("henipavirus_nucleoprotein", "henipavirus-nucleoprotein"),
)


@dataclass(frozen=True, slots=True)
class AntigenClass:
    """The antigen class of a payload, stored as ``https://vsv.bio/antigen-class/{slug}``."""

    slug: str
    filovirus_glycoprotein: ClassVar[AntigenClass]
    sudan_ebolavirus_glycoprotein: ClassVar[AntigenClass]
    bundibugyo_ebolavirus_glycoprotein: ClassVar[AntigenClass]
    reston_ebolavirus_glycoprotein: ClassVar[AntigenClass]
    taiforest_ebolavirus_glycoprotein: ClassVar[AntigenClass]
    marburgvirus_glycoprotein: ClassVar[AntigenClass]
    arenavirus_gpc: ClassVar[AntigenClass]
    henipavirus_fusion: ClassVar[AntigenClass]
    henipavirus_glycoprotein: ClassVar[AntigenClass]
    sarbecovirus_spike: ClassVar[AntigenClass]
    influenza_hemagglutinin: ClassVar[AntigenClass]
    hiv_envelope: ClassVar[AntigenClass]
    pneumovirus_fusion: ClassVar[AntigenClass]
    nairovirus_glycoprotein: ClassVar[AntigenClass]
    zika_envelope: ClassVar[AntigenClass]
    hantavirus_glycoprotein: ClassVar[AntigenClass]
    chikungunya_envelope: ClassVar[AntigenClass]
    rabies_glycoprotein: ClassVar[AntigenClass]
    bandavirus_glycoprotein: ClassVar[AntigenClass]
    herpesvirus_glycoprotein: ClassVar[AntigenClass]
    papillomavirus_capsid: ClassVar[AntigenClass]
    dengue_envelope: ClassVar[AntigenClass]
    gm_csf: ClassVar[AntigenClass]
    interferon_beta: ClassVar[AntigenClass]
    interleukin_4: ClassVar[AntigenClass]
    thymidine_kinase: ClassVar[AntigenClass]
    hepacivirus_glycoprotein: ClassVar[AntigenClass]
    hepacivirus_structural: ClassVar[AntigenClass]
    west_nile_envelope: ClassVar[AntigenClass]
    yersinia_lcrv: ClassVar[AntigenClass]
    vesiculovirus_glycoprotein: ClassVar[AntigenClass]
    henipavirus_nucleoprotein: ClassVar[AntigenClass]

    def __post_init__(self) -> None:
        if _SLUG.fullmatch(self.slug) is None:
            raise InputError(f"Invalid antigen class {self.slug!r}.")

    @property
    def role(self) -> str:
        """The SBOL role URI for this class."""
        return f"{_ANTIGEN_CLASS}/{self.slug}"


for _name, _slug in _KNOWN:
    setattr(AntigenClass, _name, AntigenClass(_slug))


@dataclass(frozen=True, slots=True)
class AntigenFeature:
    """A genome feature that carries one antigen-class role and its coding sequence."""

    feature: sbol3.SequenceFeature | sbol3.SubComponent
    antigen_class: AntigenClass
    sequence: Nucleotides
