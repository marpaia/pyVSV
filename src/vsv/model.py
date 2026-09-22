"""Parameters a VSV architecture can supply to an intracellular replication model.

The linked model is the stochastic intracellular model of King, Berezin, and Peccoud (2024).
This module records that citation and the architectural parameters the model can consume. It
does not simulate the infection.
"""

from __future__ import annotations

import sbol3

from vsv.architecture import Architecture, ModuleKind

INTRACELLULAR_MODEL_SOURCE = "https://doi.org/10.1371/journal.pcbi.1011373"
INTRACELLULAR_MODEL_LANGUAGE = "https://vsv.bio/language/stochastic-simulation"
INTRACELLULAR_MODEL_FRAMEWORK = "https://identifiers.org/SBO:0000063"


def parameters_for(architecture: Architecture) -> dict[str, object]:
    """Gene count, ordered locus names, and replication class for ``architecture``."""
    transcribed = [
        module
        for module in architecture.modules
        if module.kind not in {ModuleKind.LEADER, ModuleKind.TRAILER}
    ]
    return {
        "replication": architecture.replication.value,
        "gene_count": len(transcribed),
        "loci": [module.id for module in transcribed],
    }


def attach_model(document: sbol3.Document, architecture_id: str) -> sbol3.Model:
    """Link the King 2024 intracellular model to the architecture ``architecture_id``."""
    component = document.find(f"https://vsv.bio/{architecture_id}")
    model = sbol3.Model(
        f"{architecture_id}_intracellular",
        INTRACELLULAR_MODEL_SOURCE,
        INTRACELLULAR_MODEL_LANGUAGE,
        INTRACELLULAR_MODEL_FRAMEWORK,
        name="Stochastic intracellular VSV replication",
    )
    document.add(model)
    component.models = [model]
    return model
