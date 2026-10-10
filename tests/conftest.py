"""Shared fixtures."""

from pathlib import Path

import pytest

from tests.synthetic import SyntheticGenome, build_genome


@pytest.fixture
def synthetic() -> SyntheticGenome:
    return build_genome()


@pytest.fixture
def synthetic_fasta(tmp_path: Path, synthetic: SyntheticGenome) -> Path:
    path = tmp_path / "synthetic.fasta"
    path.write_text(synthetic.fasta)
    return path
