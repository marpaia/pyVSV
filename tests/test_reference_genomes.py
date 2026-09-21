"""Check the heuristic annotator against curated GenBank records in ``data/raw``.

The data directory holds downloaded inputs and is not part of the repository, so these tests
are skipped for any accession whose FASTA and GenBank files are not both present.
"""

from pathlib import Path

import pytest

from vsv import sbol
from vsv.annotation import annotate
from vsv.formats import read_genomes
from vsv.genes import Gene

RAW = Path(__file__).resolve().parents[1] / "data" / "raw"
ACCESSIONS = ["MF196237", "MF196238", "NC_001560"]


@pytest.fixture(params=ACCESSIONS)
def accession(request: pytest.FixtureRequest) -> str:
    name = str(request.param)
    if not (RAW / f"{name}.fasta").exists() or not (RAW / f"{name}.gb").exists():
        pytest.skip(f"{name} FASTA and GenBank files are not in data/raw")
    return name


def test_annotated_coding_sequences_match_the_curated_record(accession: str) -> None:
    _, (heuristic,) = read_genomes([RAW / f"{accession}.fasta"])
    _, (curated,) = read_genomes([RAW / f"{accession}.gb"])
    assert heuristic.bases == curated.bases
    annotate(heuristic)
    for name in Gene:
        found, expected = heuristic.cds(name), curated.cds(name)
        assert found is not None, f"{name} was not found"
        assert expected is not None, f"{name} is missing from the curated record"
        assert heuristic.span(found) == curated.span(expected), name
        assert heuristic.protein(name) == curated.protein(name), name


def test_curated_genes_enclose_the_annotated_coding_sequences(accession: str) -> None:
    _, (curated,) = read_genomes([RAW / f"{accession}.gb"])
    genes = {str(f.name): curated.span(f) for f in curated.features(sbol.SO_GENE)}
    for name, (start, end) in genes.items():
        cds_start, cds_end = curated.span(curated.cds(name))
        assert start <= cds_start < cds_end <= end, name


def test_annotated_genes_match_the_curated_genes() -> None:
    if not (RAW / "NC_001560.fasta").exists() or not (RAW / "NC_001560.gb").exists():
        pytest.skip("NC_001560 FASTA and GenBank files are not in data/raw")
    _, (heuristic,) = read_genomes([RAW / "NC_001560.fasta"])
    _, (curated,) = read_genomes([RAW / "NC_001560.gb"])
    annotate(heuristic)
    for name in Gene:
        found, expected = heuristic.gene(name), curated.gene(name)
        assert found is not None and expected is not None, name
        assert heuristic.span(found) == curated.span(expected), name
