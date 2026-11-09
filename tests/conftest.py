"""Shared paths for the test suite."""

from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"


@pytest.fixture
def chr19_fancy():
    """HCT116 NPM1-antibody nucleolar FANCY sort, replicate 1, chr19 only."""
    return str(DATA / "HCT116_NPM1_Ab_sort_rep1.chr19.w50k_b1k.bedgraph.gz")


@pytest.fixture
def chr19_input():
    """Matched input of the same replicate, chr19 only."""
    return str(DATA / "HCT116_NPM1_Ab_sort_rep1_input.chr19.w50k_b1k.bedgraph.gz")


@pytest.fixture
def expected():
    return DATA / "expected"
