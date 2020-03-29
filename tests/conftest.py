import pytest

import pyfrost


@pytest.fixture(scope="module")
def mccortex():
    return pyfrost.build_from_refs(['data/mccortex.fasta'], k=5, g=3)
