import pytest

import pyfrost


@pytest.fixture(scope="module")
def mccortex():
    return pyfrost.build_from_refs(['data/mccortex.fasta'], k=5, g=3)


@pytest.fixture(scope="module")
def mccortex2():
    return pyfrost.build_from_refs(['data/mccortex.fasta', 'data/mccortex2.fasta'],
                                   k=5, g=3)
