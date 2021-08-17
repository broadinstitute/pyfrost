import pytest

import pyfrost
from pyfrost import links


@pytest.fixture(scope="module")
def mccortex():
    return pyfrost.build_from_refs(['data/mccortex.fasta'], k=5, g=3)


@pytest.fixture(scope="module")
def mccortex2():
    return pyfrost.build_from_refs(['data/mccortex.fasta', 'data/mccortex2.fasta'],
                                   k=5, g=3)


@pytest.fixture(scope="module")
def linked_mccortex(mccortex):
    linkdb = links.MemLinkDB()
    links.add_links_from_sequence(mccortex, linkdb, "TTTCGATGCGATGCGATGCCACG")

    return mccortex, linkdb
