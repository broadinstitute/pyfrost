import pytest

import bifrost_python


@pytest.fixture
def mccortex():
    return bifrost_python.build_from_refs(['data/mccortex.fasta'], k=5, g=3)


def test_graph_structure(mccortex):
    g = mccortex

    assert len(g) == 6
    assert len(g.nodes) == 6

    n = g.nodes['ACTGA']
    assert str(n) == "ACTGATTTCGA"

    succ = [s for s in g.succ[n]]
    assert len(succ) == 2
    assert set(str(s) for s in g.succ[n]) == {"TCGAAATCAGT", "TCGAT"}

    n2 = g.nodes['TCGAA']
    assert str(n2) == "TCGAAATCAGT"

    succ2 = [s for s in g.succ[n2]]
    assert len(succ2) == 0
    assert set(str(s) for s in g.succ[n2]) == set()

