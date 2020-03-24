import pytest

import pyfrost


@pytest.fixture
def mccortex():
    return pyfrost.build_from_refs(['data/mccortex.fasta'], k=5, g=3)


def test_graph_structure(mccortex):
    g = mccortex

    assert len(g) == 6
    assert len(g.nodes) == 6

    n = g.nodes['ACTGA']
    assert str(n) == "ACTGATTTCGA"

    with pytest.raises(IndexError):
        _ = g.nodes['GGGGG']

    succ = [s for s in g.succ[n]]
    assert len(succ) == 2
    assert set(str(s) for s in g.succ[n]) == {"TCGAAATCAGT", "TCGAT"}

    n2 = g.nodes['TCGAA']
    assert str(n2) == "TCGAAATCAGT"

    assert n2 in g.succ[n]
    assert n not in g.succ[n]

    succ2 = [s for s in g.succ[n2]]
    assert len(succ2) == 0
    assert set(str(s) for s in g.succ[n2]) == set()


def test_iteration(mccortex):
    g = mccortex

    node_set = set()
    for n in g:
        node_set.add(str(n))
        node_set.add(pyfrost.reverse_complement(str(n)))

    truth = {
        "ACTGATTTCGA", "TCGAAATCAGT",
        "ATCGA", "TCGAT",
        "CGATGC", "GCATCG",
        "CCACCGTGG", "CCACGGTGG",
        "ATCGCAT", "ATGCGAT",
        "ATGCCAC", "GTGGCAT"
    }

    assert node_set == truth

    node_set = set()
    for n in g.nodes:
        node_set.add(str(n))
        node_set.add(pyfrost.reverse_complement(str(n)))

    assert node_set == truth
