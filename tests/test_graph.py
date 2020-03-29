import pytest

import pyfrost


@pytest.fixture
def mccortex():
    return pyfrost.build_from_refs(['data/mccortex.fasta'], k=5, g=3)


def test_node_access(mccortex):
    g = mccortex

    assert len(g) == 6
    assert len(g.nodes) == 6

    n = g.nodes['ACTGA']
    assert str(n) == "ACTGATTTCGA"
    assert n.is_full_mapping

    # Same node as above, but then reverse complement
    n2 = g.nodes['TCGAA']
    assert str(n2) == "TCGAAATCAGT"
    assert n2.is_full_mapping

    with pytest.raises(IndexError):
        _ = g.nodes['GGGGG']


def test_graph_attr(mccortex):
    g = mccortex

    g.graph['test'] = 1
    assert g.graph['test'] == 1


def test_successors(mccortex):
    g = mccortex

    n = g.nodes['ACTGA']
    succ = [s for s in g.succ[n]]
    assert len(succ) == 2

    succ_set = set(str(s) for s in g.succ[n])
    assert succ_set == {"TCGAAATCAGT", "TCGAT"}

    n2 = g.nodes['TCGAA']
    assert n2 in g.succ[n]
    assert n not in g.succ[n]

    succ2 = [s for s in g.succ[n2]]
    assert len(succ2) == 0
    assert set(str(s) for s in g.succ[n2]) == set()

    assert set(str(s) for s in g[n]) == succ_set
    assert set(str(s) for s in g.neighbors(n)) == succ_set

    with pytest.raises(IndexError):
        _ = g.succ['GGGGG']

    with pytest.raises(IndexError):
        _ = g['GGGGG']

    with pytest.raises(IndexError):
        _ = g.neighbors('GGGGG')

    with pytest.raises(IndexError):
        _ = g.succ[pyfrost.Kmer('GGGGG')]


def test_predecessors(mccortex):
    g = mccortex

    n = g.nodes['ACTGA']
    pred = [p for p in g.pred[n]]
    assert len(pred) == 0
    assert len(g.pred[n]) == 0

    assert set(str(p) for p in g.pred[n]) == set()

    n2 = g.nodes['TCGAA']
    print([repr(p) for p in g.pred[n2]])
    print(repr(n))

    assert n in g.pred[n2]
    assert n2 not in g.pred[n2]

    pred2 = [s for s in g.pred[n2]]
    assert len(pred2) == 2
    assert len(g.pred[n2]) == 2

    pred_set = set(str(p) for p in g.pred[n2])
    assert pred_set == {"ACTGATTTCGA", "ATCGA"}

    assert set(str(s) for s in g.predecessors(n2)) == pred_set

    with pytest.raises(IndexError):
        _ = g.pred['GGGGG']

    with pytest.raises(IndexError):
        _ = g.predecessors('GGGGG')

    with pytest.raises(IndexError):
        _ = g.pred[pyfrost.Kmer('GGGGG')]


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
