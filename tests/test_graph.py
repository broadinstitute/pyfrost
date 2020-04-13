import pytest

import pyfrost


def test_find_unitig(mccortex):
    g = mccortex

    assert len(g) == 6
    assert len(g.nodes) == 6

    # Single k-mer which is a head of a unitig
    kmer1 = g.find('ACTGA', True)

    # A single k-mer should not be a node of the graph
    assert kmer1 not in g
    assert kmer1 not in g.nodes

    # Get the full unitig
    n1 = kmer1.get_full_mapping()

    assert str(kmer1) == "ACTGA"
    assert str(n1) == "ACTGATTTCGA"

    # Same node as above, but then reverse complement
    kmer2 = g.find('TCGAA', True)
    n2 = kmer2.get_full_mapping()
    assert str(kmer2) == "TCGAA"
    assert str(n2) == "TCGAAATCAGT"

    # Single k-mer in the middle of an unitig
    kmer3 = g.find('ATTTC')
    n3 = kmer3.get_full_mapping()

    assert str(kmer3) == "ATTTC"
    assert str(n3) == "ACTGATTTCGA"

    # If extremities_only=True, then the above should return False
    assert not g.find("ATTTC", True)

    # Non-existent k-mer
    assert not g.find("GGGGG")


def test_graph_attr(mccortex):
    g = mccortex

    assert g.graph['k'] == 5
    assert isinstance(g.graph['color_names'], list)

    g.graph['test'] = 1
    assert g.graph['test'] == 1


def test_successors(mccortex):
    g = mccortex

    n = g.find('ACTGA', True).get_full_mapping()
    succ = [s for s in g.succ[n]]
    assert len(succ) == 2

    succ_set = set(str(s) for s in g.succ[n])
    assert succ_set == {"TCGAAATCAGT", "TCGAT"}

    n2 = g.find('TCGAA', True).get_full_mapping()
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

    n = g.find('ACTGA', True).get_full_mapping()
    pred = [p for p in g.pred[n]]
    assert len(pred) == 0
    assert len(g.pred[n]) == 0

    assert set(str(p) for p in g.pred[n]) == set()

    n2 = g.find('TCGAA', True).get_full_mapping()

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

    nodes = [n for n in g]
    nodes_set = set(str(n) for n in nodes)

    assert nodes_set == {
        "TCGAAATCAGT",
        "CCACGGTGG",
        "CGATGC",
        "ATGCGAT",
        "GTGGCAT",
        "ATCGA",
    }

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

    u1 = g.find('ACTGA', True).get_full_mapping()
    u2 = g.find('TCGAA', True).get_full_mapping()
    kmer1 = g.find('ATTTC')

    # Discards any UnitigMapping that's not a node (i.e. a full unitig mapping)
    assert set(g.nbunch_iter([u1, u2, kmer1])) == {u1, u2}
    assert set(g.nbunch_iter([u1, u2, 20])) == {u1, u2}

    # neighbors(), successors() and predecessors() all should return an iterator
    assert next(g.neighbors(u1))
    assert next(g.successors(u1))

    assert next(g.predecessors(u2))
