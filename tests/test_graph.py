import pytest

import pyfrost
from pyfrost import Kmer


def test_find_unitig(mccortex):
    g = mccortex

    assert len(g) == 12
    assert len(g.nodes) == 12

    # Single k-mer which is a head of a unitig
    kmer1 = Kmer('ACTGA')
    assert str(kmer1) == "ACTGA"

    # Kmer1 is the head of a unitig in the graph
    assert kmer1 in g
    assert kmer1 in g.nodes

    # This kmer is in the middle of a unitig, so should not be a node
    kmer2 = Kmer('GATTT')
    assert kmer2 not in g
    assert kmer2 not in g.nodes

    # Should also works with strings
    assert "GATTT" not in g
    assert "GATTT" not in g.nodes

    # Same node as kmer1, but then reverse complement
    kmer3 = Kmer('TCGAA')
    assert kmer3 in g
    assert kmer3 in g.nodes

    assert str(kmer3) == "TCGAA"


def test_graph_attr(mccortex):
    g = mccortex

    assert g.graph['k'] == 5

    g.graph['test'] = 1
    assert g.graph['test'] == 1


def test_successors(mccortex):
    g = mccortex

    n = Kmer('ACTGA')
    succ = set(g.succ[n])
    assert len(succ) == 2
    assert succ == {Kmer("TCGAA"), Kmer("TCGAT")}

    n2 = Kmer('TCGAA')
    assert n2 in g.succ[n]
    assert n not in g.succ[n]

    succ2 = set(g.succ[n2])
    assert len(succ2) == 0
    assert succ2 == set()

    # Different ways to access successors
    assert set(g[n]) == succ
    assert set(g.neighbors(n)) == succ

    with pytest.raises(IndexError):
        _ = g.succ['GGGGG']

    with pytest.raises(IndexError):
        _ = g['GGGGG']

    with pytest.raises(IndexError):
        _ = g.neighbors('GGGGG')

    with pytest.raises(IndexError):
        _ = g.succ[Kmer('GGGGG')]


def test_predecessors(mccortex):
    g = mccortex

    n = Kmer('ACTGA')
    pred = set(g.pred[n])
    assert len(pred) == 0
    assert len(g.pred[n]) == 0
    assert pred == set()

    n2 = Kmer('TCGAA')

    assert n in g.pred[n2]
    assert n2 not in g.pred[n2]

    pred2 = set(g.pred[n2])
    assert len(pred2) == 2
    assert len(g.pred[n2]) == 2

    assert pred2 == {Kmer("ACTGA"), Kmer("ATCGA")}
    assert set(g.predecessors(n2)) == pred2

    with pytest.raises(IndexError):
        _ = g.pred['GGGGG']

    with pytest.raises(IndexError):
        _ = g.predecessors('GGGGG')

    with pytest.raises(IndexError):
        _ = g.pred[pyfrost.Kmer('GGGGG')]


def test_iteration(mccortex):
    g = mccortex

    nodes = set(g)
    truth = {
        Kmer("ACTGA"), Kmer("TCGAA"),
        Kmer("ATCGA"), Kmer("TCGAT"),
        Kmer("CGATG"), Kmer("GCATC"),
        Kmer("CCACC"), Kmer("CCACG"),
        Kmer("ATCGC"), Kmer("ATGCG"),
        Kmer("ATGCC"), Kmer("GTGGC")
    }

    assert nodes == truth

    n1 = Kmer('ACTGA')
    n2 = Kmer('TCGAA')
    kmer = Kmer('ATTTC')

    # Ignores any value in the list that's not a node
    assert set(g.nbunch_iter(["ACTGA", n2, kmer])) == {
        Kmer("ACTGA"), Kmer("TCGAA")}
    assert set(g.nbunch_iter([n1, n2, 20])) == {n1, n2}

    # neighbors(), successors() and predecessors() all should return an iterator
    assert next(g.neighbors(n1))
    assert next(g.successors(n1))

    assert next(g.predecessors(n2))


def test_in_out_edges(mccortex):
    g = mccortex

    for s, t in g.edges:
        successors = set(g.successors(s))
        assert t in successors

    assert len(g.edges) == 16

    assert ("ACTGA", "TCGAT") in g.edges

    n1 = Kmer("ACTGA")
    n2 = Kmer("TCGAT")

    assert (n1, n2) in g.edges

    n3 = Kmer("ATCGA")

    assert (n1, n3) not in g.edges

    assert g.edges[n1, n2]['orientation'] == (g.nodes[n1]['strand'], g.nodes[n2]['strand'])
    assert g.edges[n1, n2]['label'] == "T"
