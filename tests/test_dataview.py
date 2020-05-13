from collections.abc import Set, Mapping
import pytest  # noqa

from pyfrost import Strand
from pyfrost.views import NodeView, NodeDataView


def test_node_dataview(mccortex):
    g = mccortex

    assert isinstance(g.nodes, Set)
    assert isinstance(g.nodes, Mapping)
    assert isinstance(g.nodes, NodeView)
    assert isinstance(g.nodes(data=True), NodeDataView)

    # Test call without arguments, should return the same as without call
    node_set = set(n for n in g.nodes)
    node_set2 = set(n for n in g.nodes())

    assert node_set == node_set2

    # Now with data mapping
    node_set3 = set()
    for n, data in g.nodes(data=True):
        assert isinstance(data, Mapping)
        assert data['is_full_mapping'] is True
        node_set3.add(n)

    assert node_set == node_set3

    # Check that the Set mixin works
    n = g.find('ACTGA', True).get_full_mapping()
    assert (g.nodes & {n}) == {n}

    # Now with specific key in the data dict
    for n, data in g.nodes(data="is_full_mapping"):
        # All nodes should be full mappings
        assert data is True

    # Now with given default value if key doesn't exists
    for n, data in g.nodes(data="non_existent", default=False):
        assert data is False

    # Check that the `Mapping` mixin works
    keys = set(g.nodes.keys())
    assert keys == node_set

    values = list(g.nodes.values())
    for node_data in values:
        assert node_data['is_full_mapping'] is True


def test_unitig_data(mccortex):
    g = mccortex

    kmer1 = g.find('ACTGA', True)
    kmer2 = g.find('TCGAA', True)
    kmer3 = g.find('TGATT')
    kmer4 = g.find('GAAAT')

    assert kmer1.data['is_full_mapping'] is False
    assert isinstance(kmer1.data['pos'], int)

    assert isinstance(kmer1.data['length'], int)
    assert isinstance(kmer1.data['unitig_length'], int)
    assert kmer1.data['unitig_length'] == 7

    assert kmer1.data['pos'] == 0
    assert kmer1.data['strand'] == Strand.REVERSE

    assert kmer2.data['pos'] == 0
    assert kmer2.data['strand'] == Strand.FORWARD

    assert kmer3.data['pos'] == 2
    assert kmer3.data['strand'] == Strand.REVERSE

    assert kmer4.data['pos'] == 2
    assert kmer4.data['strand'] == Strand.FORWARD

    n1 = kmer1.get_full_mapping()
    assert n1.data['is_full_mapping'] is True
    assert g.nodes[n1]['is_full_mapping'] is True
    assert kmer1.data['length'] != n1.data['length']
    assert kmer1.data['unitig_length'] == n1.data['unitig_length']
    assert n1.data['strand'] == Strand.REVERSE

    # Test user attributes
    n1.data['test'] = 1
    assert n1.data['test'] == 1

    # kmer1 points to the same unitig so should have the same value
    assert kmer1.data['test'] == 1

    with pytest.raises(KeyError):
        _ = kmer1.data['dgfdghsdfkgh']

    # Harcoded keys should be read only
    with pytest.raises(KeyError):
        kmer1.data['is_full_mapping'] = True


def test_unitig_colors(mccortex):
    g = mccortex

    kmer1 = g.find('TCGAA')
    colors = list(iter(kmer1.data['colors']))

    assert colors == [(0, 0)]

    n1 = kmer1.get_full_mapping()
    colors = list(iter(n1.data['colors']))

    assert colors == [(i, 0) for i in range(n1.data['length'])]


def test_edge_dataview(mccortex):
    g = mccortex

    for s, t in g.edges():
        successors = set(g.successors(s))
        assert t in successors

    for s, t, data in g.edges(data=True):
        assert isinstance(data, dict)
        assert 'label' in data
        assert 'orientation' in data

        orientation = data['orientation']
        assert orientation[0] == s.data['strand']
        assert orientation[1] == t.data['strand']

    for s, t, orientation in g.edges(data="orientation"):
        assert orientation[0] == s.data['strand']
        assert orientation[1] == t.data['strand']

    for s, t, non_existent in g.edges(data="dghfkjgjk"):
        assert non_existent is None

    for s, t, non_existent in g.edges(data="jfgskdfjghk", default=False):
        assert non_existent is False


def test_edge_dataview_nbunch(mccortex):
    g = mccortex

    u1 = g.find('ACTGA').get_full_mapping()
    u2 = g.find('TCGAT').get_full_mapping()

    edges = set(g.edges([u1, u2]))

    u3 = g.find('TCGAA').get_full_mapping()
    u4 = g.find('CGATG').get_full_mapping()

    assert edges == {
        (u1, u2),
        (u1, u3),
        (u2, u4)
    }

    assert len(g.edges([u1, u2])) == 3
    assert len(g.edges([u1, u2, None])) == 3

    edges = [e for e in sorted(g.edges([u1, u2], data=True),
                               key=lambda e: (e[0].head, e[1].head))]

    for n1, n2, d in edges:
        assert isinstance(d, dict)

    edges_ = [(n1, n2) for n1, n2, _ in edges]
    assert edges_ == [
        (u1, u3),
        (u1, u2),
        (u2, u4)
    ]

    edges = set(g.edges([u1, u2], data="label"))
    assert edges == {
        (u1, u2, "T"),
        (u1, u3, "A"),
        (u2, u4, "G")
    }

    edges = set(g.edges([u1, u2], data="non_existent", default="test"))
    assert edges == {
        (u1, u2, "test"),
        (u1, u3, "test"),
        (u2, u4, "test")
    }

    assert (u1, u2) in g.edges(data=True)
    assert (u1, u2, g.edges[u1, u2]) in g.edges(data=True)
    assert (u1, u2, "T") in g.edges(data="label")

    assert (u1, u4) not in g.edges(data=True)
    assert (u1, u2, "G") not in g.edges(data="label")

    assert len(g.in_edges([u3, u2])) == 4

    edges = set(g.in_edges([u3, u2]))
    u5 = g.find('ATCGA').get_full_mapping()

    assert edges == {
        (u1, u3),
        (u5, u3),
        (u1, u2),
        (u5, u2),
    }


