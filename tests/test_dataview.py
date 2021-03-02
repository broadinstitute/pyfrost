from collections.abc import Set, Mapping
import pytest  # noqa

from pyfrost import Strand, Kmer


def test_node_dataview(mccortex):
    g = mccortex

    assert isinstance(g.nodes, Set)
    assert isinstance(g.nodes, Mapping)

    # Test call without arguments, should return the same as without call
    node_set = set(g.nodes)
    node_set2 = set(g.nodes())

    assert node_set == node_set2

    # Now with data mapping
    node_set3 = set()
    for n, data in g.nodes(data=True):
        assert isinstance(data, Mapping)
        assert data['is_full_mapping'] is True
        node_set3.add(n)

    assert node_set == node_set3

    # Check that the Set mixin works
    n = Kmer('ACTGA')
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

    forward_unitigs = list(g.nodes(with_rev_compl=False))
    for kmer in forward_unitigs:
        assert g.find(kmer)['strand'] == Strand.FORWARD

    num_forward = len(list(g.nodes(with_rev_compl=False)))
    assert num_forward == int(len(g) / 2)


def test_unitig_data(mccortex):
    g = mccortex

    kmer1 = g.nodes['ACTGA']
    kmer2 = g.nodes['TCGAA']

    # These k-mers are in the middle of an unitig
    kmer3 = g.find('TGATT')
    kmer4 = g.find('GAAAT')

    assert kmer1['is_full_mapping'] is True

    assert isinstance(kmer1['length'], int)
    assert kmer1['unitig_length'] == 7
    assert str(kmer1) == "ACTGATTTCGA"

    assert kmer1['pos'] == 0
    assert kmer1['strand'] == Strand.REVERSE

    assert kmer2['pos'] == 0
    assert kmer2['strand'] == Strand.FORWARD

    assert kmer3['pos'] == 2
    assert kmer3['strand'] == Strand.REVERSE

    assert kmer4['pos'] == 2
    assert kmer4['strand'] == Strand.FORWARD

    # Test user attributes
    kmer1['test'] = 1
    assert kmer1['test'] == 1

    with pytest.raises(KeyError):
        _ = kmer1['dgfdghsdfkgh']

    # Harcoded keys should be read only
    with pytest.raises(KeyError):
        kmer1['is_full_mapping'] = True


def test_unitig_colors(mccortex2):
    g = mccortex2

    n1 = g.nodes["TCGAA"]
    colors = set(n1['colors'])

    assert colors == {0, 1}

    # Test Set mixin
    assert g.nodes['TCGAA']['colors'] & {0} == {0}

    assert set(n1['colors'][0]) == {0, 1}
    assert set(n1['colors'][-1]) == {0, 1}

    with pytest.raises(KeyError):
        _ = n1['colors'][43876]

    # Other node specific to one color
    kmer = g.find('TGGTG')
    assert set(kmer['colors']) == {1}


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
        assert orientation[0] == g.nodes[s]['strand']
        assert orientation[1] == g.nodes[t]['strand']

    for s, t, orientation in g.edges(data="orientation"):
        assert orientation[0] == g.nodes[s]['strand']
        assert orientation[1] == g.nodes[t]['strand']

    for s, t, non_existent in g.edges(data="dghfkjgjk"):
        assert non_existent is None

    for s, t, non_existent in g.edges(data="jfgskdfjghk", default=False):
        assert non_existent is False


def test_edge_dataview_nbunch(mccortex):
    g = mccortex

    kmer1 = Kmer('ACTGA')
    kmer2 = Kmer('TCGAT')

    edges = set(g.edges(['ACTGA', 'TCGAT']))

    kmer3 = Kmer('TCGAA')
    kmer4 = Kmer('CGATG')

    assert edges == {
        (kmer1, kmer2),
        (kmer1, kmer3),
        (kmer2, kmer4)
    }

    assert len(g.edges([kmer1, kmer2])) == 3
    assert len(g.edges([kmer1, kmer2, None, ()])) == 3

    edges = [e for e in sorted(g.edges([kmer1, kmer2], data=True))]

    for n1, n2, d in edges:
        assert isinstance(d, dict)

    edges_ = [(n1, n2) for n1, n2, _ in edges]
    assert edges_ == [
        (kmer1, kmer3),
        (kmer1, kmer2),
        (kmer2, kmer4)
    ]

    edges = set(g.edges([kmer1, kmer2], data="label"))
    assert edges == {
        (kmer1, kmer2, "T"),
        (kmer1, kmer3, "A"),
        (kmer2, kmer4, "G")
    }

    edges = set(g.edges([kmer1, kmer2], data="non_existent", default="test"))
    assert edges == {
        (kmer1, kmer2, "test"),
        (kmer1, kmer3, "test"),
        (kmer2, kmer4, "test")
    }

    assert (kmer1, kmer2, g.edges[kmer1, kmer2]) in g.edges(data=True)
    assert (kmer1, kmer2, "T") in g.edges(data="label")

    assert (kmer1, kmer4) not in g.edges(data=True)
    assert (kmer1, kmer2, "G") not in g.edges(data="label")

    assert len(g.in_edges([kmer3, kmer2])) == 4

    edges = set(g.in_edges([kmer3, kmer2]))
    kmer5 = Kmer('ATCGA')

    assert edges == {
        (kmer1, kmer3),
        (kmer5, kmer3),
        (kmer1, kmer2),
        (kmer5, kmer2),
    }

