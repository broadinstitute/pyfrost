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

    assert kmer1.data['is_full_mapping'] is False
    assert isinstance(kmer1.data['pos'], int)
    assert isinstance(kmer1.data['length'], int)
    assert isinstance(kmer1.data['unitig_length'], int)
    assert kmer1.data['strand'] == Strand.REVERSE

    n1 = kmer1.get_full_mapping()
    assert n1.data['is_full_mapping'] is True
    assert g.nodes[n1]['is_full_mapping'] is True
    assert kmer1.data['length'] != n1.data['length']
    assert kmer1.data['unitig_length'] == n1.data['unitig_length']
    assert n1.data['strand'] == Strand.REVERSE
