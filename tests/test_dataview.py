import pytest  # noqa


def test_node_dataview(mccortex):
    g = mccortex

    # Test call without arguments, should return the same as without call
    node_set = set(n for n in g.nodes)
    node_set2 = set(n for n in g.nodes())

    assert node_set == node_set2

    # Now with data dict
    node_set3 = set()
    for n, data in g.nodes(data=True):
        assert isinstance(data, dict)

        node_set3.add(n)

    assert node_set == node_set3

    # Now with specific key in the data dict
    for n, data in g.nodes(data="colors"):
        # Right now the data dict is still empty
        assert data is None

    # Now with given default value if key doesn't exists
    for n, data in g.nodes(data="non_existent", default=False):
        assert data is False
