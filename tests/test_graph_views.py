import pytest  # noqa

from pyfrost import Kmer


def test_subgraph(mccortex):
    g = mccortex

    nodes = [
        Kmer('CGATG'),
        Kmer('ATGCC'),
        Kmer('CCACC')
    ]

    sg = g.subgraph(nodes)

    for n in nodes:
        assert n in sg

    assert Kmer('ACTGA') not in sg
    assert set(sg) == set(nodes)

    assert (Kmer('ATGCC'), Kmer('CCACC')) in sg.edges
    assert (Kmer('ATGCC'), Kmer('CCACG')) not in sg.edges

    print(list(sg.edges))

    assert list(sg.edges) == [
        (Kmer('CGATG'), Kmer('ATGCC')),
        (Kmer('ATGCC'), Kmer('CCACC')),
    ]



