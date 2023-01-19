"""
Test modification of graphs (add/removing nodes etc.)
"""

from pyfrost import Kmer


def test_remove_unitig(mccortex):
    kmer = Kmer('TCGAA')

    mccortex.remove_node(kmer)

    assert kmer not in mccortex
    assert kmer.twin() not in mccortex
