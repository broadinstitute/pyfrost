import pytest  # noqa

from pyfrost import Kmer, path_sequence, path_nucleotide_length


def test_build_sequence(mccortex):
    g = mccortex

    path = [Kmer('ACTGA'), Kmer('TCGAT'), Kmer('CGATG')]
    assert path_sequence(g, path) == "ACTGATTTCGATGC"

    assert path_sequence(g, []) == ""


def test_compute_length(mccortex):
    g = mccortex

    path = [Kmer('ACTGA'), Kmer('TCGAT'), Kmer('CGATG')]
    assert path_nucleotide_length(g, path) == 14
