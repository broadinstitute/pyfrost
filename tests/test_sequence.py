import pytest  # noqa

from pyfrost import Kmer, sequence_for_path


def test_build_sequence(mccortex):
    g = mccortex

    path = [Kmer('ACTGA'), Kmer('TCGAT'), Kmer('CGATG')]
    assert sequence_for_path(g, path) == "ACTGATTTCGATGC"

    assert sequence_for_path(g, []) == ""
