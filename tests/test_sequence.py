import pickle
import pytest  # noqa

from pyfrost import Kmer, path_sequence, path_nucleotide_length, path_kmers, kmerize_str, set_k


def test_kmer_pickle(mccortex, tmp_path):
    node_list = list(mccortex.nodes)
    with open(tmp_path / "test.pickle", "wb") as o:
        pickle.dump(node_list, o)

    with open(tmp_path / "test.pickle", "rb") as f:
        node_list2 = pickle.load(f)

    assert node_list == node_list2


def test_build_sequence(mccortex):
    g = mccortex

    path = [Kmer('ACTGA'), Kmer('TCGAT'), Kmer('CGATG')]
    assert path_sequence(g, path) == "ACTGATTTCGATGC"

    assert path_sequence(g, []) == ""


def test_compute_length(mccortex):
    g = mccortex

    path = [Kmer('ACTGA'), Kmer('TCGAT'), Kmer('CGATG')]
    assert path_nucleotide_length(g, path) == 14


def test_kmerize_str():
    set_k(5)

    kmers = list(kmerize_str("ACTGATTTCGATGC"))
    assert kmers == [
        Kmer('ACTGA'),
        Kmer('CTGAT'),
        Kmer('TGATT'),
        Kmer('GATTT'),
        Kmer('ATTTC'),
        Kmer('TTTCG'),
        Kmer('TTCGA'),
        Kmer('TCGAT'),
        Kmer('CGATG'),
        Kmer('GATGC')
    ]


def test_kmerize_path(mccortex):
    set_k(5)
    g = mccortex

    kmers = list(path_kmers(g, [Kmer('ACTGA'), Kmer('TCGAT')]))
    assert kmers == [
        Kmer('ACTGA'),
        Kmer('CTGAT'),
        Kmer('TGATT'),
        Kmer('GATTT'),
        Kmer('ATTTC'),
        Kmer('TTTCG'),
        Kmer('TTCGA'),
        Kmer('TCGAT')
    ]
