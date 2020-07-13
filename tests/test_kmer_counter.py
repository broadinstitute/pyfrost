from collections import Counter

import pytest  # noqa

from pyfrost import Kmer, KmerCounter


def test_kmer_counter():
    test_str = "ACTGATTTCGATGCGATGCGATGCCACGGTGG"
    truth_counter = Counter(Kmer(test_str[i:i+5]) for i in range(len(test_str) - 5))

    counter = KmerCounter(5, 3).count_kmers(test_str)

    for kmer, truth in truth_counter.items():
        assert counter.query(kmer) == truth

    counter2 = KmerCounter(5, 3).count_kmers_files(['data/mccortex.fasta'])

    for kmer, truth in truth_counter.items():
        assert counter2.query(kmer) == truth

    truth_counter_rep = Counter(Kmer(test_str[i:i+5]).rep() for i in range(len(test_str) - 5))

    counter3 = KmerCounter(5, 3).count_kmers(test_str)

    for kmer, truth in truth_counter_rep.items():
        assert counter3.query(kmer) == truth

    counter4 = KmerCounter(5, 3).count_kmers_files(['data/mccortex.fasta'])

    for kmer, truth in truth_counter_rep.items():
        assert counter4.query(kmer) == truth
