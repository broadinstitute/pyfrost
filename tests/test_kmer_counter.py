from collections import Counter

import pytest  # noqa

from pyfrost import Kmer, KmerCounter


def test_kmer_counter():
    test_str = "ACTGATTTCGATGCGATGCGATGCCACGGTGG"

    counter = KmerCounter(5, 3, canonical=False).count_kmers(test_str)
    truth_counter = Counter(Kmer(test_str[i:i+5]) for i in range(len(test_str) - 5 + 1))

    for kmer, truth in truth_counter.items():
        assert counter.query(kmer) == truth

    assert set(counter.keys()) == set(truth_counter.keys())

    for kmer, count in counter.items():
        assert count == truth_counter[kmer]

    assert counter.max_count == max(truth_counter.values())

    counter2 = KmerCounter(5, 3, canonical=False).count_kmers_files(['data/mccortex.fasta'])

    for kmer, truth in truth_counter.items():
        assert counter2.query(kmer) == truth

    truth_counter_rep = Counter(Kmer(test_str[i:i+5]).rep() for i in range(len(test_str) - 5 + 1))

    counter3 = KmerCounter(5, 3).count_kmers(test_str)

    for kmer, truth in truth_counter_rep.items():
        assert counter3.query(kmer) == truth

    counter4 = KmerCounter(5, 3).count_kmers_files(['data/mccortex.fasta'])

    for kmer, truth in truth_counter_rep.items():
        assert counter4.query(kmer) == truth


def test_kmer_counter_load_save(tmp_path):
    test_str = "ACTGATTTCGATGCGATGCGATGCCACGGTGG"
    truth_counter = Counter(Kmer(test_str[i:i+5]) for i in range(len(test_str) - 5 + 1))

    counter = KmerCounter(5, 3).count_kmers(test_str)

    file_path = str(tmp_path / "test.counts")
    counter.save(file_path)

    counter_loaded = KmerCounter.from_file(file_path)
    for kmer, truth in truth_counter.items():
        assert counter_loaded.query(kmer) == truth


def test_freq_spectrum():
    test_str = "AAAAAAAA"
    counter = KmerCounter(5, 3).count_kmers(test_str)
    spectrum = counter.frequency_spectrum()

    assert len(spectrum) == counter.max_count
    assert spectrum[0] == 0
    assert spectrum[1] == 0
    assert spectrum[2] == 0
    assert spectrum[3] == 1

    truth_counter = Counter(Kmer(test_str[i:i+5]) for i in range(len(test_str) - 5 + 1))
    truth_spectrum = Counter(truth_counter.values())

    for freq in sorted(truth_spectrum):
        num_kmers = truth_spectrum[freq]

        assert spectrum[freq-1] == num_kmers

    # More complicated test sequence
    test_str = "ACTGATTTCGATGCGATGCGATGCCACGGTGG"

    counter = KmerCounter(5, 3).count_kmers(test_str)
    spectrum = counter.frequency_spectrum()

    truth_counter = Counter(Kmer(test_str[i:i+5]) for i in range(len(test_str) - 5 + 1))
    truth_spectrum = Counter(truth_counter.values())

    for freq in sorted(truth_spectrum):
        num_kmers = truth_spectrum[freq]

        assert spectrum[freq-1] == num_kmers
