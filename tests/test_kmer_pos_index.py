import pyfrost
import pytest  # noqa

import skbio

from pyfrost.kmer_index import Kmer, KmerPosIndex


def test_kmer_pos_index(mccortex):
    index = KmerPosIndex()
    k = pyfrost.k_g()[0]

    with open("data/mccortex.fasta") as f:
        index.index_fasta(mccortex, f)

    with open("data/mccortex.fasta") as f:
        seqs = [str(r) for r in skbio.io.read(f, "fasta")]

    for kmer, positions in index.index.items():
        for pos in positions:
            assert kmer == Kmer(seqs[pos.contig_id][pos.pos:pos.pos+k]).rep()


def test_save_load_pos_index(mccortex, tmp_path):
    index = KmerPosIndex()

    with open("data/mccortex.fasta") as f:
        index.index_fasta(mccortex, f)

    with open(tmp_path / "index.json", "w") as o:
        index.save(o)

    index2 = KmerPosIndex()
    with open(tmp_path / "index.json") as in_file:
        index2.load(in_file)

    assert index.contig_ids == index2.contig_ids
    assert index.index == index2.index
