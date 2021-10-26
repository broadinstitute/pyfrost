import sys

import numpy
import pytest

from pyfrost import Kmer, links, reverse_complement


def get_known_links():
   return {
       Kmer("TTCGA"): {b"TGGCG"},
       Kmer("TCGAT"): {b"GGCG"},
       Kmer("GCGAT"): {b"GCG", b"CG"},
       Kmer("GCCAC"): {b"G"},
    }


def get_known_links_paired_end_same_unitig():
    return {
        Kmer("TTCGA"): {b"TGGCG"},
        Kmer("TCGAT"): {b"GGCG"},
        Kmer("GCGAT"): {b"GCG", b"CG"},
        Kmer("CGTGG"): {b"CCAA"},
        Kmer("GGCAT"): {b"CCAA"},
        Kmer("CGCAT"): {b"CAA", b"AA"},
        Kmer("ATCGA"): {b"A"},
        Kmer("GCCAC"): {b"G"},
    }


def get_known_links_paired_end_on_repeat():
    return {
        Kmer("TTCGA"): {b"TG"},
        Kmer("TCGAT"): {b"G"},
        Kmer("GCGAT"): {b"CG"},
        Kmer("CGTGG"): {b"C"},
        Kmer("GGCAT"): {b"C"},
        Kmer("CGCAT"): {b"AA"},
        Kmer("ATCGA"): {b"A"},
        Kmer("GCCAC"): {b"G"},
    }


def test_add_links(linked_mccortex):
    g, mem_link_db = linked_mccortex
    assert mem_link_db.color is None

    known_links = get_known_links()

    for kmer, tree in mem_link_db.items():
        all_choices = set()
        for node in links.jt.postorder(tree):
            if node.is_leaf():
                choices = node.junction_choices()
                all_choices.add(choices)

        assert all_choices == known_links.get(kmer, set())


def test_mapping_result(mccortex):
    g = mccortex
    linkdb = links.MemLinkDB()
    result = links.add_links_from_single_sequence(g, linkdb, "TTTCGATGCGATGCGATGCCACG")
    matches = result.matching_kmers()
    assert numpy.all(matches)
    assert result.junctions == b"TGGCG"
    assert result.start_unitig() == Kmer("ACTGA")
    assert result.end_unitig() == Kmer("CCACG")
    assert result.mapping_start == 0
    assert result.mapping_end == len(matches) - 1
    assert result.unitig_visits[Kmer("ATGCG")] == 2
    assert result.unitig_visits[Kmer("CGATG")] == 3

    # The link below has a 'sequencing error' near the end of the sequence and thus does not completely match with
    # the unitig
    linkdb = links.MemLinkDB()
    result = links.add_links_from_single_sequence(g, linkdb, "TTTCGATGCGATGCGATGCCACGGAGG")
    matches = result.matching_kmers()
    assert numpy.all(matches[:-3])
    assert not numpy.any(matches[-3:])
    assert result.junctions == b"TGGCG"
    assert result.start_unitig() == Kmer("ACTGA")
    assert result.end_unitig() == Kmer("CCACG")
    assert result.mapping_start == 0
    assert result.mapping_end == len(matches) - 4


def test_link_paired_read_same_unitig(mccortex):
    g = mccortex
    READ1 = "CTGATTTCGAT"
    READ2 = "CGTGGCATCGCATCGCATCGA"
    first_pass_db = links.MemLinkDB()
    links.add_links_from_single_sequence(g, first_pass_db, READ1)
    links.add_links_from_single_sequence(g, first_pass_db, reverse_complement(READ1))
    links.add_links_from_single_sequence(g, first_pass_db, READ2)
    links.add_links_from_single_sequence(g, first_pass_db, reverse_complement(READ2))

    linkdb = links.MemLinkDB()
    result = links.add_links_from_paired_read(g, linkdb, READ1, READ2, first_pass_db=first_pass_db)

    known_links = get_known_links_paired_end_same_unitig()
    for kmer, jt in linkdb.items():
        all_choices = set()
        for node in links.jt.preorder(jt):
            if node.is_leaf():
                choices = node.junction_choices()
                all_choices.add(choices)
        assert all_choices == known_links.get(kmer, set())


def test_link_paired_read_on_repeat(mccortex):
    # This read 1 ends on a unitig traversed multiple times, so links don't get extended by read 2
    g = mccortex
    READ1 = "CTGATTTCGATGCGATGC"
    READ2 = "CCACCGTGGCATCGCATC"
    first_pass_db = links.MemLinkDB()
    links.add_links_from_single_sequence(g, first_pass_db, READ1)
    links.add_links_from_single_sequence(g, first_pass_db, reverse_complement(READ1))
    links.add_links_from_single_sequence(g, first_pass_db, READ2)
    links.add_links_from_single_sequence(g, first_pass_db, reverse_complement(READ2))

    linkdb = links.MemLinkDB()
    result = links.add_links_from_paired_read(g, linkdb, READ1, READ2, first_pass_db=first_pass_db)

    known_links = get_known_links_paired_end_on_repeat()
    for kmer, jt in linkdb.items():
        all_choices = set()
        for node in links.jt.preorder(jt):
            if node.is_leaf():
                choices = node.junction_choices()
                all_choices.add(choices)
        print(kmer, all_choices, file=sys.stderr)
        assert all_choices == known_links.get(kmer, set())


def test_navigate_with_links(linked_mccortex):
    g, mem_link_db = linked_mccortex

    path = list(links.link_supported_path_from(g, mem_link_db, Kmer("ACTGA")))
    assert path == [
        Kmer("TCGAT"),
        Kmer("CGATG"),
        Kmer("ATGCG"),
        Kmer("CGATG"),
        Kmer("ATGCG"),
        Kmer("CGATG"),
        Kmer("ATGCC"),
        Kmer("CCACG"),
        Kmer("GTGGC"),
        Kmer("GCATC"),
    ]


def test_navigation_conflicting_links(linked_mccortex):
    g, mem_link_db = linked_mccortex
    links.add_links_from_single_sequence(g, mem_link_db, "TTTCGATGCCACG")

    path = list(links.link_supported_path_from(g, mem_link_db, Kmer("ACTGA")))
    assert path == [
        Kmer("TCGAT"),
        Kmer("CGATG")
    ]


def test_memlinkdb_file_load_save(linked_mccortex, tmp_path):
    g, mem_link_db = linked_mccortex
    mem_link_db.color = 0

    file_path = str(tmp_path / "mccortex.links")
    mem_link_db.save(file_path)

    mem_link_db2 = links.MemLinkDB.from_file(file_path)
    assert mem_link_db2.color == 0

    known_links = get_known_links()

    for kmer, tree in mem_link_db2.items():
        if kmer in known_links:
            assert len(tree) == len(known_links[kmer])
        else:
            assert len(tree) == 0

    for kmer, tree in mem_link_db2.items():
        all_choices = set()
        for node in links.jt.postorder(tree):
            if node.is_leaf():
                choices = node.junction_choices()
                all_choices.add(choices)

        assert all_choices == known_links.get(kmer, set())


def test_links_from_file(mccortex):
    g = mccortex
    mem_link_db = links.MemLinkDB()

    links.add_links_from_fasta(g, mem_link_db, "data/mccortex_links.fasta")

    known_links = get_known_links()
    for kmer, tree in mem_link_db.items():
        if kmer in known_links:
            assert len(tree) == len(known_links[kmer])
        else:
            assert len(tree) == 0

    for kmer, tree in mem_link_db.items():
        all_choices = set()
        for node in links.jt.postorder(tree):
            if node.is_leaf():
                choices = node.junction_choices()
                all_choices.add(choices)

        assert all_choices == known_links.get(kmer, set())
