from collections import deque
import pytest

from pyfrost import Kmer, links


def test_add_links(linked_mccortex):
    g, mem_link_db = linked_mccortex

    known_lengths = {
        Kmer("TTCGA"): 1,
        Kmer("TCGAT"): 1,
        Kmer("GCGAT"): 2
    }

    for kmer, tree in mem_link_db.items():
        if kmer in known_lengths:
            assert len(tree) == known_lengths[kmer]
        else:
            assert len(tree) == 0

    known_choices = {
        Kmer("TTCGA"): {"TGGCG"},
        Kmer("TCGAT"): {"GGCG"},
        Kmer("GCGAT"): {"GCG", "CG"}
    }

    for kmer, tree in mem_link_db.items():
        all_choices = set()
        for node in links.jt.postorder(tree):
            if node.is_leaf():
                choices = links.jt.junction_choices(node)
                all_choices.add(choices)

        assert all_choices == known_choices.get(kmer, set())


def test_navigate_with_links(linked_mccortex):
    g, mem_link_db = linked_mccortex

    path = list(links.link_supported_path_from(g, mem_link_db, Kmer("ACTGA")))
    assert path == [
        Kmer("ACTGA"),
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
    links.add_links_from_sequence(g, mem_link_db, "TTTCGATGCCACG")

    path = list(links.link_supported_path_from(g, mem_link_db, Kmer("ACTGA")))
    assert path == [
        Kmer("ACTGA"),
        Kmer("TCGAT"),
        Kmer("CGATG")
    ]
