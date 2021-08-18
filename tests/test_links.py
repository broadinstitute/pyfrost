from collections import deque
import pytest

from pyfrost import Kmer, links


def get_known_links():
   return {
        Kmer("TTCGA"): {"TGGCG"},
        Kmer("TCGAT"): {"GGCG"},
        Kmer("GCGAT"): {"GCG", "CG"}
    }


def test_add_links(linked_mccortex):
    g, mem_link_db = linked_mccortex
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
                choices = links.jt.junction_choices(node)
                all_choices.add(choices)

        assert all_choices == known_links.get(kmer, set())


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


def test_memlinkdb_file_load_save(linked_mccortex, tmp_path):
    g, mem_link_db = linked_mccortex

    file_path = str(tmp_path / "mccortex.links")
    mem_link_db.save(file_path)

    mem_link_db2 = links.MemLinkDB.from_file(file_path)
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
                choices = links.jt.junction_choices(node)
                all_choices.add(choices)

        assert all_choices == known_links.get(kmer, set())
