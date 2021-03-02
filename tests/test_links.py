from collections import deque
import pytest

from pyfrost import LinkedDBG, Kmer, jt


def test_add_links(mccortex):
    linked_dbg = LinkedDBG(mccortex)

    linked_dbg.add_links_from_seq("TTTCGATGCGATGCGATGCCACG")

    known_lengths = {
        Kmer("TTCGA"): 1,
        Kmer("TCGAT"): 1,
        Kmer("GCGAT"): 2
    }

    for kmer, tree in linked_dbg.items():
        if kmer in known_lengths:
            assert len(tree) == known_lengths[kmer]
        else:
            assert len(tree) == 0

    known_choices = {
        Kmer("TTCGA"): {"TGGCG"},
        Kmer("TCGAT"): {"GGCG"},
        Kmer("GCGAT"): {"GCG", "CG"}
    }

    for kmer, tree in linked_dbg.items():
        all_choices = set()
        for node in jt.postorder(tree):
            if node.is_terminal():
                choices = deque()
                curr = node
                while curr.parent_edge:
                    choices.appendleft(curr.parent_edge)
                    curr = curr.parent

                all_choices.add("".join(choices))

        assert all_choices == known_choices.get(kmer, set())
