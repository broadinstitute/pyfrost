"""
:mod:`pyfrost.seq` - Utilities for DNA sequences, kmers, graphs, etc.
=====================================================================
"""

from typing import List

from bifrost_python import reverse_complement, Kmer, Strand, set_k

__all__ = ['reverse_complement', 'Kmer', 'Strand', 'set_k', 'path_sequence', 'path_nucleotide_length']


def path_sequence(g: 'bifrost_python.BifrostDiGraph', path: List[Kmer]) -> str:
    """Build the DNA sequence that is spelled by the given path."""

    if not path:
        return ""

    node_iter = iter(path)
    start = next(node_iter)

    k = g.graph['k']

    seqs = [g.nodes[start]['unitig_sequence']]
    prev = start
    for n in node_iter:
        if (prev, n) not in g.edges:
            raise ValueError(f"Invalid path specified, ({prev}, {n}) is not an edge.")

        seqs.append(g.nodes[n]['unitig_sequence'][k-1:])
        prev = n

    return "".join(seqs)


def path_nucleotide_length(g: 'bifrost_python.BifrostDiGraph', path: List[Kmer]) -> int:
    """Compute the length of a path in nucleotides."""

    if not path:
        return 0

    node_iter = iter(path)
    start = next(node_iter)

    k = g.graph['k']

    length = g.nodes[start]['length'] + k - 1
    prev = start
    for n in node_iter:
        if (prev, n) not in g.edges:
            raise ValueError(f"Invalid path specified, ({prev}, {n}) is not an edge.")

        length += g.nodes[n]['length']
        prev = n

    return length

