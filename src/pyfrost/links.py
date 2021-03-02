from __future__ import annotations
from collections.abc import Mapping

from pyfrostcpp import LinkedDBG as _LinkedDBG
from pyfrost.seq import Kmer

__all__ = ['LinkedDBG']


class LinkedDBG(Mapping):
    """
    Add long-range connectivity information to the De Bruijn graph using Junction Trees

    Notes
    -----
    See also Turner I, Garimella KV, Iqbal Z, McVean G. Integrating long-range connectivity information into
    de Bruijn graphs. Bioinformatics. 2018 Aug 1;34(15):2556–65.
    """

    def __init__(self, graph: pyfrost.BifrostDiGraph):
        if graph._ccdbg is None:
            raise ValueError("Trying to instantiate LinkedDBG on an empty graph.")

        self._links = _LinkedDBG(graph._ccdbg)

    def add_links_from_seq(self, sequence: str):
        return self._links.add_links_from_seq(sequence)

    def get_links(self, kmer):
        return self._links.get_links(kmer)

    def __getitem__(self, kmer):
        if isinstance(kmer, str):
            kmer = Kmer(kmer)

        return self._links[kmer]

    def __len__(self):
        return len(self._links)

    def __iter__(self):
        return iter(self._links)