"""
Tools to aid in link-informed navigation of a ccDBG.

Think of building paths to the graph, depth-first search, etc.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, NamedTuple, List, Mapping

if TYPE_CHECKING:
    from pyfrost import BifrostDiGraph, Kmer
    from pyfrost.links.jt import JunctionTreeNode


class PickedUpLink(NamedTuple):
    link: str
    pos: int
    age: int


class LinkNavHistory:
    def __init__(self):
        self.links: List[PickedUpLink] = []


def link_informed_navigation(g: BifrostDiGraph, linkdb: Mapping[Kmer, JunctionTreeNode], strategy=None):
    pass
