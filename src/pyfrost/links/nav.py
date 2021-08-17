"""
Tools to aid in link-informed navigation of a ccDBG.

Think of building paths through the graph, depth-first search, etc.
"""

from __future__ import annotations
import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Iterable, Union, Optional
from pyfrost.links import jt

import networkx

if TYPE_CHECKING:
    from collections.abc import Callable
    from pyfrost import BifrostDiGraph, Kmer
    from pyfrost.links.db import LinkDB
    from pyfrost.links.jt import JunctionTreeNode

__all__ = ['oldest_link', 'link_supported_path_from']
logger = logging.getLogger(__name__)


@dataclass(order=True)
class PickedUpLink:
    dist: int
    choices: str = field()
    pos: int = field(compare=False)

    def add_dist(self, dist: int):
        self.dist += dist

    def inc_pos(self):
        self.pos += 1

    def current_choice(self) -> str:
        return self.choices[self.pos]

    def __iter__(self):
        return iter((self.dist, self.choices, self.pos))


def oldest_link(picked_up_links: list[PickedUpLink]) -> Optional[str]:
    if not picked_up_links:
        return

    oldest_links = set()
    oldest = picked_up_links[0].dist
    for dist, choices, pos in picked_up_links:
        if dist < oldest:
            break

        oldest_links.add(choices[pos])

    if len(oldest_links) == 1:
        return next(iter(oldest_links))


def link_supported_path_from(G: BifrostDiGraph, linkdb: LinkDB, source: Union[Kmer, str],
                             strategy: Callable[[list[PickedUpLink]], Optional[str]]=oldest_link,
                             distance_limit: int=None) -> Iterable[Kmer]:
    """
    Builds the longest contiguous path possible as supported by the links. Halts at branch point where the links don't
    provide a clear choice.

    Optionally, you can specify the maximum length of the path (in number of *kmers*, not unitigs), and the
    algorithm will halt when that distance is reached.

    Parameters
    ----------
    G : BifrostDiGraph
        The compacted colored De Bruijn graph to navigate
    linkdb : LinkDB
        Database with link information
    source : Kmer, str
        Source node to build path from
    strategy : callable
        A function that decides which junction choices are supported by the links. The default is to pick the branch
        that is supported by the oldest link.
    distance_limit : int
        Maximum path length in k-mers

    Yields
    ------
    Kmer
        Nodes along the path from source, including the source node itself
    """

    if isinstance(source, str):
        source = Kmer(source)

    if not source in G:
        raise networkx.NetworkXError(f"Node '{source}' not found in the graph")

    picked_up_links: list[PickedUpLink] = []

    next_node = source, 0
    while next_node:
        n, distance = next_node
        ndata = G.nodes[n]
        unitig_size = ndata['length']
        tail = ndata['tail']
        yield n

        # Pick up any links associated with this node
        if tail in linkdb:
            picked_up_links.extend(
                PickedUpLink(0, jt.junction_choices(jt_node), 0)
                for jt_node in jt.preorder(linkdb[tail]) if jt_node.is_leaf()
            )

        if not picked_up_links:
            # No link support anymore, quit
            return

        # Key neighbors by the nucleotide added
        neighbors_dict = {
            str(neighbor)[-1]: neighbor for neighbor in G[n]
        }

        # Which neighbors are supported by the links?
        if len(neighbors_dict) == 0:
            next_node = None  # No neighbors anymore, quit
        elif len(neighbors_dict) == 1:
            neighbor = next(iter(neighbors_dict.values()))
            next_node = neighbor, distance + unitig_size
        else:
            choice = strategy(picked_up_links)

            if choice:
                if choice not in neighbors_dict:
                    raise ValueError(
                        f"Mismatched graph and links! Current node: {n!r} (tail: {tail})\nNeighbors: {neighbors_dict}\n"
                        f"Link choice: {choice}\nPicked up links: {picked_up_links}"
                    )

                neighbor = neighbors_dict[choice]
                next_node = neighbor, distance + unitig_size
                picked_up_links = _move_and_prune_links(picked_up_links, choice, unitig_size)
            else:
                # No oldest link found, or multiple supported choices. We stop.
                next_node = None

        # Quit if too far from source
        if distance_limit is not None and next_node is not None:
            if next_node[1] >= distance_limit:
                next_node = None


def _move_and_prune_links(picked_up_links: list[PickedUpLink], choice: str, unitig_size: int) -> list[PickedUpLink]:
    new_picked_up_links: list[PickedUpLink] = []
    for link in picked_up_links:
        if link.pos + 1 >= len(link.choices) or link.current_choice() != choice:
            # End of link or incompatible with junction choice
            continue

        link.inc_pos()
        link.add_dist(unitig_size)
        new_picked_up_links.append(link)

    return new_picked_up_links
