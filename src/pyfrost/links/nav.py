"""
Tools to aid in link-informed navigation of a ccDBG.

Think of building paths through the graph, depth-first search, etc.
"""

from __future__ import annotations
import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Iterable, Union, Optional

import networkx

from pyfrost.links import jt

if TYPE_CHECKING:
    from collections.abc import Callable
    from pyfrost import BifrostDiGraph, Kmer
    from pyfrost.links.db import LinkDB

__all__ = ['oldest_link', 'link_supported_path_from', 'PickedUpLink']

logger = logging.getLogger(__name__)


@dataclass
class PickedUpLink:
    dist: int
    link: jt.Link
    pos: int

    def add_dist(self, dist: int):
        self.dist += dist

    def inc_pos(self):
        self.pos += 1

    def current_choice(self) -> str:
        return self.link.choices[self.pos]

    def current_cov(self) -> int:
        return self.link.coverage[self.pos]

    def remaining(self) -> int:
        return len(self.link.choices) - self.pos

    def __iter__(self):
        return iter((self.dist, self.link, self.pos))

    def __eq__(self, o):
        if not isinstance(o, PickedUpLink):
            return NotImplemented

        return self.dist == o.dist and self.pos == o.pos and self.link == o.link

    def __lt__(self, o):
        if not isinstance(o, PickedUpLink):
            return NotImplemented

        return (self.dist, self.current_cov(), self.remaining()) < (o.dist, o.current_cov(), o.remaining())

    def __le__(self, o):
        if not isinstance(o, PickedUpLink):
            return NotImplemented

        return (self.dist, self.current_cov(), self.remaining()) <= (o.dist, o.current_cov(), o.remaining())

    def __gt__(self, o):
        if not isinstance(o, PickedUpLink):
            return NotImplemented

        return not self.__le__(o)

    def __ge__(self, o):
        if not isinstance(o, PickedUpLink):
            return NotImplemented

        return not self.__lt__(o)


def oldest_link(picked_up_links: list[PickedUpLink]) -> Optional[str]:
    if not picked_up_links:
        return

    oldest_links = set()
    oldest = picked_up_links[0].dist
    logger.debug("Oldest link distance: %d", oldest)
    for dist, link, pos in picked_up_links:
        if dist < oldest:
            break

        oldest_links.add(link.choices[pos])

    if len(oldest_links) == 1:
        logger.debug("%s%s, cov: %d, remaining: %d", "-" * picked_up_links[0].pos, "V",
                     picked_up_links[0].current_cov(), picked_up_links[0].remaining())
        logger.debug(picked_up_links[0].link.choices)

        return next(iter(oldest_links))
    else:
        logger.debug("Multiple oldest links: %s, halting navigation.", oldest_links)


def link_supported_path_from(G: BifrostDiGraph, linkdb: LinkDB, source: Union[Kmer, str],
                             strategy: Callable[[list[PickedUpLink]], Optional[str]]=oldest_link,
                             distance_limit: int=None, stop_unitig: Optional[Kmer]=None) -> Iterable[Kmer]:
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
    distance_limit : int, optional
        Maximum path length in k-mers
    stop_unitig : Kmer, optional
        Stop traversal when reaching the given unitig.

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
        logger.debug("Current node: %s (tail: %s)", n, tail)
        yield n

        if n == stop_unitig:
            break

        # Pick up any links associated with this node
        if tail in linkdb:
            # New links should at least provide more information than the current oldest
            min_length = picked_up_links[0].remaining() + 1 if picked_up_links else None

            to_add = [PickedUpLink(0, link, 0) for link in jt.pick_up_links(linkdb[tail], min_length)]
            picked_up_links.extend(sorted(to_add, reverse=True))

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
            increase_distances(picked_up_links, unitig_size)
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
                picked_up_links = move_and_prune_links(picked_up_links, choice)
                increase_distances(picked_up_links, unitig_size)
            else:
                # No oldest link found, or multiple supported choices. We stop.
                next_node = None

        # Quit if too far from source
        if distance_limit is not None and next_node is not None:
            if next_node[1] >= distance_limit:
                next_node = None


def increase_distances(picked_up_links: list[PickedUpLink], dist: int):
    for picked_up_link in picked_up_links:
        picked_up_link.add_dist(dist)


def move_and_prune_links(picked_up_links: list[PickedUpLink], choice: str) -> list[PickedUpLink]:
    new_picked_up_links: list[PickedUpLink] = []
    for picked_up_link in picked_up_links:
        if picked_up_link.pos + 1 >= len(picked_up_link.link.choices) or picked_up_link.current_choice() != choice:
            # End of link or incompatible with junction choice
            continue

        picked_up_link.inc_pos()
        new_picked_up_links.append(picked_up_link)

    return new_picked_up_links
