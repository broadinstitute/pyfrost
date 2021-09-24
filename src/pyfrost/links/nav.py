"""
Tools to aid in link-informed navigation of a ccDBG.

Think of building paths through the graph, depth-first search, etc.
"""

from __future__ import annotations
import logging
import sys
from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Union, Optional
from collections.abc import Callable, Sequence

import networkx

from pyfrost.seq import Kmer
from pyfrost.links import jt
from pyfrost.exceptions import PyfrostInvalidNodeError, PyfrostInvalidPathError, PyfrostLinkMismatchError

if TYPE_CHECKING:
    from pyfrost.graph import BifrostDiGraph, Node
    from pyfrost.links.db import LinkDB

__all__ = ['follow_link', 'link_junction_edges', 'NavigationEngine', 'oldest_link', 'link_supported_path_from',
           'PickedUpLink']

logger = logging.getLogger(__name__)


def follow_link(g: BifrostDiGraph, start_kmer: Kmer, link: jt.Link) -> Iterable[Kmer]:
    """
    Returns the path through the graph represented by `link`.

    This differs from `link_supported_path` in that it just follows one link, and stops when `link` doesn't provide
    any junction choices anymore. `link_supported_path` can pickup multiple other links on the way, and only stops
    when there's ambiguity.

    Parameters
    ----------
    g : BifrostDiGraph
        The graph to navigate
    start_kmer : Kmer
        From which node the link starts
    link : jt.Link
        The link to follow

    Yields
    ------
    Kmer
        Nodes visited by following `link`.
    """

    start_unitig = g.find(start_kmer)
    if not start_unitig:
        raise networkx.NodeNotFound(f"Could not find start kmer '{start_kmer}'")

    curr_node = start_unitig['head']
    pos = 0
    while pos < len(link.choices):
        yield curr_node

        # Find neighbor
        if g.out_degree[curr_node] > 1:
            neighbors = {
                str(neighbor)[-1]: neighbor for neighbor in g.neighbors(curr_node)
            }
            choice = link.choices[pos]
            if choice not in neighbors:
                raise PyfrostLinkMismatchError(
                    f"Link junction choice is not a valid neighbor in the graph! Are links and graph mismatched?\n"
                    f"Current node: {curr_node}, junction choice: {choice}, neighbors: {neighbors}"
                )

            curr_node = neighbors[choice]
            pos += 1
        elif g.out_degree[curr_node] == 1:
            curr_node = next(g.neighbors(curr_node))
        else:
            if len(link.choices) - pos > 1:
                logger.warning(f"Link has remaining junction choices but no outgoing edges available for node "
                               f"{curr_node}")

            return


EdgeWithCov = tuple[tuple[Kmer, Kmer], int]


def link_junction_edges(g: BifrostDiGraph, start_kmer: Kmer, link: jt.Link) -> Iterable[EdgeWithCov]:
    """
    Returns the edges taken represented by the junction choices encoded in a link.

    Parameters
    ----------
    g : BifrostDiGraph
        The graph to navigate
    start_kmer : Kmer
        From which node the link starts
    link : jt.Link
        The link to follow

    Yields
    ------
    tuple[tuple[Kmer, Kmer], int]
        Tuple with the edge (head, tail), and the link coverage of the edge.
    """

    start_unitig = g.find(start_kmer)
    if not start_unitig:
        raise PyfrostInvalidNodeError(f"Could not find start kmer '{start_kmer}'")

    curr_node = start_unitig['head']
    pos = 0
    while pos < len(link.choices):
        print("Current node:", curr_node, file=sys.stderr)
        # Find neighbor
        if g.out_degree[curr_node] > 1:
            neighbors = {
                str(neighbor)[-1]: neighbor for neighbor in g.neighbors(curr_node)
            }
            choice = link.choices[pos]
            if choice not in neighbors:
                raise PyfrostLinkMismatchError(
                    f"Link junction choice is not a valid neighbor in the graph! Are links and graph mismatched?\n"
                    f"Current node: {curr_node}, junction choice: {choice}, neighbors: {neighbors}"
                )

            yield (curr_node, neighbors[choice]), link.coverage[pos]

            curr_node = neighbors[choice]
            pos += 1
        elif g.out_degree[curr_node] == 1:
            curr_node = next(g.neighbors(curr_node))
        else:
            if len(link.choices) - pos > 1:
                logger.warning(f"Link has remaining junction choices but no outgoing edges available for node "
                               f"{curr_node}")

            return


@dataclass
class PickedUpLink:
    source: Kmer
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


Strategy = Callable[[list[PickedUpLink]], Optional[str]]


class NavigationEngine:
    """
    This class implements all the logic to traverse a graph with support from links
    """

    def __init__(self, g: BifrostDiGraph, linkdb: LinkDB, link_color: int = None, strategy: Strategy = None):
        self.g: BifrostDiGraph = g
        self.linkdb: LinkDB = linkdb
        self.link_color: int = link_color

        self.strategy = strategy if strategy else oldest_link

        self.picked_up_links: list[PickedUpLink] = []
        self.used_links: set[tuple[Kmer, str]] = set()

    def _pick_up_links(self, node: Kmer, tail: Kmer):
        if tail not in self.linkdb:
            return

        # New links should at least provide more information than the current oldest
        min_length = self.picked_up_links[0].remaining() + 1 if self.picked_up_links else None
        to_add = []
        for link in jt.get_all_links(self.linkdb[tail], min_length):
            # Only use each link once
            if (tail, link.choices) in self.used_links:
                continue

            to_add.append(PickedUpLink(node, 0, link, 0))
            self.used_links.add((tail, link.choices))

        logger.debug("Picked up %d links", len(to_add))
        self.picked_up_links.extend(sorted(to_add, reverse=True))

    def _find_next_node(self, curr_node: Kmer, neighbors_dict: dict[str, Kmer]):
        # Which neighbors are supported by the links?
        if len(neighbors_dict) == 0:
            return
        elif len(neighbors_dict) == 1:
            neighbor = next(iter(neighbors_dict.values()))
            if neighbor == curr_node:
                logger.debug("The only successor of %s is itself, stopping to prevent infinite loop", curr_node)
                return
            else:
                return neighbor
        else:
            choice = self.strategy(self.picked_up_links)

            if choice:
                if choice not in neighbors_dict:
                    raise PyfrostLinkMismatchError(
                        f"Mismatched graph and links! Current node: {curr_node!r}\nNeighbors:"
                        f" {neighbors_dict}\nLink choice: {choice}\nPicked up links: {self.picked_up_links}"
                    )

                neighbor = neighbors_dict[choice]
                return neighbor

    def traverse(self, source: Union[Node, Sequence[Node]], distance_limit: int = None,
                 stop_unitig: Union[Kmer, set[Kmer]] = None) -> Iterable[Kmer]:
        """
        Builds the longest contiguous path possible as supported by the links. Halts at branch point where the links
        don't provide a clear choice.

        Optionally, you can specify the maximum length of the path (in number of *kmers*, not unitigs), and the
        algorithm will halt when that distance is reached.

        Parameters
        ----------
        source : Kmer, str, list[Kmer], list[str]
            Source node to build path from. Additionally, you can specify an existing path as source. The algorithm
            will pickup any links along that path, to gain as much "link context", which will aid in traversal beyond
            the end of the given path.
        distance_limit : int, optional
            Maximum path length in k-mers, excludes any source nodes
        stop_unitig : Kmer, set[Kmer], optional
            Stop traversal when reaching the given unitig.

        Yields
        ------
        Kmer
            Nodes along the path from source, excluding the source node(s) itself
        """

        if not isinstance(source, Sequence):
            source = [source]

        source = [Kmer(n) if not isinstance(n, Kmer) else n for n in source]

        if not stop_unitig:
            stop_unitig = set()
        elif not isinstance(stop_unitig, set):
            stop_unitig = {stop_unitig}

        self.picked_up_links = []
        self.used_links = set()

        n = source[0]
        ndata = self.g.nodes[n]
        source_ix = 0
        curr_dist = 0
        while n:
            tail = ndata['tail']
            logger.debug("Current node: %s (tail: %s)", n, tail)

            # Pick up any links associated with this node
            self._pick_up_links(n, tail)

            # Key neighbors by the nucleotide added
            neighbors_dict = {
                str(neighbor)[-1]: neighbor for neighbor in self.g.color_restricted_successors(n, self.link_color)
            }

            # Follow source path first before using the links
            if source_ix < len(source) - 1:
                source_ix += 1
                choice = str(source[source_ix])[-1]

                if choice not in neighbors_dict:
                    raise PyfrostInvalidPathError(
                        f"Given source path is not a valid path through the graph!\n"
                        f"Current node: {n}, neighbors: {neighbors_dict}, expected neighbor: {choice} "
                        f"({source[source_ix]})."
                    )

                n = source[source_ix]
                ndata = self.g.nodes[n]

                if len(neighbors_dict) > 1:
                    self._move_and_prune_links(choice)

                self._increase_distances(ndata['length'])
            else:
                n = self._find_next_node(n, neighbors_dict)
                if n:
                    ndata = self.g.nodes[n]
                    curr_dist += ndata['length']

                    if len(neighbors_dict) > 1:
                        self._move_and_prune_links(str(n)[-1])

                    self._increase_distances(ndata['length'])

                    if distance_limit and curr_dist > distance_limit:
                        n = None
                    else:
                        yield n

                    if n in stop_unitig:
                        n = None

    def _increase_distances(self, dist: int):
        for picked_up_link in self.picked_up_links:
            picked_up_link.add_dist(dist)

    def _move_and_prune_links(self, choice: str):
        new_picked_up_links: list[PickedUpLink] = []
        for picked_up_link in self.picked_up_links:
            if picked_up_link.pos + 1 >= len(picked_up_link.link.choices) or picked_up_link.current_choice() != choice:
                # End of link or incompatible with junction choice
                continue

            picked_up_link.inc_pos()
            new_picked_up_links.append(picked_up_link)

        self.picked_up_links = new_picked_up_links


def link_supported_path_from(g: BifrostDiGraph, linkdb: LinkDB, source: Union[Node, Sequence[Node]],
                             link_color: int = None, distance_limit: int = None,
                             stop_unitig: Union[Kmer, set[Kmer]] = None) -> Iterable[Kmer]:
    """
    Shortcut for building a link supported path with default navigation settings and strategy.
    """

    engine = NavigationEngine(g, linkdb, link_color)
    yield from engine.traverse(source, distance_limit, stop_unitig)
