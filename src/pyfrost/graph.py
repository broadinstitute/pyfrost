"""
:mod:`pyfrost.graph` - Efficient colored compacted De Bruijn graph datastructures
=================================================================================

"""

from __future__ import annotations
from typing import Union, Iterable
from collections import deque, abc

import networkx

from pyfrostcpp import PyfrostCCDBG, NodesDict, NodeDataDict, AdjacencyOuterDict, AdjacencyType, Kmer
from pyfrost.views import PyfrostAdjacencyView, PyfrostNodeView

__all__ = ['BifrostDiGraph', 'Node', 'NodesDict', 'NodeDataDict', 'AdjacencyOuterDict', 'AdjacencyType',
           'get_neighborhood']


def disable_factory_func():
    raise networkx.NetworkXError("NetworkX factory functions are disabled for Bifrost graphs, adding nodes is done "
                                 "using different functions!")


Node = Union[Kmer, str]


class BifrostDiGraph(networkx.DiGraph):
    """
    A class representing a compacted, colored De Bruijn graph, with Bifrost as backend.
    """

    node_dict_factory = disable_factory_func
    node_attr_dict_factory = disable_factory_func
    adjlist_outer_dict_factory = disable_factory_func
    adjlist_inner_dict_factory = disable_factory_func
    edge_attr_dict_factory = disable_factory_func

    def __init__(self, bifrost_ccdbg=None, **attr):
        self.graph = attr

        if bifrost_ccdbg is None:
            # Initialize with empty dicts if no PyfrostCCDBG object given.
            #
            # Many NetworkX functions (e.g. subgraph_view) rely on being able to construct a new graph object without
            # arguments to the constructor. Those functions often assign _adj, _pred, and _succ afterwards so these
            # empty dicts will often be replaced immediately.
            self._node = {}
            self._adj = {}
            self._pred = {}
            self._succ = self._adj
        else:
            self._ccdbg: PyfrostCCDBG = bifrost_ccdbg

            self.graph['color_names'] = list(bifrost_ccdbg.color_names())
            self.graph['k'] = bifrost_ccdbg.get_k()
            self.graph['g'] = bifrost_ccdbg.get_g()

            self._node = NodesDict(bifrost_ccdbg)
            self._adj = AdjacencyOuterDict(bifrost_ccdbg, AdjacencyType.SUCCESSORS)
            self._pred = AdjacencyOuterDict(bifrost_ccdbg, AdjacencyType.PREDECESSORS)
            self._succ = self._adj

    def add_node(self, node_for_adding, **attr):
        raise networkx.NetworkXError("Manually adding nodes is not supported for BifrostDiGraph")

    def add_nodes_from(self, nodes_for_adding, **attr):
        raise networkx.NetworkXError("Manually adding nodes is not supported for BifrostDiGraph")

    def remove_node(self, n: Node):
        try:
            self._ccdbg.remove(n)
        except KeyError:
            raise networkx.NetworkXError(f"Node {n} does not exists")

    def remove_nodes_from(self, nodes: Iterable[Node]):
        for n in nodes:
            try:
                self.remove_node(n)
            except networkx.NetworkXError:
                pass

    def add_edge(self, u_of_edge, v_of_edge, **attr):
        raise networkx.NetworkXError("Manually adding edges is not supported for BifrostDiGraph")

    def add_edges_from(self, ebunch_to_add, **attr):
        raise networkx.NetworkXError("Manually adding edges is not supported for BifrostDiGraph")

    def remove_edge(self, u, v):
        raise networkx.NetworkXError("Manually removing edges is not supported for BifrostDiGraph")

    def remove_edges_from(self, ebunch):
        raise networkx.NetworkXError("Manually removing edges is not supported for BifrostDiGraph")

    def find(self, kmer: Node, extremities_only=False):
        """
        Find the unitig (node) containing a given k-mer
        """

        return self._ccdbg.find(kmer, extremities_only)

    def nbunch_iter(self, nbunch=None):
        # Transform any nodes given as str to Kmer objects
        yield from (
            Kmer(n) if isinstance(n, str) else n for n in super().nbunch_iter(nbunch)
        )

    @property
    def nodes(self):
        nodes = PyfrostNodeView(self)
        # Lazy View creation: overload the (class) property on the instance
        # Then future G.nodes use the existing View
        # setattr doesn't work because attribute already exists
        self.__dict__['nodes'] = nodes
        return nodes

    @property
    def adj(self):
        return PyfrostAdjacencyView(self._succ)

    @property
    def succ(self):
        return PyfrostAdjacencyView(self._succ)

    @property
    def pred(self):
        return PyfrostAdjacencyView(self._pred)

    def color_restricted_successors(self, node: Node, allowed_colors: Union[set, int] = None):
        """Get succesor nodes with certain colors.

        If leaving `allowed_colors` to None, it will return all neighbors."""

        if allowed_colors is None:
            yield from self.successors(node)
        else:
            if not isinstance(allowed_colors, set):
                allowed_colors = {allowed_colors}

            if not isinstance(node, Kmer):
                node = Kmer(node)

            yield from (n.head for n in self._ccdbg.color_restricted_successors(node, allowed_colors))

    color_restricted_neighbors = color_restricted_successors

    def color_restricted_predecessors(self, node: Node, allowed_colors: Union[set, int] = None):
        """Get predecessor nodes with certain colors.

        If leaving `allowed_colors` to None, it will return all predecessors."""

        if allowed_colors is None:
            yield from self.predecessors(node)
        else:
            if not isinstance(allowed_colors, set):
                allowed_colors = {allowed_colors}

            if not isinstance(node, Kmer):
                node = Kmer(node)

            yield from (n.head for n in self._ccdbg.color_restricted_predecessors(node, allowed_colors))


def get_neighborhood(g: BifrostDiGraph, source, radius=3, both_directions=True, ext_colors=None):
    """
    Get the neighborhood around one or more given source nodes.

    Returns a subgraph with all neighboring nodes reachable within a maximum
    number of edge traversals (`radius`) from a source node.

    If `both_directions` is True, then it also considers predecessors of the nodes,
    otherwise it will only traverse edges in the forward direction.

    Optionally, you can give a set of colors for which the neighborhood is extended. If a node has a color other than
    the ones given in `ext_colors` it is not further extended.
    """

    neighborhood = networkx.DiGraph()
    neighborhood_nodes = set()

    queue = deque()
    if isinstance(source, abc.Iterable):
        for src_node in source:
            queue.append((0, Kmer(src_node)))
    else:
        queue.append((0, Kmer(source)))

    if ext_colors is not None and not isinstance(ext_colors, set):
        ext_colors = {ext_colors}

    while queue:
        level, node = queue.popleft()

        if node in neighborhood_nodes:
            continue

        ndata = g.nodes[node]
        neighborhood_nodes.add(node)

        # Don't extend any further if ext_colors is given and the node has a different color
        if level < radius and (ext_colors is None or ext_colors & ndata['colors']):
            for succ in g.successors(node):
                if succ in neighborhood:
                    continue

                queue.append((level + 1, succ))

            if both_directions:
                for pred in g.predecessors(node):
                    if pred in neighborhood:
                        continue

                    queue.append((level + 1, pred))

    return g.subgraph(neighborhood_nodes)
