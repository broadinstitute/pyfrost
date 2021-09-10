"""
:mod:`pyfrost.graph` - Efficient colored compacted De Bruijn graph datastructures
=================================================================================

"""

from __future__ import annotations
from typing import Union, Iterable

import networkx

from pyfrostcpp import PyfrostCCDBG, NodesDict, NodeDataDict, AdjacencyOuterDict, AdjacencyType, Kmer
from pyfrost.views import PyfrostAdjacencyView, PyfrostNodeView

__all__ = ['BifrostDiGraph', 'NodesDict', 'NodeDataDict', 'AdjacencyOuterDict', 'AdjacencyType']


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

    def find(self, kmer: Union[str, Kmer], extremities_only=False):
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
