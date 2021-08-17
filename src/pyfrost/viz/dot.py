"""
:mod:`pyfrost.viz.dot` - Visualize Pyfrost graph using Dot
==========================================================

"""

from __future__ import annotations
from typing import TYPE_CHECKING

import pydot

if TYPE_CHECKING:
    from pyfrost.graph import BifrostDiGraph
    from IPython.display import Image


def to_pydot(dbg: BifrostDiGraph, full_sequence: bool=False) -> pydot.Dot:
    """
    Transform a Pyfrost compacted colored De Bruijn graph to a PyDot graph data structure.

    :param dbg: The Pyfrost Graph to convert
    :param full_sequence: If true, use full unitig sequence as node label, otherwise abbreviate
    :return: PyDot graph data structure
    """

    pg = pydot.Dot('', graph_type='digraph', rankdir='LR')
    pg.set_node_defaults(shape='record')

    for n, data in dbg.nodes(data=True):
        if full_sequence:
            label = data['unitig_sequence']
        else:
            if data['unitig_length'] > (dbg.graph['k'] * 2):
                label = f"{data['head']}...{data['tail']}"
            else:
                label = data['unitig_sequence']

        pn = pydot.Node(str(n), label=label)
        pg.add_node(pn)

    for n1, n2 in dbg.edges:
        pe = pydot.Edge(str(n1), str(n2))
        pg.add_edge(pe)

    return pg


def display_graph(dbg: BifrostDiGraph, *args, **kwargs) -> Image:
    """
    Display the graph
    :param dbg:
    :param args:
    :param kwargs:
    :return:
    """
    from IPython.display import Image

    pydot_graph = to_pydot(dbg, *args, **kwargs)
    return Image(pydot_graph.create_png(prog='dot'))
