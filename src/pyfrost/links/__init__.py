"""
Link De Bruijn Graph support on top of Bifrost
"""

from __future__ import annotations
from typing import TYPE_CHECKING

import pyfrostcpp

import pyfrost.links.jt
from pyfrost.links.db import *
from pyfrost.links.nav import *

if TYPE_CHECKING:
    from pyfrost.graph import BifrostDiGraph


def add_links_from_sequence(linkdb: LinkDB, graph: BifrostDiGraph, sequence: str):
    """
    Thread `sequence` through the ccDBG (`graph`), and annotate which branch was taken at junctions in `linkdb`.

    Parameters
    ----------
    linkdb : LinkDB
        The link database object
    graph : BifrostDiGraph
        The ccDBG
    sequence : str
        The sequence to thread to the graph
    """

    return pyfrostcpp.add_links_from_sequence(linkdb, graph._ccdbg, sequence)
