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


def add_links_from_sequence(graph: BifrostDiGraph, linkdb: LinkDB, sequence: str):
    """
    Thread `sequence` through the ccDBG (`graph`), and annotate which branch was taken at junctions in `linkdb`.

    Parameters
    ----------
    graph : BifrostDiGraph
        The ccDBG
    linkdb : LinkDB
        The link database object
    sequence : str
        The sequence to thread to the graph
    """

    return pyfrostcpp.add_links_from_sequence(graph._ccdbg, linkdb, sequence)
