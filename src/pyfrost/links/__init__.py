"""
Link De Bruijn Graph support on top of Bifrost
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Union
from pathlib import Path

import pyfrostcpp

import pyfrost.links.jt
from pyfrost.links.db import *
from pyfrost.links.nav import *

if TYPE_CHECKING:
    from pyfrost.graph import BifrostDiGraph


def add_links_from_sequence(graph: BifrostDiGraph, linkdb: LinkDB, sequence: str):
    """
    Thread `sequence` through the ccDBG (`graph`), and annotate which branches were taken at junctions in `linkdb`.

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


def add_links_from_fasta(graph: BifrostDiGraph, linkdb: LinkDB, files: Union[str, Path, list[str], list[Path]],
                         batch_size: int=1000):
    """
    Thread sequences read from the given file(s) through the graph, and annotate which branches were taken as links
    in `linkdb`. Supports compressed FASTA files too.

    Parameters
    ----------
    graph : BifrostDiGraph
        The ccDBG
    linkdb : LinkDB
        The link database to store link annotations
    files : str, Path, list
        One or more paths to FASTA files to read
    batch_size : int
        Number of sequences to read per batch. Defaults to 1000.
    """

    if not isinstance(files, list):
        files = [files]

    files = [p if isinstance(p, str) else str(p) for p in files]

    return pyfrostcpp.add_links_from_files(graph._ccdbg, linkdb, files, batch_size)
