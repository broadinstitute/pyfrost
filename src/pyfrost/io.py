"""
:mod:`pyfrost.io` - Load/save/build colored compacted De Bruijn graphs
======================================================================
"""

from __future__ import annotations
from typing import Union
from pathlib import Path

import pyfrostcpp

from pyfrost.graph import BifrostDiGraph


def load(graph: Union[str, Path], **kwargs):
    """
    Load a Bifrost graph from a file.

    Automatically determines the bfg_colors path from the given GFA file.

    Parameters
    ----------
    graph : str, Path
        Path to the Bifrost GFA. Path to colors file is inferred
    kwargs
        Additional options, e.g. `nb_threads` or `verbose`

    Returns
    -------
    BifrostDiGraph
    """

    if not isinstance(graph, Path):
        graph = Path(graph)

    if not graph.is_file():
        raise IOError(f"Could not find graph GFA file {graph}")

    colors_file = graph.with_suffix(".bfg_colors")

    if not colors_file.is_file():
        raise IOError(f"Could not find graph colors file {colors_file}")

    return BifrostDiGraph(pyfrostcpp.load(str(graph), str(colors_file), **kwargs))


def build(refs, samples, **kwargs):
    """
    Build a Bifrost graph from multiple FASTA and/or FASTQ files.

    Parameters
    ----------
    refs
    samples
    kwargs

    Returns
    -------

    """
    return BifrostDiGraph(pyfrostcpp.build(refs, samples, **kwargs))


def build_from_refs(refs, **kwargs):
    return build(refs, [], **kwargs)


def build_from_samples(samples, **kwargs):
    return build([], samples, **kwargs)


def dump(g: BifrostDiGraph, fname_prefix, num_threads=2):
    return pyfrostcpp.dump(g._ccdbg, fname_prefix, num_threads)


__all__ = ['build', 'build_from_refs', 'build_from_samples', 'load', 'dump']
