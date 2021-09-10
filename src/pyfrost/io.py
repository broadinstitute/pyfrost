"""
:mod:`pyfrost.io` - Load/save/build colored compacted De Bruijn graphs
======================================================================
"""

from __future__ import annotations
import bz2
import gzip
import lzma
from pathlib import Path
from itertools import zip_longest
from typing import Union, TextIO, BinaryIO, Iterable, Optional, NamedTuple

import pyfrostcpp

from pyfrost.graph import BifrostDiGraph

__all__ = ['build', 'build_from_refs', 'build_from_samples', 'load', 'dump', 'open_compressed',
           'read_fastq', 'read_paired_fastq']


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


def build(refs: list[str], samples: list[str], **kwargs):
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


def build_from_refs(refs: list[str], **kwargs):
    return build(refs, [], **kwargs)


def build_from_samples(samples: list[str], **kwargs):
    return build([], samples, **kwargs)


def dump(g: BifrostDiGraph, fname_prefix: str, num_threads: int=2):
    return pyfrostcpp.dump(g._ccdbg, fname_prefix, num_threads)


def open_compressed(filename, *args, **kwargs) -> Union[TextIO, BinaryIO]:
    if not isinstance(filename, Path):
        filename = Path(filename)

    fopen = {
        '.gz': gzip.open,
        '.bz2': bz2.open,
        '.lz': lzma.open
    }.get(filename.suffix, open)

    return fopen(filename, *args, **kwargs)


class FastxRecord(NamedTuple):
    name: str
    seq: str
    qual: Optional[str]


def read_fastq(fp: TextIO) -> Iterable[FastxRecord]:
    """Heng Li's fast FASTQ reader."""

    last = None
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break

        if not last:
            break

        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])

        if not last or last[0] != '+':  # this is a fasta record
            yield FastxRecord(name, ''.join(seqs), None)  # yield a fasta record

            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield FastxRecord(name, seq, ''.join(seqs))  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield FastxRecord(name, seq, None)  # yield a fasta record instead
                break


def read_paired_fastq(fp1: TextIO, fp2: TextIO) -> Iterable[tuple[FastxRecord, FastxRecord]]:
    for read1, read2 in zip_longest(read_fastq(fp1), read_fastq(fp2)):
        if not read1 or not read2:
            raise IOError("FASTQ files contain different number of reads.")

        read1_name = read1[0]
        read2_name = read2[0]
        if read1_name.endswith('/1'):
            read1_name = read1_name[:-2]

        if read2_name.endswith('/2'):
            read2_name = read2_name[:-2]

        if read1_name != read2_name:
            raise ValueError(f"Mismatched reads! Read ID {read1_name} doesn't match {read2_name}!")

        yield read1, read2
