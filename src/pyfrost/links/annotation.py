"""
:mod:`pyfrost.links.annotation` - Link annotation generation and utilities.
===========================================================================

This module contains functions and classes that annotate a (compacted, colored) de Bruijn graph with links [TURNER]_.
The main work horse is `LinkAnnotator`, which k-merizes a sequence and adds link annotations to the right nodes in a
graph. Additionally, there's `RefLinkAnnator`, which is almost the same as `LinkAnnotator`, except it always
annotates the node containing the first k-mer of a sequence. This makes it possible to always reconstruct the
original genome if you know the start and end k-mer.

Other helper functions include adding links from FASTA/FASTQ files, and more.
"""

from __future__ import annotations

import logging
from pathlib import Path
from enum import Enum, auto
from typing import TYPE_CHECKING, Union, Type, NamedTuple, TextIO

import numpy
import pyfrostcpp
from pyfrostcpp import LinkAnnotator, ColorAssociatedAnnotator, MappingResult, Strand
from pyfrost.links.db import LinkDB, MemLinkDBWithCov, prune_db
from pyfrost.links.nav import link_supported_path_from
from pyfrost.io import read_fastq, read_paired_fastq, open_compressed

if TYPE_CHECKING:
    from pyfrost import BifrostDiGraph

__all__ = ['add_links_from_single_sequence', 'add_links_from_fasta',
           'add_links_from_fastq', 'add_links_from_paired_read',
           'LinkAnnotator', 'ColorAssociatedAnnotator', 'MappingResult']

logger = logging.getLogger(__name__)


class PairedAnnotationResult(Enum):
    UNKNOWN = auto()

    SAME_UNITIG = auto()
    PATH_FOUND = auto()

    REPETITIVE_UNITIG = auto()
    NO_PATH_FOUND = auto()
    READ1_MISMATCH = auto()
    READ2_MISMATCH = auto()


class ReadMappingResult(NamedTuple):
    forward: MappingResult
    reverse: MappingResult


class PairedReadMappingResult(NamedTuple):
    forward: tuple[MappingResult, MappingResult, PairedAnnotationResult]
    reverse: tuple[MappingResult, MappingResult, PairedAnnotationResult]


def add_links_from_single_sequence(graph: BifrostDiGraph, db: LinkDB, seq: str,
                                   annotator_cls: Type[LinkAnnotator] = LinkAnnotator,
                                   max_link_length: int = None) -> MappingResult:
    """
    Thread `seq` through the ccDBG (`graph`), and annotate which branches were taken at junctions in `linkdb`.

    Parameters
    ----------
    graph : BifrostDiGraph
        The ccDBG
    db : LinkDB
        The link database object
    seq : str
        The sequence to thread to the graph
    annotator_cls : LinkAnnotator
        Which annotator class to use, options are `LinkAnnotator` or `ColorAssociatedAnnotator`, which only adds links
        for a specific color in the graph.
    max_link_length : int
        Stop adding junction choices to nodes to be annotated after a certain distance in k-mers.

    Returns
    -------
    MappingResult
        An object describing how well the given sequence mapped to the graph.
    """
    annotator = annotator_cls(graph._ccdbg, db)

    if max_link_length:
        annotator.max_link_length = max_link_length

    return annotator.add_links_from_sequence(seq)


def add_links_from_paired_read(g: BifrostDiGraph, db: LinkDB, read1: str, read2: str,
                               first_pass_db: LinkDB = None,
                               read2_orientation: Strand = Strand.REVERSE):
    read1_rc = pyfrostcpp.reverse_complement(read1)
    read2_rc = pyfrostcpp.reverse_complement(read2)

    if read2_orientation == Strand.REVERSE:
        # First: R1 ----> <---- R2
        logger.debug("Add links R1 ---> <--- R2")
        result1 = add_links_from_paired_read_with_extension(g, db, read1, read2_rc, first_pass_db)

        # Flip them around: R2 ----> <---- R1
        logger.debug("Add links R2 ---> <--- R1")
        result2 = add_links_from_paired_read_with_extension(g, db, read2, read1_rc, first_pass_db)
    else:
        # First: R1 ----> ----> R2
        logger.debug("Add links R1 ---> ---> R2")
        result1 = add_links_from_paired_read_with_extension(g, db, read1, read2, first_pass_db)

        # Flip them around: R2 <---- <---- R1
        logger.debug("Add links R2 <--- <--- R1")
        result2 = add_links_from_paired_read_with_extension(g, db, read2_rc, read1_rc, first_pass_db)

    return PairedReadMappingResult(result1, result2)


def add_links_from_paired_read_with_extension(g: BifrostDiGraph, db: LinkDB, read1: str, read2: str,
                                              first_pass_db: LinkDB = None):
    """
    This function adds links from the given reads. It attempts to extend the links generated from read 1 with
    junction choices from read 2. If there's any ambiguity about how to extend them (k-mer mismatches, repeats,
    different alleles), links will be generated independently of each other.

    Assumes R1 and R2 are oriented in the same direction.
    """
    annotator = LinkAnnotator(g._ccdbg, db) if db.color is None else ColorAssociatedAnnotator(g._ccdbg, db)
    mapping_result = annotator.add_links_from_sequence(read1)
    read1_end_unitig = mapping_result.end_unitig()
    read1_kmer_matches = mapping_result.matching_kmers()
    if read1_end_unitig:
        logger.debug("Read 1 end unitig: %s, last k-mer match: %s", mapping_result.end_unitig(),
                     bool(read1_kmer_matches[-1]))

    continue_add_links = False

    # Let's see if we can find where read 2 starts
    try:
        read2_start_kmer = next(iter(pyfrostcpp.kmerize_str(read2)))
        umap = g.find(read2_start_kmer)
    except StopIteration:
        # No valid k-mers found
        umap = None

    if umap and db.color is not None and db.color not in umap['colors']:
        umap = None

    if umap and read1_end_unitig:
        read2_start_unitig = umap['head']
        logger.debug("Read 2 start unitig: %s", read2_start_unitig)

        # Ensure we don't end on a repetitive unitig, which could result in ambiguity when extending links from read 2
        end_unitig_visits = mapping_result.unitig_visits[read1_end_unitig]
        result = PairedAnnotationResult.UNKNOWN
        if end_unitig_visits == 1:
            if read2_start_unitig == read1_end_unitig:
                # Read 2 starts at the unitig where read 1 ended, continue adding links from read2 to the links
                # generated from read1
                continue_add_links = True
                result = PairedAnnotationResult.SAME_UNITIG
                logger.debug("End unitig is same and non-repetitive.")
            elif first_pass_db is not None:
                # Let's see if we can find a path from the unitig where read 1 ends to the unitig where read 2 starts,
                # utilizing the links from a first pass that added links without pair information.
                # However, we do start building the path from the first mapping k-mer of read 1, to get as much
                # context of where this read maps on to the graph, and pick up any links from this read that will aid
                # navigating the right path.
                # TODO: distance limit configurable
                logger.debug("search path...")
                logger.debug("Read 1 mapped path: %s", mapping_result.paths[-1])
                path = list(link_supported_path_from(g, first_pass_db, mapping_result.paths[-1],
                                                     distance_limit=1500, stop_unitig=read2_start_unitig))

                if path and path[-1] == read2_start_unitig:
                    # Found path to start of read 2
                    logger.debug("Found path to read 2 start unitig: %s->%s", read1_end_unitig, path)
                    continue_add_links = True
                    result = PairedAnnotationResult.PATH_FOUND

                    # Make sure we add junction choices of the found path to our junction trees
                    annotator.add_links_from_path([read1_end_unitig, *path])

                if not continue_add_links:
                    result = PairedAnnotationResult.NO_PATH_FOUND
        elif end_unitig_visits > 1:
            result = PairedAnnotationResult.REPETITIVE_UNITIG
    elif umap:
        result = PairedAnnotationResult.READ1_MISMATCH
    else:
        result = PairedAnnotationResult.READ2_MISMATCH

    logger.debug("Continue adding links? %s", continue_add_links)
    mapping_result2 = annotator.add_links_from_sequence(read2, keep_nodes=continue_add_links)

    return mapping_result, mapping_result2, result


def add_links_from_fasta(graph: BifrostDiGraph, linkdb: LinkDB, files: Union[str, Path, list[str], list[Path]],
                         annotator_cls: Type[LinkAnnotator] = LinkAnnotator, batch_size: int = 1000,
                         both_strands: bool = False, max_link_length: int = None):
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
    color : int
        Associate links generated from this FASTA file with a specific color. Optional.
    annotator_cls : Type[LinkAnnotator]
        The class to use for creating links. Defaults to `LinkAnnotator`, but you could for example specify
        `ColorAssociatedAnnotator` here, to only add links for a certain color.
    batch_size : int
        Number of sequences to read per batch. Defaults to 1000.
    both_strands : bool
        Whether to add links from both the forward and reverse orientations of each sequence in the FASTA file.
    max_link_length : int
        Stop adding junction choices to nodes to be annotated after a certain distance in k-mers.
    """

    if not isinstance(files, list):
        files = [files]

    files = [p if isinstance(p, str) else str(p) for p in files]

    annotator = annotator_cls(graph._ccdbg, linkdb)

    if max_link_length:
        annotator.max_link_length = max_link_length

    return pyfrostcpp.add_links_from_fasta(annotator, files, batch_size, both_strands)


def add_links_from_fastq_single(g: BifrostDiGraph, linkdb: LinkDB, file_path: Union[str, Path],
                                mapping_results_out: TextIO = None):
    annotator = LinkAnnotator(g._ccdbg, linkdb) if linkdb.color is None else ColorAssociatedAnnotator(g._ccdbg, linkdb)

    # TODO: parallelize? Need to make LinkDB thread-safe then
    with open_compressed(file_path, "rt") as fp:
        for read in read_fastq(fp):
            if len(read.seq) < g.graph['k']:
                continue

            read_name = read.name
            if read_name.endswith('/1'):
                read_name = read_name[:-2]

            logger.debug("Processing read %s...", read_name)

            result = ReadMappingResult(
                annotator.add_links_from_sequence(read.seq),
                annotator.add_links_from_sequence(pyfrostcpp.reverse_complement(read.seq))
            )

            if mapping_results_out:
                print(read_name, "F", *log_mapping_result(result.forward), sep='\t', file=mapping_results_out)
                print(read_name, "R", *log_mapping_result(result.reverse), sep='\t', file=mapping_results_out)


def add_links_from_fastq(g: BifrostDiGraph, linkdb: LinkDB,
                         files: Union[str, Path, list[str], list[Path]], twopass: bool = True,
                         read2_orientation: Strand = Strand.REVERSE, mapping_results_out: TextIO = None,
                         *, existing_first_pass: str = None, first_pass_ofile: str = None,
                         prune_first_pass: int = 0):
    if not isinstance(files, list):
        files = [files]

    if not files:
        return

    if linkdb.color is None:
        logger.info("Generating links for the whole graph (no color selected)...")
    else:
        logger.info("Generating links for color %s...", g.graph['color_names'][linkdb.color])

    if len(files) == 1:
        # Single read file
        logger.info("PASS 1 / 1: Adding links from %s (single-end mode)...", files[0])
        add_links_from_fastq_single(g, linkdb, files[0], mapping_results_out=mapping_results_out)
    else:
        # Paired-end reads
        first_pass_db = None
        if twopass:
            if existing_first_pass:
                logger.info("PASS 1 / 2: Loading existing first-pass database from %s", existing_first_pass)
                first_pass_db = MemLinkDBWithCov.from_file(existing_first_pass)
            else:
                # In the first pass we just add links independently without pair information
                logger.info("PASS 1.1 / 2: Building first-pass link database...")
                first_pass_db = MemLinkDBWithCov()
                first_pass_db.color = linkdb.color

                add_links_from_fastq_single(g, first_pass_db, files[0])
                add_links_from_fastq_single(g, first_pass_db, files[1])

                if prune_first_pass:
                    logger.info("PASS 1.2 / 2: Pruning first-pass link database...")
                    prune_db(first_pass_db, prune_first_pass)

                if first_pass_ofile:
                    logger.info("PASS 1.3 / 2: Saving first-pass link database...")
                    first_pass_db.save(first_pass_ofile)

            logger.info("PASS 2 / 2: Adding links from %s and %s (paired-end mode), utilizing the first-pass "
                        "database...", *files)
        else:
            logger.info("PASS 1 / 1: Adding links from %s and %s (single pass paired-end mode)...", *files)

        with open_compressed(files[0], "rt") as fp1, open_compressed(files[1], "rt") as fp2:
            for read1, read2 in read_paired_fastq(fp1, fp2):
                if len(read1.seq) < g.graph['k'] or len(read2.seq) < g.graph['k']:
                    continue

                read_name = read1.name
                if read_name.endswith('/1'):
                    read_name = read_name[:-2]

                logger.debug("Processing read pair %s...", read_name)

                result = add_links_from_paired_read(g, linkdb, read1.seq, read2.seq,
                                                    first_pass_db=first_pass_db,
                                                    read2_orientation=read2_orientation)

                if mapping_results_out:
                    print(read_name, "F", str(result.forward[2]),
                          *log_mapping_result(result.forward[0]), *log_mapping_result(result.forward[1]),
                          sep='\t', file=mapping_results_out)
                    print(read_name, "R", str(result.reverse[2]),
                          *log_mapping_result(result.reverse[0]), *log_mapping_result(result.reverse[1]),
                          sep='\t', file=mapping_results_out)


def log_mapping_result(result: MappingResult):
    matches = result.matching_kmers()
    matches_str = "".join(map(str, matches))

    pct_match = (numpy.count_nonzero(matches) * 100) / len(matches)

    start_unitig = result.start_unitig()
    end_unitig = result.end_unitig()
    return [start_unitig if start_unitig else ".", end_unitig if end_unitig else ".",
            len(matches), result.mapping_start, result.mapping_end, matches_str, f"{pct_match:.2f}"]
