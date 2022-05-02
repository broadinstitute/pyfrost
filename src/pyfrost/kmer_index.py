"""
Index the unitig head k-mers of reference sequences and map their locations.

This data structure provides thus quick look-up for where unitigs are located within a
reference.
"""

from __future__ import annotations
import json
import logging
from collections import defaultdict
from typing import TextIO, NamedTuple

import skbio

from pyfrost import kmerize_str, Kmer, BifrostDiGraph

logger = logging.getLogger(__name__)


class KmerPosition(NamedTuple):
    contig_id: int
    pos: int


def json_decode(kv_pairs):
    if kv_pairs[0][0] in {"contig_ids", "index"}:
        return dict(kv_pairs)
    else:
        return {Kmer(k): KmerPosition(*v) for k, v in kv_pairs.items()}


class KmerPosIndex:
    def __init__(self):
        self.contig_ids = []
        self.index: dict[Kmer, list[KmerPosition]] = defaultdict(list)

    def index_fasta(self, g: BifrostDiGraph, f: TextIO):
        for r in skbio.io.read(f, "fasta"):
            current_contig_id = len(self.contig_ids)
            self.contig_ids.append(r.metadata['id'])

            for pos, kmer in enumerate(kmerize_str(str(r))):
                unitig = g.find(kmer, extremities_only=True)
                canonical = kmer.rep()

                if unitig:
                    self.index[canonical].append(KmerPosition(current_contig_id, pos))

    def get_positions(self, kmer: Kmer) -> list[KmerPosition]:
        return self.index.get(kmer.rep(), [])

    def save(self, ofile: TextIO):
        print(*self.contig_ids, sep='\t', file=ofile)
        for kmer, positions in self.index.items():
            pos_str = [f"{c}:{pos}" for c, pos in positions]
            print(str(kmer), *pos_str, sep='\t', file=ofile)

    def load(self, ifile: TextIO):
        first = True
        for line in ifile:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if first:
                self.contig_ids = list(parts)
                first = False
            else:
                kmer = Kmer(parts[0])
                for contig_pos in parts[1:]:
                    contig_id, pos = contig_pos.split(':')
                    self.index[kmer].append(KmerPosition(int(contig_id), int(pos)))
