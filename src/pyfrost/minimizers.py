"""
Tools to obtain minimizers from a sequence
"""

from __future__ import annotations
from typing import Optional, NamedTuple, Iterable

from pyfrostcpp import Strand, Minimizer, minhash_iter, k_g

__all__ = ['MinimizerResult', 'Minimizer', 'all_minimizers']


class MinimizerResult(NamedTuple):
    """A single identified minimizer in a given sequence, with additional information
    like the current k-mer position, and the minimizer position in the sequence."""

    minimizer: Minimizer
    kmer_pos: int
    minimizer_pos: int
    strand: Strand


def all_minimizers(sequence: str, k: Optional[int] = None, g: Optional[int] = None) -> Iterable[MinimizerResult]:
    """
    Yield all unique minimizers in a sequence.

    Always yields the lexigraphically smallest minimizer, i.e., it could be located on
    the forward or reverse strand.
    """

    # Because minhash_iter is a C++ function, just passing None values for k or g
    # will result in a type error (the type is size_t). So that's why we do more
    # convoluted overloading.
    curr_k, curr_g = k_g()
    if k is None and g is None:
        mh_iter = iter(minhash_iter(sequence))
    elif k is not None and g is None:
        mh_iter = iter(minhash_iter(sequence, k))
    elif k is None and g is not None:
        mh_iter = iter(minhash_iter(sequence, curr_k, g))
    else:
        mh_iter = iter(minhash_iter(sequence, k, g))

    curr_k, curr_g = k_g()
    last_pos = None

    for kmer_pos, minhash_results in enumerate(mh_iter):
        for result in minhash_results:
            if last_pos is not None and result.pos <= last_pos:
                continue

            minimizer = Minimizer(sequence[result.pos:result.pos+curr_g])
            minimizer_rep = minimizer.rep()
            strand = Strand.FORWARD if minimizer == minimizer_rep else Strand.REVERSE

            yield MinimizerResult(minimizer_rep, kmer_pos, result.pos, strand)

            last_pos = result.pos
