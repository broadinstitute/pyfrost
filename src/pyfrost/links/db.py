"""
Classes and functions to store and create links
"""

from pyfrostcpp import LinkDB, MemLinkDB, MemLinkDBWithCov

__all__ = ['LinkDB', 'MemLinkDB', 'MemLinkDBWithCov', 'prune_db']


def prune_db(db: LinkDB, threshold: int):
    """
    Prune links with a coverage lower than a given threshold.

    Modifies the given database in-place.
    """

    for kmer, tree in db.items():
        tree.prune(threshold)
