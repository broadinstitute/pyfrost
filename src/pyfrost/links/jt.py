"""
This module implements some basic tree traversal algorithms,
to be used with our JunctionTreeNode data structure.

Notes
-----
For reference see https://en.wikipedia.org/wiki/Tree_traversal.
The algorithms have been modified to support non-binary trees.
"""

from __future__ import annotations
from functools import total_ordering
from collections import deque
from typing import TYPE_CHECKING, Iterable, NamedTuple

import numpy

if TYPE_CHECKING:
    from pyfrostcpp import JunctionTreeNode

__all__ = ['preorder', 'postorder', 'bfs', 'junction_choices', 'junction_coverages', 'Link', 'pick_up_links']


def preorder(jt: JunctionTreeNode, with_depth=False) -> Iterable[JunctionTreeNode]:
    if not jt:
        return

    stack = deque([(jt, 0)])

    while stack:
        node, depth = stack.popleft()

        if with_depth:
            yield node, depth
        else:
            yield node

        # Reversed to that the right most child gets processed last
        for child in reversed(list(node.values())):
            stack.appendleft((child, depth+1))


def postorder(jt: JunctionTreeNode) -> Iterable[JunctionTreeNode]:
    if not jt:
        return

    stack = deque()
    child_ix = deque()
    node = jt
    while stack or node:
        if node is not None:
            stack.appendleft(node)
            child_ix.appendleft(0)

            if len(node) > 0:
                node = next(iter(node.values()))
            else:
                node = None
        else:
            peek = stack[0]
            child = child_ix[0]
            children = list(peek.values())

            if child < len(children):
                node = children[child]
            else:
                yield peek

                stack.popleft()
                child_ix.popleft()

                if child_ix:
                    child_ix[0] += 1  # Update child ix of parent


def bfs(jt: JunctionTreeNode) -> Iterable[JunctionTreeNode]:
    if not jt:
        return

    queue = deque([jt])
    while queue:
        node = queue.popleft()
        yield node

        for c in node.values():
            queue.append(c)


def junction_choices(jt: JunctionTreeNode) -> str:
    """
    Get junction choices that lead to the given tree node represented as a string.

    Parameters
    ----------
    jt : JunctionTreeNode
        The node in the tree for which to obtain all junction choices. Usually this is a leaf node.

    Returns
    -------
    str
        Each nucleotide in the returned string represents a junction choice.
    """

    choices = []
    curr = jt
    while curr.parent_edge:
        choices.append(curr.parent_edge)
        curr = curr.parent

    return "".join(reversed(choices))


def junction_coverages(jt: JunctionTreeNode) -> numpy.ndarray:
    covs = []
    curr = jt
    while curr.parent_edge:
        covs.append(curr.count)
        curr = curr.parent

    covs.reverse()
    return numpy.array(covs, dtype=numpy.uint16)


@total_ordering
class Link(NamedTuple):
    choices: str
    coverage: numpy.ndarray

    def __lt__(self, other: Link):
        if not isinstance(other, Link):
            raise NotImplemented

        return (self.coverage, self.choices) < (other.coverage, other.choices)

    def __le__(self, other: Link):
        if not isinstance(other, Link):
            raise NotImplemented

        return (self.coverage, self.choices) <= (other.coverage, other.choices)


def pick_up_links(jt: JunctionTreeNode, min_length: int=None) -> Iterable[Link]:
    """
    Yield each link in a junction tree.

    This function traverses the tree (in preorder), and yields the list of junction choices at each leaf node.

    Parameters
    ----------
    jt : JunctionTreeNode
        Root (or other internal junction tree node) representing the junction tree
    min_length : int
        Minimum number of junction choices of a link to report it

    Yields
    ------
    Link
        Single Link object representing all junction choices and corresponding coverages
    """

    with_depth = min_length is not None
    if with_depth:
        for node, depth in preorder(jt, with_depth):
            if node.is_leaf() and depth >= min_length:
                yield Link(junction_choices(node), junction_coverages(node))
    else:
        for node in preorder(jt):
            if node.is_leaf():
                yield Link(junction_choices(node), junction_coverages(node))
