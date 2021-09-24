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

__all__ = ['preorder', 'postorder', 'bfs',  'Link', 'get_all_links']


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


@total_ordering
class Link(NamedTuple):
    node: JunctionTreeNode
    choices: str
    coverage: numpy.ndarray

    def __eq__(self, other):
        if not isinstance(other, Link):
            raise NotImplemented

        return self.choices == other.choices

    def __lt__(self, other: Link):
        if not isinstance(other, Link):
            raise NotImplemented

        return (self.coverage[0], len(self.choices)) < (other.coverage[0], len(other.choices))

    def __le__(self, other: Link):
        if not isinstance(other, Link):
            raise NotImplemented

        return (self.coverage[0], len(self.choices)) <= (other.coverage[0], len(other.choices))


def get_all_links(jt: JunctionTreeNode, min_length: int=None) -> Iterable[Link]:
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
                yield Link(node, node.junction_choices(), node.coverages())
    else:
        for node in preorder(jt):
            if node.is_leaf():
                yield Link(node, node.junction_choices(), node.coverages())
