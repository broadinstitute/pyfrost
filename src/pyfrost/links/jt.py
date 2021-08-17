"""
This module implements some basic tree traversal algorithms,
to be used with our JunctionTreeNode data structure.

Notes
-----
For reference see https://en.wikipedia.org/wiki/Tree_traversal.
The algorithms have been modified to support non-binary trees.
"""

from __future__ import annotations
from collections import deque
from typing import TYPE_CHECKING, Iterable

if TYPE_CHECKING:
    from pyfrostcpp import JunctionTreeNode


def preorder(jt: JunctionTreeNode) -> Iterable[JunctionTreeNode]:
    if not jt:
        return

    stack = deque([jt])

    while stack:
        node = stack.popleft()

        yield node

        # Reversed to that the right most child gets processed last
        for child in reversed(list(node.values())):
            stack.appendleft(child)


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