from collections import defaultdict

import pytest  # noqa
import networkx as nx

from pyfrost import Kmer


def test_dfs(mccortex):
    g = mccortex

    start = Kmer('ACTGA')

    visited = defaultdict(int)
    visited[start] = 1

    for edge in nx.dfs_edges(g, start):
        visited[edge[1]] += 1

    for n in visited:
        assert visited[n] == 1


def test_bfs(mccortex):
    g = mccortex

    start = Kmer('ACTGA')

    visited = defaultdict(int)
    visited[start] = 1

    for edge in nx.bfs_edges(g, start):
        visited[edge[1]] += 1

    for n in visited:
        assert visited[n] == 1

