from collections import defaultdict

import pytest  # noqa
import networkx as nx


def test_dfs(mccortex):
    g = mccortex

    start = g.nodes['ACTGA']

    visited = defaultdict(int)
    visited[start] = 1

    for edge in nx.dfs_edges(g, start):
        visited[edge[1]] += 1

    for n in visited:
        assert visited[n] == 1

    visited_set = set(str(n) for n in visited.keys())
    truth = {
        "ACTGATTTCGA",
        "TCGAAATCAGT",
        "TCGAT",
        "CGATGC",
        "ATGCCAC",
        "CCACCGTGG",
        "CCACGGTGG",
        "GCATCG",
        "GTGGCAT",
        "ATGCGAT"
    }

    assert visited_set == truth


def test_bfs(mccortex):
    g = mccortex

    start = g.nodes['ACTGA']

    visited = defaultdict(int)
    visited[start] = 1

    for edge in nx.bfs_edges(g, start):
        visited[edge[1]] += 1

    for n in visited:
        assert visited[n] == 1

    visited_set = set(str(n) for n in visited.keys())
    truth = {
        "ACTGATTTCGA",
        "TCGAAATCAGT",
        "TCGAT",
        "CGATGC",
        "ATGCCAC",
        "CCACCGTGG",
        "CCACGGTGG",
        "GCATCG",
        "GTGGCAT",
        "ATGCGAT"
    }

    assert visited_set == truth
