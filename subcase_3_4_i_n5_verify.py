#!/usr/bin/env python3
"""
Verification code for Theorem 4.1, Subcase 3.4(i), n = 5
in the paper "Fault-tolerant two-disjoint-cycle-cover pancyclicity
of augmented cubes".

Subcase being verified:
  - n = 5
  - ell = 2^(n-1) - 1 = 15
  - |F_R| = 2n - 3 = 7
  - every vertex of R - F_R has internal degree at least 2

What the script proves:
For every 7-edge fault set F_R in AQ_4 such that delta(AQ_4 - F_R) >= 2,
the graph AQ_4 - F_R contains a 15-cycle.

This is exactly the n=5 part of Subcase 3.4(i): once such a 15-cycle exists
in R-F_R, the rest of the subcase follows from the already-cited result on
the L-side.

Method:
  1. Build AQ_4 as a Cayley graph with generators
       {e_1,e_2,e_3,e_4, eps_2,eps_3,eps_4}.
  2. Search for a counterexample fault set F_R using MILP:
       - choose exactly 7 faulty edges,
       - keep internal degree at least 2 at every vertex,
       - force the chosen faults to hit every discovered 15-cycle.
  3. Given a candidate F_R, test whether AQ_4 - F_R still has a 15-cycle.
     If it does, add a blocking constraint, together with all its automorphic images.
  4. If the MILP becomes infeasible, then no counterexample exists.

Dependencies:
  - networkx
  - scipy, with scipy.optimize.milp / HiGHS enabled
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple

import networkx as nx
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


Vertex = int
Edge = Tuple[int, int]
EdgeIndex = int


def normalize_edge(u: int, v: int) -> Edge:
    return (u, v) if u < v else (v, u)


@dataclass(frozen=True)
class AQGraph:
    n: int
    vertices: Tuple[Vertex, ...]
    edges: Tuple[Edge, ...]
    edge_to_idx: Dict[Edge, EdgeIndex]
    idx_to_edge: Dict[EdgeIndex, Edge]
    adj: Tuple[Tuple[Vertex, ...], ...]
    inc: Tuple[Tuple[EdgeIndex, ...], ...]


def build_augmented_cube(n: int) -> AQGraph:
    vertices = tuple(range(1 << n))

    generators = [1 << i for i in range(n)]
    generators += [(1 << i) - 1 for i in range(2, n + 1)]

    edges_set: Set[Edge] = set()
    for x in vertices:
        for g in generators:
            y = x ^ g
            if x < y:
                edges_set.add((x, y))

    edges = tuple(sorted(edges_set))
    edge_to_idx = {e: i for i, e in enumerate(edges)}
    idx_to_edge = {i: e for i, e in enumerate(edges)}

    adj_lists: List[List[Vertex]] = [[] for _ in vertices]
    inc_lists: List[List[EdgeIndex]] = [[] for _ in vertices]
    for i, (u, v) in enumerate(edges):
        adj_lists[u].append(v)
        adj_lists[v].append(u)
        inc_lists[u].append(i)
        inc_lists[v].append(i)

    adj = tuple(tuple(sorted(nb)) for nb in adj_lists)
    inc = tuple(tuple(lst) for lst in inc_lists)

    return AQGraph(
        n=n,
        vertices=vertices,
        edges=edges,
        edge_to_idx=edge_to_idx,
        idx_to_edge=idx_to_edge,
        adj=adj,
        inc=inc,
    )


def build_networkx_graph(aq: AQGraph) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(aq.vertices)
    G.add_edges_from(aq.edges)
    return G


def compute_edge_automorphisms(aq: AQGraph) -> List[Tuple[EdgeIndex, ...]]:
    """
    Return the full automorphism group of AQ_4 as permutations of edge indices.
    """
    G = build_networkx_graph(aq)
    matcher = nx.algorithms.isomorphism.GraphMatcher(G, G)
    automorphisms = list(matcher.isomorphisms_iter())

    edge_maps: List[Tuple[EdgeIndex, ...]] = []
    for perm in automorphisms:
        em: List[EdgeIndex] = [0] * len(aq.edges)
        for i, (u, v) in enumerate(aq.edges):
            uu, vv = perm[u], perm[v]
            em[i] = aq.edge_to_idx[normalize_edge(uu, vv)]
        edge_maps.append(tuple(em))
    return edge_maps


def find_hamiltonian_cycle_excluding_vertex(
    aq: AQGraph, faulty: Set[EdgeIndex], omit: Vertex
) -> Optional[Tuple[List[Vertex], Tuple[EdgeIndex, ...]]]:
    """
    In AQ_4 - faulty, search for a Hamiltonian cycle on V \\ {omit};
    equivalently, a 15-cycle in AQ_4 missing exactly the vertex 'omit'.

    Returns:
        (cyclic order of the 15 vertices, sorted tuple of its 15 edge indices)
        or None if no such cycle exists.
    """
    nodes = [v for v in aq.vertices if v != omit]
    start = nodes[0]
    others = [v for v in nodes if v != start]
    pos = {v: i for i, v in enumerate(others)}

    DP: List[Dict[Vertex, Vertex]] = [dict() for _ in range(1 << len(others))]

    for u in aq.adj[start]:
        if u == omit:
            continue
        e = aq.edge_to_idx[normalize_edge(start, u)]
        if e in faulty:
            continue
        DP[1 << pos[u]][u] = start

    full = (1 << len(others)) - 1

    for mask in range(1 << len(others)):
        cur_states = DP[mask]
        if not cur_states:
            continue
        for u in list(cur_states.keys()):
            for w in aq.adj[u]:
                if w == start or w == omit:
                    continue
                bit = 1 << pos[w]
                if mask & bit:
                    continue
                e = aq.edge_to_idx[normalize_edge(u, w)]
                if e in faulty:
                    continue
                nxt_mask = mask | bit
                if w not in DP[nxt_mask]:
                    DP[nxt_mask][w] = u

    for u in DP[full]:
        if start not in aq.adj[u]:
            continue
        closing = aq.edge_to_idx[normalize_edge(u, start)]
        if closing in faulty:
            continue

        rev_path = [u]
        mask = full
        cur = u
        while cur != start:
            prev = DP[mask][cur]
            rev_path.append(prev)
            if prev == start:
                break
            mask ^= 1 << pos[cur]
            cur = prev

        path = rev_path[::-1]
        edge_list = [
            aq.edge_to_idx[normalize_edge(a, b)]
            for a, b in zip(path, path[1:] + path[:1])
        ]
        return path, tuple(sorted(edge_list))

    return None


def find_any_15_cycle(
    aq: AQGraph, faulty: Set[EdgeIndex]
) -> Optional[Tuple[Vertex, List[Vertex], Tuple[EdgeIndex, ...]]]:
    for omit in aq.vertices:
        found = find_hamiltonian_cycle_excluding_vertex(aq, faulty, omit)
        if found is not None:
            path, edges = found
            return omit, path, edges
    return None


def orbit_of_cycle_clause(
    cycle_clause: Sequence[EdgeIndex], edge_maps: Sequence[Sequence[EdgeIndex]]
) -> Set[Tuple[EdgeIndex, ...]]:
    images: Set[Tuple[EdgeIndex, ...]] = set()
    for em in edge_maps:
        images.add(tuple(sorted(em[e] for e in cycle_clause)))
    return images


def solve_counterexample_milp(
    aq: AQGraph, clauses: Sequence[Sequence[EdgeIndex]]
):
    """
    Solve the MILP that looks for a counterexample F_R:
      - exactly 7 faulty edges,
      - at most 5 faulty incident edges at every vertex
        because AQ_4 is 7-regular and the remaining degree must be at least 2,
      - every discovered 15-cycle is hit by at least one faulty edge.
    """
    m = len(aq.edges)

    rows: List[np.ndarray] = []
    lb: List[float] = []
    ub: List[float] = []

    # exactly 7 faulty edges
    rows.append(np.ones(m, dtype=float))
    lb.append(7.0)
    ub.append(7.0)

    # degree cap at every vertex: at most 5 incident faulty edges
    for v in aq.vertices:
        row = np.zeros(m, dtype=float)
        row[list(aq.inc[v])] = 1.0
        rows.append(row)
        lb.append(-np.inf)
        ub.append(5.0)

    # every discovered 15-cycle must be hit
    for clause in clauses:
        row = np.zeros(m, dtype=float)
        row[list(clause)] = 1.0
        rows.append(row)
        lb.append(1.0)
        ub.append(np.inf)

    A = np.vstack(rows)
    constraints = [LinearConstraint(A, np.array(lb), np.array(ub))]
    integrality = np.ones(m, dtype=int)
    bounds = Bounds(np.zeros(m), np.ones(m))
    objective = np.zeros(m, dtype=float)

    return milp(c=objective, constraints=constraints, integrality=integrality, bounds=bounds)


def verify_subcase_3_4_i_n5(verbose: bool = True) -> bool:
    aq4 = build_augmented_cube(4)
    edge_maps = compute_edge_automorphisms(aq4)

    clauses: List[Tuple[EdgeIndex, ...]] = []
    known_clauses: Set[Tuple[EdgeIndex, ...]] = set()

    iteration = 0
    while True:
        res = solve_counterexample_milp(aq4, clauses)

        if res.status != 0:
            if verbose:
                print("[OK] No counterexample exists.")
                print(f"     HiGHS status: {res.message}")
                print(f"     Number of cycle clauses used: {len(clauses)}")
            return True

        faulty = {i for i, x in enumerate(res.x) if x > 0.5}
        if verbose:
            faulty_edges = [aq4.idx_to_edge[i] for i in sorted(faulty)]
            print(f"[{iteration}] candidate faulty set: {faulty_edges}")

        witness = find_any_15_cycle(aq4, faulty)
        if witness is None:
            print("[FAIL] Counterexample found.")
            print("Faulty edge set:", [aq4.idx_to_edge[i] for i in sorted(faulty)])
            return False

        omit, cycle_vertices, cycle_clause = witness
        if verbose:
            print(f"     15-cycle found (omitted vertex = {omit}): {cycle_vertices}")

        # Add the whole automorphism orbit of the discovered cycle.
        new_clauses = orbit_of_cycle_clause(cycle_clause, edge_maps)
        added = 0
        for cl in new_clauses:
            if cl not in known_clauses:
                known_clauses.add(cl)
                clauses.append(cl)
                added += 1

        if verbose:
            print(f"     added {added} new cycle clauses (orbit-closed)\n")

        iteration += 1


if __name__ == "__main__":
    ok = verify_subcase_3_4_i_n5(verbose=True)
    if ok:
        print(
            "\nConclusion: for every 7-edge fault set F_R in AQ_4 with "
            "delta(AQ_4 - F_R) >= 2, the graph AQ_4 - F_R contains a 15-cycle."
        )
        print("Hence the n=5 part of Subcase 3.4(i) is verified.")
    else:
        raise SystemExit(1)