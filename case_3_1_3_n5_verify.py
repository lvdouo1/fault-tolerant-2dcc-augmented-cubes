#!/usr/bin/env python3
"""
Computer-assisted verification for Theorem 4.1, Case 3.1.3, n = 5.

Statement verified
------------------
For every set F_R of seven faulty edges in AQ_4 such that
minimum degree delta(AQ_4 - F_R) >= 2, the graph AQ_4 - F_R
contains a cycle of length 15.

Method
------
The script searches for a counterexample by constraint generation.
The master MILP selects seven faulty edges, enforces minimum remaining
degree at least two, and requires every previously discovered 15-cycle
to contain a faulty edge. For each candidate fault set, an exact
Held--Karp dynamic program searches for a 15-cycle. If one is found,
the corresponding cycle cut and all its automorphic images are added.
If the master MILP is proved infeasible, no counterexample exists.

A solver run is accepted as conclusive only when SciPy/HiGHS returns
status 0 (optimal) or status 2 (infeasible). Time limits and all other
termination statuses raise an exception and are never interpreted as
verification success.

Dependencies: networkx, numpy, scipy
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.optimize import Bounds, LinearConstraint, OptimizeResult, milp

Vertex = int
Edge = Tuple[int, int]
EdgeIndex = int

FAULT_COUNT = 7
MIN_REMAINING_DEGREE = 2


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
    """Build AQ_n as Cay(Z_2^n, {e_i} union {epsilon_i : 2 <= i <= n})."""
    vertices = tuple(range(1 << n))
    generators = [1 << i for i in range(n)]
    generators += [(1 << i) - 1 for i in range(2, n + 1)]

    edge_set: Set[Edge] = set()
    for x in vertices:
        for generator in generators:
            edge_set.add(normalize_edge(x, x ^ generator))

    edges = tuple(sorted(edge_set))
    edge_to_idx = {edge: idx for idx, edge in enumerate(edges)}
    idx_to_edge = {idx: edge for idx, edge in enumerate(edges)}

    adj_lists: List[List[Vertex]] = [[] for _ in vertices]
    inc_lists: List[List[EdgeIndex]] = [[] for _ in vertices]
    for idx, (u, v) in enumerate(edges):
        adj_lists[u].append(v)
        adj_lists[v].append(u)
        inc_lists[u].append(idx)
        inc_lists[v].append(idx)

    return AQGraph(
        n=n,
        vertices=vertices,
        edges=edges,
        edge_to_idx=edge_to_idx,
        idx_to_edge=idx_to_edge,
        adj=tuple(tuple(sorted(neighbors)) for neighbors in adj_lists),
        inc=tuple(tuple(indices) for indices in inc_lists),
    )


def build_networkx_graph(aq: AQGraph) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(aq.vertices)
    graph.add_edges_from(aq.edges)
    return graph


def compute_edge_automorphisms(aq: AQGraph) -> List[Tuple[EdgeIndex, ...]]:
    """Return the full automorphism group as permutations of edge indices."""
    graph = build_networkx_graph(aq)
    matcher = nx.algorithms.isomorphism.GraphMatcher(graph, graph)

    edge_maps: List[Tuple[EdgeIndex, ...]] = []
    for vertex_map in matcher.isomorphisms_iter():
        image = [0] * len(aq.edges)
        for idx, (u, v) in enumerate(aq.edges):
            mapped_edge = normalize_edge(vertex_map[u], vertex_map[v])
            image[idx] = aq.edge_to_idx[mapped_edge]
        edge_maps.append(tuple(image))
    return edge_maps


def find_hamiltonian_cycle_excluding_vertex(
    aq: AQGraph,
    faulty: Set[EdgeIndex],
    omitted: Vertex,
) -> Optional[Tuple[List[Vertex], Tuple[EdgeIndex, ...]]]:
    """Find a Hamiltonian cycle of AQ_4 - F - {omitted} by exact DP."""
    nodes = [vertex for vertex in aq.vertices if vertex != omitted]
    start = nodes[0]
    other_nodes = [vertex for vertex in nodes if vertex != start]
    bit_position = {vertex: idx for idx, vertex in enumerate(other_nodes)}

    # predecessor[mask][u] is the predecessor of u on one path from start
    # visiting exactly the vertices represented by mask among other_nodes.
    predecessor: List[Dict[Vertex, Vertex]] = [
        {} for _ in range(1 << len(other_nodes))
    ]

    for neighbor in aq.adj[start]:
        if neighbor == omitted:
            continue
        edge_idx = aq.edge_to_idx[normalize_edge(start, neighbor)]
        if edge_idx not in faulty:
            predecessor[1 << bit_position[neighbor]][neighbor] = start

    full_mask = (1 << len(other_nodes)) - 1
    for mask in range(1 << len(other_nodes)):
        if not predecessor[mask]:
            continue
        for endpoint in tuple(predecessor[mask]):
            for neighbor in aq.adj[endpoint]:
                if neighbor in (start, omitted):
                    continue
                bit = 1 << bit_position[neighbor]
                if mask & bit:
                    continue
                edge_idx = aq.edge_to_idx[normalize_edge(endpoint, neighbor)]
                if edge_idx in faulty:
                    continue
                next_mask = mask | bit
                predecessor[next_mask].setdefault(neighbor, endpoint)

    for endpoint in predecessor[full_mask]:
        if start not in aq.adj[endpoint]:
            continue
        closing_idx = aq.edge_to_idx[normalize_edge(endpoint, start)]
        if closing_idx in faulty:
            continue

        reverse_path = [endpoint]
        mask = full_mask
        current = endpoint
        while current != start:
            previous = predecessor[mask][current]
            reverse_path.append(previous)
            if previous == start:
                break
            mask ^= 1 << bit_position[current]
            current = previous

        cycle_vertices = reverse_path[::-1]
        cycle_edges = tuple(
            sorted(
                aq.edge_to_idx[normalize_edge(u, v)]
                for u, v in zip(
                    cycle_vertices,
                    cycle_vertices[1:] + cycle_vertices[:1],
                )
            )
        )
        return cycle_vertices, cycle_edges

    return None


def find_any_15_cycle(
    aq: AQGraph,
    faulty: Set[EdgeIndex],
) -> Optional[Tuple[Vertex, List[Vertex], Tuple[EdgeIndex, ...]]]:
    """Return one 15-cycle, or None after checking all omitted vertices."""
    for omitted in aq.vertices:
        result = find_hamiltonian_cycle_excluding_vertex(aq, faulty, omitted)
        if result is not None:
            vertices, edge_indices = result
            return omitted, vertices, edge_indices
    return None


def orbit_of_cycle_clause(
    cycle_clause: Sequence[EdgeIndex],
    edge_maps: Sequence[Sequence[EdgeIndex]],
) -> Set[Tuple[EdgeIndex, ...]]:
    return {
        tuple(sorted(edge_map[idx] for idx in cycle_clause))
        for edge_map in edge_maps
    }


def solve_counterexample_milp(
    aq: AQGraph,
    cycle_clauses: Sequence[Sequence[EdgeIndex]],
) -> OptimizeResult:
    """Solve the master MILP for a candidate counterexample."""
    edge_count = len(aq.edges)
    rows: List[int] = []
    cols: List[int] = []
    data: List[float] = []
    lower: List[float] = []
    upper: List[float] = []
    row = 0

    # Exactly FAULT_COUNT faulty edges.
    for edge_idx in range(edge_count):
        rows.append(row)
        cols.append(edge_idx)
        data.append(1.0)
    lower.append(float(FAULT_COUNT))
    upper.append(float(FAULT_COUNT))
    row += 1

    # AQ_4 is 7-regular. Minimum remaining degree >= 2 means at most
    # five incident faulty edges at every vertex.
    max_faulty_incident = len(aq.inc[0]) - MIN_REMAINING_DEGREE
    for vertex in aq.vertices:
        for edge_idx in aq.inc[vertex]:
            rows.append(row)
            cols.append(edge_idx)
            data.append(1.0)
        lower.append(-np.inf)
        upper.append(float(max_faulty_incident))
        row += 1

    # Every previously discovered 15-cycle must be hit by a fault.
    for clause in cycle_clauses:
        for edge_idx in clause:
            rows.append(row)
            cols.append(edge_idx)
            data.append(1.0)
        lower.append(1.0)
        upper.append(np.inf)
        row += 1

    matrix = sp.coo_array(
        (data, (rows, cols)), shape=(row, edge_count)
    ).tocsr()
    return milp(
        c=np.zeros(edge_count),
        constraints=LinearConstraint(matrix, lower, upper),
        integrality=np.ones(edge_count, dtype=int),
        bounds=Bounds(np.zeros(edge_count), np.ones(edge_count)),
        options={"disp": False, "mip_rel_gap": 0.0},
    )


def decode_fault_set(aq: AQGraph, result: OptimizeResult) -> Set[EdgeIndex]:
    """Validate and decode an optimal binary master-MILP solution."""
    if result.status != 0 or not result.success or result.x is None:
        raise RuntimeError("Attempted to decode a non-optimal MILP result.")

    rounded = np.rint(np.asarray(result.x)).astype(int)
    if np.max(np.abs(np.asarray(result.x) - rounded)) > 1e-6:
        raise RuntimeError("MILP returned a nonintegral candidate solution.")

    faulty = {idx for idx, value in enumerate(rounded) if value == 1}
    if len(faulty) != FAULT_COUNT:
        raise RuntimeError("Decoded candidate has the wrong number of faults.")

    max_faulty_incident = len(aq.inc[0]) - MIN_REMAINING_DEGREE
    for vertex in aq.vertices:
        if sum(idx in faulty for idx in aq.inc[vertex]) > max_faulty_incident:
            raise RuntimeError("Decoded candidate violates the degree constraint.")
    return faulty


def verify_case_3_1_3_n5(verbose: bool = True) -> bool:
    aq4 = build_augmented_cube(4)
    edge_maps = compute_edge_automorphisms(aq4)

    clauses: List[Tuple[EdgeIndex, ...]] = []
    known_clauses: Set[Tuple[EdgeIndex, ...]] = set()
    iteration = 0

    while True:
        result = solve_counterexample_milp(aq4, clauses)

        if result.status == 2:
            if verbose:
                print("[VERIFIED] The counterexample MILP is infeasible.")
                print(f"HiGHS message: {result.message}")
                print(f"Iterations: {iteration}")
                print(f"Cycle cuts: {len(clauses)}")
            return True

        if result.status != 0 or not result.success or result.x is None:
            raise RuntimeError(
                "The master MILP did not terminate conclusively: "
                f"status={result.status}, message={result.message}"
            )

        faulty = decode_fault_set(aq4, result)
        if verbose:
            print(
                f"[{iteration}] candidate faults: "
                f"{[aq4.idx_to_edge[idx] for idx in sorted(faulty)]}"
            )

        witness = find_any_15_cycle(aq4, faulty)
        if witness is None:
            print("[COUNTEREXAMPLE] No 15-cycle exists for this fault set.")
            print([aq4.idx_to_edge[idx] for idx in sorted(faulty)])
            return False

        omitted, cycle_vertices, cycle_clause = witness
        new_clauses = orbit_of_cycle_clause(cycle_clause, edge_maps)
        added = 0
        for clause in new_clauses:
            if clause not in known_clauses:
                known_clauses.add(clause)
                clauses.append(clause)
                added += 1

        if added == 0:
            raise RuntimeError(
                "A fault-free 15-cycle was found, but no new cycle cut was generated."
            )

        if verbose:
            print(f"  omitted vertex: {omitted}")
            print(f"  15-cycle: {cycle_vertices}")
            print(f"  added {added} orbit-closed cycle cuts\n")
        iteration += 1


if __name__ == "__main__":
    verified = verify_case_3_1_3_n5(verbose=True)
    if not verified:
        raise SystemExit(1)
    print(
        "\nConclusion: for every seven-edge fault set F_R in AQ_4 with "
        "delta(AQ_4 - F_R) >= 2, AQ_4 - F_R contains a 15-cycle."
    )
    print("Hence the n=5 part of Case 3.1.3 is verified.")
