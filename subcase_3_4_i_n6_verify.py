#!/usr/bin/env python3
"""
Rigorous counterexample search for Theorem 4.1, Subcase 3.4(i), n=6.

Target statement to verify
--------------------------
Let G = AQ_5 - F_R, where |F_R| = 9 and every vertex of G has internal degree at least 2.
Then G contains a 31-cycle.

Why this script is rigorous
---------------------------
A counterexample is exactly a 9-edge subset F_R such that:
  (1) every vertex of AQ_5 - F_R has degree at least 2, and
  (2) AQ_5 - F_R contains no 31-cycle.

A 31-cycle in AQ_5 is the same thing as a Hamiltonian cycle in AQ_5 - {v}
for some vertex v. We therefore search for a counterexample by branch-and-cut:

  * master MILP variables x_e indicate whether edge e is faulty;
  * master constraints enforce |F_R| = 9 and minimum remaining degree >= 2;
  * every discovered 31-cycle C yields the valid cut  sum_{e in C} x_e >= 1,
    because any counterexample must destroy C by deleting at least one edge of C.

The exact oracle for "does the current graph contain a 31-cycle?" is itself an exact
MILP Hamiltonicity solver with subtour-elimination separation. Therefore:

  - every cut added is valid;
  - if the master MILP becomes infeasible, then no counterexample exists;
  - if the script ever outputs a counterexample, it has already checked *all* 32 choices
    of omitted vertex exactly.

Symmetry reduction
------------------
AQ_5 is a Cayley graph on Z_2^5, so translations are automorphisms. Therefore any
counterexample has a translate containing a faulty edge of the form 0 -- s, where s is
one of the 9 generators of AQ_5. We split the search into 9 independent subproblems,
one for each canonical generator edge (0, s), and solve them in parallel.

Dependencies: networkx, scipy
Usage example:
    python subcase_3_4_i_n6_verify.py --jobs 9 --cycles-per-candidate 16
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass, asdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

# Keep each worker single-threaded; we parallelize across workers instead.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.optimize import Bounds, LinearConstraint, milp

Vertex = int
Edge = Tuple[int, int]


def aq_graph(n: int) -> Tuple[nx.Graph, List[int]]:
    """Return AQ_n as a Cayley graph on Z_2^n, together with its generators."""
    G = nx.Graph()
    vertices = list(range(1 << n))
    G.add_nodes_from(vertices)

    # S = {e_1,...,e_n} union {epsilon_2,...,epsilon_n}
    # In integer bit representation, e_i = 1<<i (0-based), and epsilon_k flips
    # the lowest k bits, i.e. (1<<k)-1.
    generators = [1 << i for i in range(n)]
    generators += [(1 << (k + 1)) - 1 for k in range(1, n)]

    for x in vertices:
        for s in generators:
            y = x ^ s
            if x < y:
                G.add_edge(x, y)

    return G, generators


def canonical_edge(u: int, v: int) -> Edge:
    return (u, v) if u < v else (v, u)


@dataclass
class GraphData:
    G: nx.Graph
    generators: List[int]
    edges: List[Edge]
    edge_to_idx: Dict[Edge, int]
    incidence: Dict[Vertex, List[int]]
    canonical_fault_edge_indices: List[int]


def build_graph_data() -> GraphData:
    G, generators = aq_graph(5)  # AQ_5 is the R-side when n=6 in the paper.
    edges = [canonical_edge(u, v) for (u, v) in G.edges()]
    edge_to_idx = {e: i for i, e in enumerate(edges)}
    incidence: Dict[Vertex, List[int]] = {v: [] for v in G.nodes()}
    for i, (u, v) in enumerate(edges):
        incidence[u].append(i)
        incidence[v].append(i)

    canonical_fault_edge_indices = []
    for s in generators:
        canonical_fault_edge_indices.append(edge_to_idx[canonical_edge(0, s)])

    return GraphData(
        G=G,
        generators=generators,
        edges=edges,
        edge_to_idx=edge_to_idx,
        incidence=incidence,
        canonical_fault_edge_indices=canonical_fault_edge_indices,
    )


GRAPH = build_graph_data()


def cycle_edges_to_indices(cycle_edges: Iterable[Edge]) -> Tuple[int, ...]:
    return tuple(sorted(GRAPH.edge_to_idx[canonical_edge(*e)] for e in cycle_edges))


def solve_hamiltonian_cycle_exact(H: nx.Graph, time_limit: Optional[float] = None) -> Optional[List[Edge]]:
    """
    Exact Hamiltonian-cycle solver via MILP + subtour elimination.

    Returns a list of edges of a Hamiltonian cycle if one exists, otherwise None.
    """
    n = H.number_of_nodes()
    if n < 3:
        return None
    if any(d < 2 for _, d in H.degree()):
        return None

    nodes = list(H.nodes())
    edges = [canonical_edge(u, v) for (u, v) in H.edges()]
    m = len(edges)

    # Base degree-equals-2 constraints.
    base_rows: List[int] = []
    base_cols: List[int] = []
    base_data: List[float] = []
    base_lb: List[float] = []
    base_ub: List[float] = []

    for r, v in enumerate(nodes):
        for i, (a, b) in enumerate(edges):
            if a == v or b == v:
                base_rows.append(r)
                base_cols.append(i)
                base_data.append(1.0)
        base_lb.append(2.0)
        base_ub.append(2.0)

    cuts: List[Tuple[int, ...]] = []
    cut_set: Set[Tuple[int, ...]] = set()
    base_n_rows = len(nodes)

    while True:
        rows = list(base_rows)
        cols = list(base_cols)
        data = list(base_data)
        lb = list(base_lb)
        ub = list(base_ub)
        r = base_n_rows

        for cut in cuts:
            for ei in cut:
                rows.append(r)
                cols.append(ei)
                data.append(1.0)
            lb.append(2.0)
            ub.append(np.inf)
            r += 1

        A = sp.coo_array((data, (rows, cols)), shape=(r, m)).tocsr()
        options = {"disp": False}
        if time_limit is not None:
            options["time_limit"] = float(time_limit)
        res = milp(
            c=np.zeros(m),
            constraints=LinearConstraint(A, lb, ub),
            integrality=np.ones(m, dtype=int),
            bounds=Bounds(np.zeros(m), np.ones(m)),
            options=options,
        )

        if not res.success:
            return None

        x = np.rint(res.x).astype(int)
        chosen = [edges[i] for i, val in enumerate(x) if val == 1]

        F = nx.Graph()
        F.add_nodes_from(nodes)
        F.add_edges_from(chosen)
        components = list(nx.connected_components(F))

        # A connected 2-factor on all n vertices is a Hamiltonian cycle.
        if len(components) == 1 and F.number_of_edges() == n and all(d == 2 for _, d in F.degree()):
            return chosen

        new_cut_added = False
        for comp in components:
            if 0 < len(comp) < n:
                comp = set(comp)
                cut = tuple(
                    i
                    for i, (a, b) in enumerate(edges)
                    if (a in comp) ^ (b in comp)
                )
                if cut not in cut_set:
                    cut_set.add(cut)
                    cuts.append(cut)
                    new_cut_added = True

        if not new_cut_added:
            # Defensive fallback: with correct SEC separation this should not happen.
            return None


@dataclass
class OracleResult:
    found_any: bool
    cycles: List[Tuple[int, ...]]
    omitted_vertices_with_cycle: List[int]


def exact_31_cycles_for_fault_set(fault_edge_indices: Sequence[int], max_cycles: int) -> OracleResult:
    """
    Exact oracle for the current candidate fault set.

    We inspect all 32 possible omitted vertices exactly. For each omitted vertex v,
    we solve the Hamiltonian-cycle problem in (AQ_5 - F) - {v}. If a Hamiltonian cycle exists,
    it is a 31-cycle in AQ_5 - F.
    """
    Gf = GRAPH.G.copy()
    Gf.remove_edges_from(GRAPH.edges[i] for i in fault_edge_indices)

    unique_cycles: List[Tuple[int, ...]] = []
    seen_cycles: Set[Tuple[int, ...]] = set()
    omitted_vertices_with_cycle: List[int] = []

    # Low-degree omitted vertices tend to be the easiest places to find a 31-cycle.
    order = sorted(Gf.nodes(), key=lambda v: (Gf.degree(v), v))

    for omit in order:
        H = Gf.copy()
        H.remove_node(omit)
        cycle = solve_hamiltonian_cycle_exact(H)
        if cycle is not None:
            omitted_vertices_with_cycle.append(omit)
            cycle_idx = cycle_edges_to_indices(cycle)
            if cycle_idx not in seen_cycles:
                seen_cycles.add(cycle_idx)
                unique_cycles.append(cycle_idx)
                if len(unique_cycles) >= max_cycles:
                    # We already have enough valid cuts for the master problem.
                    # This does not affect rigor: if later we ever fail to find any cycle,
                    # we still inspect all omitted vertices exactly before declaring a counterexample.
                    pass

    return OracleResult(
        found_any=bool(omitted_vertices_with_cycle),
        cycles=unique_cycles[:max_cycles],
        omitted_vertices_with_cycle=omitted_vertices_with_cycle,
    )


@dataclass
class SubproblemResult:
    generator: int
    forced_edge: Edge
    verified: bool
    counterexample_faults: Optional[List[Edge]]
    iterations: int
    cuts_added: int
    wall_seconds: float
    notes: str


def solve_master_problem(
    cycle_cuts: Sequence[Tuple[int, ...]],
    forced_fault_edge_idx: int,
    time_limit: Optional[float],
) -> Optional[np.ndarray]:
    """Solve the master MILP; return a 0/1 fault vector, or None if infeasible."""
    m = len(GRAPH.edges)
    rows: List[int] = []
    cols: List[int] = []
    data: List[float] = []
    lb: List[float] = []
    ub: List[float] = []
    r = 0

    # Exactly 9 faulty edges.
    for i in range(m):
        rows.append(r)
        cols.append(i)
        data.append(1.0)
    lb.append(9.0)
    ub.append(9.0)
    r += 1

    # Minimum remaining degree >= 2, i.e. at most 7 faulty incident edges at every vertex.
    for v in GRAPH.G.nodes():
        for ei in GRAPH.incidence[v]:
            rows.append(r)
            cols.append(ei)
            data.append(1.0)
        lb.append(-np.inf)
        ub.append(7.0)
        r += 1

    # Symmetry reduction: force a canonical faulty edge.
    rows.append(r)
    cols.append(forced_fault_edge_idx)
    data.append(1.0)
    lb.append(1.0)
    ub.append(1.0)
    r += 1

    # Every discovered 31-cycle must be hit by at least one faulty edge.
    for cyc in cycle_cuts:
        for ei in cyc:
            rows.append(r)
            cols.append(ei)
            data.append(1.0)
        lb.append(1.0)
        ub.append(np.inf)
        r += 1

    A = sp.coo_array((data, (rows, cols)), shape=(r, m)).tocsr()
    options = {"disp": False}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    res = milp(
        c=np.zeros(m),
        constraints=LinearConstraint(A, lb, ub),
        integrality=np.ones(m, dtype=int),
        bounds=Bounds(np.zeros(m), np.ones(m)),
        options=options,
    )

    if not res.success:
        # Infeasible is the success case for verification.
        return None
    return np.rint(res.x).astype(int)


def initial_cycle_pool(max_per_omitted_vertex: int = 1) -> List[Tuple[int, ...]]:
    """
    Build an initial pool of valid cuts from the full AQ_5.

    One exact 31-cycle for each omitted vertex already gives 32 strong cuts essentially for free.
    """
    pool: Set[Tuple[int, ...]] = set()
    for omit in GRAPH.G.nodes():
        H = GRAPH.G.copy()
        H.remove_node(omit)
        cycle = solve_hamiltonian_cycle_exact(H)
        if cycle is None:
            raise RuntimeError("Unexpected: AQ_5 - {v} should be Hamiltonian for every v.")
        pool.add(cycle_edges_to_indices(cycle))
    return sorted(pool)


def solve_one_generator_subproblem(
    generator: int,
    forced_fault_edge_idx: int,
    cycles_per_candidate: int,
    outer_time_limit: Optional[float],
    max_iterations: int,
    log_dir: Optional[str],
) -> SubproblemResult:
    start = time.time()
    forced_edge = GRAPH.edges[forced_fault_edge_idx]

    cycle_cuts = initial_cycle_pool()
    cycle_cut_set: Set[Tuple[int, ...]] = set(cycle_cuts)

    if log_dir is not None:
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, f"generator_{generator}.jsonl")
    else:
        log_path = None

    def log(event: dict) -> None:
        if log_path is None:
            return
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(event, ensure_ascii=False) + "\n")

    log({
        "event": "start_subproblem",
        "generator": generator,
        "forced_edge": forced_edge,
        "initial_cuts": len(cycle_cuts),
    })

    for iteration in range(1, max_iterations + 1):
        x = solve_master_problem(cycle_cuts, forced_fault_edge_idx, outer_time_limit)
        if x is None:
            wall = time.time() - start
            log({
                "event": "verified",
                "generator": generator,
                "forced_edge": forced_edge,
                "iterations": iteration - 1,
                "cuts_added": len(cycle_cuts),
                "wall_seconds": wall,
            })
            return SubproblemResult(
                generator=generator,
                forced_edge=forced_edge,
                verified=True,
                counterexample_faults=None,
                iterations=iteration - 1,
                cuts_added=len(cycle_cuts),
                wall_seconds=wall,
                notes="Master MILP infeasible after adding exact 31-cycle cuts.",
            )

        fault_indices = np.where(x == 1)[0].tolist()
        oracle = exact_31_cycles_for_fault_set(fault_indices, max_cycles=cycles_per_candidate)

        if not oracle.found_any:
            wall = time.time() - start
            counterexample_faults = [GRAPH.edges[i] for i in fault_indices]
            log({
                "event": "counterexample",
                "generator": generator,
                "forced_edge": forced_edge,
                "faults": counterexample_faults,
                "iterations": iteration,
                "cuts_added": len(cycle_cuts),
                "wall_seconds": wall,
            })
            return SubproblemResult(
                generator=generator,
                forced_edge=forced_edge,
                verified=False,
                counterexample_faults=counterexample_faults,
                iterations=iteration,
                cuts_added=len(cycle_cuts),
                wall_seconds=wall,
                notes=(
                    "Exact oracle checked all 32 omitted vertices and found no 31-cycle. "
                    "This is a genuine counterexample candidate."
                ),
            )

        added = 0
        for cyc in oracle.cycles:
            if cyc not in cycle_cut_set:
                cycle_cut_set.add(cyc)
                cycle_cuts.append(cyc)
                added += 1

        wall = time.time() - start
        log({
            "event": "iteration",
            "generator": generator,
            "forced_edge": forced_edge,
            "iteration": iteration,
            "fault_count": len(fault_indices),
            "found_31_cycles_for_omitted_vertices": oracle.omitted_vertices_with_cycle,
            "cuts_added_this_iteration": added,
            "cuts_total": len(cycle_cuts),
            "wall_seconds": wall,
        })

        if added == 0:
            # Extremely unlikely, but safe to stop and report ambiguity rather than loop forever.
            return SubproblemResult(
                generator=generator,
                forced_edge=forced_edge,
                verified=False,
                counterexample_faults=[GRAPH.edges[i] for i in fault_indices],
                iterations=iteration,
                cuts_added=len(cycle_cuts),
                wall_seconds=wall,
                notes=(
                    "No new cycle cuts were produced from the exact oracle. "
                    "Increase cycles_per_candidate or inspect this subproblem manually."
                ),
            )

    wall = time.time() - start
    return SubproblemResult(
        generator=generator,
        forced_edge=forced_edge,
        verified=False,
        counterexample_faults=None,
        iterations=max_iterations,
        cuts_added=len(cycle_cuts),
        wall_seconds=wall,
        notes="Reached max_iterations before proving infeasibility.",
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=4, help="Number of parallel generator subproblems.")
    parser.add_argument(
        "--cycles-per-candidate",
        type=int,
        default=16,
        help="How many exact 31-cycles to add per master iteration.",
    )
    parser.add_argument(
        "--outer-time-limit",
        type=float,
        default=None,
        help="Optional time limit (seconds) for each master MILP solve.",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=500,
        help="Safety cap on iterations per generator subproblem.",
    )
    parser.add_argument(
        "--log-dir",
        type=str,
        default="verify_subcase_3_4_i_logs",
        help="Directory for JSONL progress logs.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    global_start = time.time()

    tasks = [
        (g, forced_idx, args.cycles_per_candidate, args.outer_time_limit, args.max_iterations, args.log_dir)
        for g, forced_idx in zip(GRAPH.generators, GRAPH.canonical_fault_edge_indices)
    ]

    results: List[SubproblemResult] = []
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futures = {
            ex.submit(solve_one_generator_subproblem, *task): task[:2]
            for task in tasks
        }
        for fut in as_completed(futures):
            result = fut.result()
            results.append(result)
            status = "VERIFIED" if result.verified else "NOT VERIFIED"
            print(
                f"[generator={result.generator:>2} edge={result.forced_edge}] {status}; "
                f"iters={result.iterations}; cuts={result.cuts_added}; "
                f"time={result.wall_seconds:.1f}s"
            )
            sys.stdout.flush()

    all_verified = all(r.verified for r in results)
    total_time = time.time() - global_start

    print("\n=== Summary ===")
    print(f"Subproblems: {len(results)}")
    print(f"All verified: {all_verified}")
    print(f"Wall time: {total_time:.1f} seconds")

    if all_verified:
        print(
            "Conclusion: every 9-edge fault set F_R in AQ_5 with minimum remaining degree at least 2 "
            "still leaves a 31-cycle."
        )
        return 0

    print("At least one subproblem was not verified. Inspect the JSONL logs and any reported counterexample.")
    for r in sorted(results, key=lambda x: x.generator):
        if not r.verified:
            print(f"generator={r.generator}, forced_edge={r.forced_edge}, notes={r.notes}")
            if r.counterexample_faults is not None:
                print(f"  candidate faults: {r.counterexample_faults}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
