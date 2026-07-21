#!/usr/bin/env python3
"""
Computer-assisted verification for Theorem 4.1, Case 3.1.3, n = 6.

Statement verified
------------------
Let F_R be a set of nine faulty edges in AQ_5. If
minimum degree delta(AQ_5 - F_R) >= 2, then AQ_5 - F_R
contains a cycle of length 31.

Verification method
-------------------
A counterexample is a nine-edge set F_R satisfying the degree condition
and meeting every 31-cycle of AQ_5. The script uses an exact
branch-and-cut search:

1. A master MILP selects the faulty edges.
2. An exact Hamiltonian-cycle MILP checks whether the current graph has
   a 31-cycle, equivalently a Hamiltonian cycle after omitting one of
   the 32 vertices.
3. Every discovered 31-cycle gives a valid master cut requiring at
   least one of its edges to be faulty.
4. Translation symmetry reduces the search to nine subproblems, one
   for each generator edge (0, s) of AQ_5.

The script treats only SciPy/HiGHS status 0 (optimal) and status 2
(infeasible) as conclusive. A time limit, iteration limit, unbounded
status, numerical problem, or any other solver termination raises an
exception and is never reported as a successful verification.

Dependencies: networkx, numpy, scipy

Example:
    python case_3_1_3_n6_verify.py --jobs 9 --cycles-per-candidate 16
"""

from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

# Parallelize across symmetry subproblems, not inside BLAS/HiGHS workers.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.optimize import Bounds, LinearConstraint, OptimizeResult, milp

Vertex = int
Edge = Tuple[int, int]
EdgeIndex = int
CycleCut = Tuple[EdgeIndex, ...]

FAULT_COUNT = 9
MIN_REMAINING_DEGREE = 2


class InconclusiveSolverError(RuntimeError):
    """Raised whenever a MILP terminates without a rigorous conclusion."""


def canonical_edge(u: int, v: int) -> Edge:
    return (u, v) if u < v else (v, u)


def aq_graph(n: int) -> Tuple[nx.Graph, List[int]]:
    """Return AQ_n as a Cayley graph on Z_2^n and its generators."""
    graph = nx.Graph()
    vertices = list(range(1 << n))
    graph.add_nodes_from(vertices)

    generators = [1 << i for i in range(n)]
    generators += [(1 << k) - 1 for k in range(2, n + 1)]

    for vertex in vertices:
        for generator in generators:
            graph.add_edge(vertex, vertex ^ generator)
    return graph, generators


@dataclass(frozen=True)
class GraphData:
    graph: nx.Graph
    generators: Tuple[int, ...]
    edges: Tuple[Edge, ...]
    edge_to_idx: Dict[Edge, EdgeIndex]
    incidence: Dict[Vertex, Tuple[EdgeIndex, ...]]
    canonical_fault_edge_indices: Tuple[EdgeIndex, ...]


def build_graph_data() -> GraphData:
    graph, generators = aq_graph(5)
    edges = tuple(sorted({canonical_edge(u, v) for u, v in graph.edges()}))
    edge_to_idx = {edge: idx for idx, edge in enumerate(edges)}

    incidence_lists: Dict[Vertex, List[EdgeIndex]] = {
        vertex: [] for vertex in graph.nodes()
    }
    for idx, (u, v) in enumerate(edges):
        incidence_lists[u].append(idx)
        incidence_lists[v].append(idx)

    canonical_indices = tuple(
        edge_to_idx[canonical_edge(0, generator)] for generator in generators
    )
    return GraphData(
        graph=graph,
        generators=tuple(generators),
        edges=edges,
        edge_to_idx=edge_to_idx,
        incidence={
            vertex: tuple(indices)
            for vertex, indices in incidence_lists.items()
        },
        canonical_fault_edge_indices=canonical_indices,
    )


GRAPH = build_graph_data()


def cycle_edges_to_indices(cycle_edges: Iterable[Edge]) -> CycleCut:
    return tuple(
        sorted(
            GRAPH.edge_to_idx[canonical_edge(u, v)]
            for u, v in cycle_edges
        )
    )


def solver_options(time_limit: Optional[float]) -> Dict[str, float | bool]:
    options: Dict[str, float | bool] = {
        "disp": False,
        "mip_rel_gap": 0.0,
    }
    if time_limit is not None:
        if time_limit <= 0:
            raise ValueError("A MILP time limit must be positive.")
        options["time_limit"] = float(time_limit)
    return options


def require_optimal_or_infeasible(
    result: OptimizeResult,
    problem_name: str,
) -> bool:
    """Return True iff infeasibility is proved; return False iff optimal."""
    if result.status == 2:
        return True
    if result.status == 0 and result.success and result.x is not None:
        return False
    raise InconclusiveSolverError(
        f"{problem_name} did not terminate conclusively: "
        f"status={result.status}, message={result.message}"
    )


def rounded_binary_solution(
    result: OptimizeResult,
    problem_name: str,
) -> np.ndarray:
    """Round and validate an optimal binary MILP solution."""
    if result.status != 0 or not result.success or result.x is None:
        raise RuntimeError(f"Cannot decode a non-optimal {problem_name} result.")
    values = np.asarray(result.x, dtype=float)
    rounded = np.rint(values).astype(int)
    if np.max(np.abs(values - rounded)) > 1e-6:
        raise RuntimeError(f"{problem_name} returned a nonintegral solution.")
    if np.any((rounded < 0) | (rounded > 1)):
        raise RuntimeError(f"{problem_name} returned a nonbinary solution.")
    return rounded


def solve_hamiltonian_cycle_exact(
    graph: nx.Graph,
    time_limit: Optional[float] = None,
) -> Optional[List[Edge]]:
    """Solve Hamiltonicity exactly by MILP with subtour separation."""
    node_count = graph.number_of_nodes()
    if node_count < 3:
        return None
    if any(degree < 2 for _, degree in graph.degree()):
        return None

    nodes = list(graph.nodes())
    edges = sorted({canonical_edge(u, v) for u, v in graph.edges()})
    edge_count = len(edges)
    local_incidence: Dict[Vertex, List[EdgeIndex]] = {
        vertex: [] for vertex in nodes
    }
    for edge_idx, (u, v) in enumerate(edges):
        local_incidence[u].append(edge_idx)
        local_incidence[v].append(edge_idx)

    subtour_cuts: List[Tuple[int, ...]] = []
    known_cuts: Set[Tuple[int, ...]] = set()

    while True:
        rows: List[int] = []
        cols: List[int] = []
        data: List[float] = []
        lower: List[float] = []
        upper: List[float] = []
        row = 0

        # Degree two at every vertex.
        for vertex in nodes:
            for edge_idx in local_incidence[vertex]:
                rows.append(row)
                cols.append(edge_idx)
                data.append(1.0)
            lower.append(2.0)
            upper.append(2.0)
            row += 1

        # For every proper component S found previously, at least two
        # selected edges must cross delta(S).
        for cut in subtour_cuts:
            for edge_idx in cut:
                rows.append(row)
                cols.append(edge_idx)
                data.append(1.0)
            lower.append(2.0)
            upper.append(np.inf)
            row += 1

        matrix = sp.coo_array(
            (data, (rows, cols)), shape=(row, edge_count)
        ).tocsr()
        result = milp(
            c=np.zeros(edge_count),
            constraints=LinearConstraint(matrix, lower, upper),
            integrality=np.ones(edge_count, dtype=int),
            bounds=Bounds(np.zeros(edge_count), np.ones(edge_count)),
            options=solver_options(time_limit),
        )

        if require_optimal_or_infeasible(result, "Hamiltonian-cycle MILP"):
            return None

        solution = rounded_binary_solution(result, "Hamiltonian-cycle MILP")
        chosen_edges = [
            edges[idx] for idx, value in enumerate(solution) if value == 1
        ]

        two_factor = nx.Graph()
        two_factor.add_nodes_from(nodes)
        two_factor.add_edges_from(chosen_edges)
        if any(degree != 2 for _, degree in two_factor.degree()):
            raise RuntimeError("The Hamiltonian-cycle MILP solution is not 2-regular.")

        components = list(nx.connected_components(two_factor))
        if len(components) == 1:
            if len(chosen_edges) != node_count:
                raise RuntimeError("Connected 2-factor has an unexpected edge count.")
            return chosen_edges

        added = 0
        for component in components:
            if len(component) == node_count:
                continue
            component_set = set(component)
            cut = tuple(
                idx
                for idx, (u, v) in enumerate(edges)
                if (u in component_set) ^ (v in component_set)
            )
            if cut not in known_cuts:
                known_cuts.add(cut)
                subtour_cuts.append(cut)
                added += 1

        if added == 0:
            raise RuntimeError(
                "A disconnected 2-factor was found, but no new "
                "subtour-elimination constraint was generated."
            )


@dataclass(frozen=True)
class OracleResult:
    cycles: Tuple[Tuple[EdgeIndex, ...], ...]
    omitted_vertices_with_cycle: Tuple[Vertex, ...]
    checked_all_omitted_vertices: bool

    @property
    def found_any(self) -> bool:
        return bool(self.cycles)


def exact_31_cycle_oracle(
    fault_edge_indices: Sequence[EdgeIndex],
    max_cycles: int,
    time_limit: Optional[float] = None,
) -> OracleResult:
    """Search exactly for fault-free 31-cycles for a candidate fault set.

    The search may stop after max_cycles witnesses are found. If no
    witness is found, all 32 omitted vertices have been checked exactly.
    """
    if max_cycles < 1:
        raise ValueError("max_cycles must be positive.")

    remaining = GRAPH.graph.copy()
    remaining.remove_edges_from(GRAPH.edges[idx] for idx in fault_edge_indices)

    unique_cycles: List[Tuple[EdgeIndex, ...]] = []
    seen_cycles: Set[Tuple[EdgeIndex, ...]] = set()
    omitted_with_cycle: List[Vertex] = []
    omitted_order = sorted(remaining.nodes(), key=lambda v: (remaining.degree(v), v))

    checked_all = True
    for omitted in omitted_order:
        graph = remaining.copy()
        graph.remove_node(omitted)
        cycle = solve_hamiltonian_cycle_exact(graph, time_limit=time_limit)
        if cycle is None:
            continue

        omitted_with_cycle.append(omitted)
        cycle_indices = cycle_edges_to_indices(cycle)
        if cycle_indices not in seen_cycles:
            seen_cycles.add(cycle_indices)
            unique_cycles.append(cycle_indices)

        if len(unique_cycles) >= max_cycles:
            checked_all = False
            break

    return OracleResult(
        cycles=tuple(unique_cycles),
        omitted_vertices_with_cycle=tuple(omitted_with_cycle),
        checked_all_omitted_vertices=checked_all,
    )


def solve_master_problem(
    cycle_cuts: Sequence[Tuple[EdgeIndex, ...]],
    forced_fault_edge_idx: EdgeIndex,
    time_limit: Optional[float],
) -> Optional[np.ndarray]:
    """Return a candidate fault vector, or None if infeasibility is proved."""
    edge_count = len(GRAPH.edges)
    rows: List[int] = []
    cols: List[int] = []
    data: List[float] = []
    lower: List[float] = []
    upper: List[float] = []
    row = 0

    # Exactly nine faulty edges.
    for edge_idx in range(edge_count):
        rows.append(row)
        cols.append(edge_idx)
        data.append(1.0)
    lower.append(float(FAULT_COUNT))
    upper.append(float(FAULT_COUNT))
    row += 1

    # AQ_5 is 9-regular. Minimum remaining degree >= 2 means at most
    # seven incident faulty edges at every vertex.
    max_faulty_incident = len(GRAPH.incidence[0]) - MIN_REMAINING_DEGREE
    for vertex in GRAPH.graph.nodes():
        for edge_idx in GRAPH.incidence[vertex]:
            rows.append(row)
            cols.append(edge_idx)
            data.append(1.0)
        lower.append(-np.inf)
        upper.append(float(max_faulty_incident))
        row += 1

    # Symmetry reduction.
    rows.append(row)
    cols.append(forced_fault_edge_idx)
    data.append(1.0)
    lower.append(1.0)
    upper.append(1.0)
    row += 1

    # Every previously discovered 31-cycle must be hit by a fault.
    for cycle in cycle_cuts:
        for edge_idx in cycle:
            rows.append(row)
            cols.append(edge_idx)
            data.append(1.0)
        lower.append(1.0)
        upper.append(np.inf)
        row += 1

    matrix = sp.coo_array(
        (data, (rows, cols)), shape=(row, edge_count)
    ).tocsr()
    result = milp(
        c=np.zeros(edge_count),
        constraints=LinearConstraint(matrix, lower, upper),
        integrality=np.ones(edge_count, dtype=int),
        bounds=Bounds(np.zeros(edge_count), np.ones(edge_count)),
        options=solver_options(time_limit),
    )

    if require_optimal_or_infeasible(result, "master MILP"):
        return None

    solution = rounded_binary_solution(result, "master MILP")
    fault_indices = np.where(solution == 1)[0]
    if len(fault_indices) != FAULT_COUNT:
        raise RuntimeError("Master MILP candidate has the wrong fault count.")

    max_faulty_incident = len(GRAPH.incidence[0]) - MIN_REMAINING_DEGREE
    for vertex in GRAPH.graph.nodes():
        incident_faults = sum(solution[idx] for idx in GRAPH.incidence[vertex])
        if incident_faults > max_faulty_incident:
            raise RuntimeError("Master MILP candidate violates the degree condition.")
    if solution[forced_fault_edge_idx] != 1:
        raise RuntimeError("Master MILP candidate violates the forced-edge condition.")
    return solution


def initial_cycle_pool(
    oracle_time_limit: Optional[float],
) -> List[Tuple[EdgeIndex, ...]]:
    """Generate one valid 31-cycle cut for every omitted vertex."""
    pool: Set[Tuple[EdgeIndex, ...]] = set()
    for omitted in GRAPH.graph.nodes():
        graph = GRAPH.graph.copy()
        graph.remove_node(omitted)
        cycle = solve_hamiltonian_cycle_exact(graph, time_limit=oracle_time_limit)
        if cycle is None:
            raise RuntimeError(
                f"Unexpected failure: AQ_5 - {{{omitted}}} has no Hamiltonian cycle."
            )
        pool.add(cycle_edges_to_indices(cycle))
    return sorted(pool)


@dataclass(frozen=True)
class SubproblemResult:
    generator: int
    forced_edge: Edge
    verified: bool
    counterexample_faults: Optional[Tuple[Edge, ...]]
    iterations: int
    cycle_cuts: int
    wall_seconds: float
    notes: str


def solve_one_generator_subproblem(
    generator: int,
    forced_fault_edge_idx: EdgeIndex,
    cycles_per_candidate: int,
    master_time_limit: Optional[float],
    oracle_time_limit: Optional[float],
    max_iterations: int,
    log_dir: Optional[str],
) -> SubproblemResult:
    start_time = time.time()
    forced_edge = GRAPH.edges[forced_fault_edge_idx]
    cycle_cuts = initial_cycle_pool(oracle_time_limit)
    known_cuts: Set[Tuple[EdgeIndex, ...]] = set(cycle_cuts)

    log_path: Optional[str] = None
    if log_dir is not None:
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, f"generator_{generator}.jsonl")
        with open(log_path, "w", encoding="utf-8") as handle:
            handle.write("")

    def log_event(event: dict) -> None:
        if log_path is None:
            return
        with open(log_path, "a", encoding="utf-8") as handle:
            handle.write(json.dumps(event, ensure_ascii=False) + "\n")

    log_event(
        {
            "event": "start",
            "generator": generator,
            "forced_edge": forced_edge,
            "initial_cycle_cuts": len(cycle_cuts),
        }
    )

    for iteration in range(1, max_iterations + 1):
        solution = solve_master_problem(
            cycle_cuts,
            forced_fault_edge_idx,
            master_time_limit,
        )

        if solution is None:
            wall_seconds = time.time() - start_time
            result = SubproblemResult(
                generator=generator,
                forced_edge=forced_edge,
                verified=True,
                counterexample_faults=None,
                iterations=iteration - 1,
                cycle_cuts=len(cycle_cuts),
                wall_seconds=wall_seconds,
                notes="The master MILP was proved infeasible.",
            )
            log_event({"event": "verified", **asdict(result)})
            return result

        fault_indices = np.where(solution == 1)[0].tolist()
        oracle = exact_31_cycle_oracle(
            fault_indices,
            max_cycles=cycles_per_candidate,
            time_limit=oracle_time_limit,
        )

        if not oracle.found_any:
            if not oracle.checked_all_omitted_vertices:
                raise RuntimeError(
                    "The oracle reported no cycle without checking all omitted vertices."
                )
            wall_seconds = time.time() - start_time
            counterexample = tuple(GRAPH.edges[idx] for idx in fault_indices)
            result = SubproblemResult(
                generator=generator,
                forced_edge=forced_edge,
                verified=False,
                counterexample_faults=counterexample,
                iterations=iteration,
                cycle_cuts=len(cycle_cuts),
                wall_seconds=wall_seconds,
                notes=(
                    "All 32 omitted vertices were checked exactly and no "
                    "31-cycle was found."
                ),
            )
            log_event({"event": "counterexample", **asdict(result)})
            return result

        added = 0
        for cycle in oracle.cycles:
            if cycle not in known_cuts:
                known_cuts.add(cycle)
                cycle_cuts.append(cycle)
                added += 1

        if added == 0:
            raise RuntimeError(
                "The oracle found a fault-free 31-cycle, but no new master cut "
                "was generated."
            )

        log_event(
            {
                "event": "iteration",
                "generator": generator,
                "forced_edge": forced_edge,
                "iteration": iteration,
                "faults": [GRAPH.edges[idx] for idx in fault_indices],
                "omitted_vertices_with_cycle": oracle.omitted_vertices_with_cycle,
                "checked_all_omitted_vertices": oracle.checked_all_omitted_vertices,
                "cuts_added": added,
                "cycle_cuts_total": len(cycle_cuts),
                "wall_seconds": time.time() - start_time,
            }
        )

    raise RuntimeError(
        f"Generator {generator} reached max_iterations={max_iterations} "
        "without a conclusive result."
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify the n=6 exceptional case in Theorem 4.1, Case 3.1.3."
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(9, os.cpu_count() or 1),
        help="Number of parallel symmetry subproblems (default: min(9, CPU count)).",
    )
    parser.add_argument(
        "--cycles-per-candidate",
        type=int,
        default=16,
        help="Maximum number of 31-cycle cuts obtained from each candidate.",
    )
    parser.add_argument(
        "--master-time-limit",
        type=float,
        default=None,
        help=(
            "Optional HiGHS time limit in seconds for each master MILP. "
            "Reaching the limit is treated as inconclusive and raises an error."
        ),
    )
    parser.add_argument(
        "--oracle-time-limit",
        type=float,
        default=None,
        help=(
            "Optional HiGHS time limit in seconds for each Hamiltonian-cycle MILP. "
            "Reaching the limit is treated as inconclusive and raises an error."
        ),
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=100000,
        help="Safety cap for each branch-and-cut subproblem.",
    )
    parser.add_argument(
        "--log-dir",
        default="verify_case_3_1_3_logs",
        help="Directory for JSONL progress logs; use an empty string to disable.",
    )
    parser.add_argument(
        "--summary-json",
        default="case_3_1_3_n6_summary.json",
        help="Path for the final machine-readable summary.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise ValueError("--jobs must be positive.")
    if args.cycles_per_candidate < 1:
        raise ValueError("--cycles-per-candidate must be positive.")
    if args.max_iterations < 1:
        raise ValueError("--max-iterations must be positive.")

    log_dir = args.log_dir or None
    tasks = list(zip(GRAPH.generators, GRAPH.canonical_fault_edge_indices))
    results: List[SubproblemResult] = []

    print("Verifying nine translation-symmetry subproblems for AQ_5...")
    print(f"Parallel jobs: {args.jobs}")

    try:
        if args.jobs == 1:
            for generator, forced_idx in tasks:
                result = solve_one_generator_subproblem(
                    generator,
                    forced_idx,
                    args.cycles_per_candidate,
                    args.master_time_limit,
                    args.oracle_time_limit,
                    args.max_iterations,
                    log_dir,
                )
                results.append(result)
                label = "VERIFIED" if result.verified else "COUNTEREXAMPLE"
                print(
                    f"generator={generator:2d}, forced_edge={result.forced_edge}: "
                    f"{label}, iterations={result.iterations}, "
                    f"cuts={result.cycle_cuts}, time={result.wall_seconds:.2f}s"
                )
        else:
            with ProcessPoolExecutor(max_workers=min(args.jobs, len(tasks))) as executor:
                future_to_generator = {
                    executor.submit(
                        solve_one_generator_subproblem,
                        generator,
                        forced_idx,
                        args.cycles_per_candidate,
                        args.master_time_limit,
                        args.oracle_time_limit,
                        args.max_iterations,
                        log_dir,
                    ): generator
                    for generator, forced_idx in tasks
                }
                for future in as_completed(future_to_generator):
                    result = future.result()
                    results.append(result)
                    label = "VERIFIED" if result.verified else "COUNTEREXAMPLE"
                    print(
                        f"generator={result.generator:2d}, "
                        f"forced_edge={result.forced_edge}: {label}, "
                        f"iterations={result.iterations}, cuts={result.cycle_cuts}, "
                        f"time={result.wall_seconds:.2f}s"
                    )
    except Exception as exc:
        print(f"[INCONCLUSIVE] Verification stopped: {exc}")
        return 2

    results.sort(key=lambda item: item.generator)
    summary = {
        "theorem": "Theorem 4.1",
        "case": "Case 3.1.3",
        "n": 6,
        "statement": (
            "Every nine-edge fault set F_R in AQ_5 with "
            "delta(AQ_5-F_R)>=2 leaves a 31-cycle."
        ),
        "verified": all(result.verified for result in results) and len(results) == 9,
        "results": [asdict(result) for result in results],
    }
    with open(args.summary_json, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    counterexamples = [result for result in results if not result.verified]
    if counterexamples:
        print("\n[COUNTEREXAMPLE] At least one symmetry subproblem failed.")
        for result in counterexamples:
            print(
                f"generator={result.generator}, faults={result.counterexample_faults}"
            )
        return 1

    if len(results) != 9:
        print("[INCONCLUSIVE] Not all nine subproblems were completed.")
        return 2

    print(
        "\n[VERIFIED] All nine symmetry subproblems are infeasible after "
        "adding exact 31-cycle cuts."
    )
    print(
        "Conclusion: every nine-edge fault set F_R in AQ_5 with "
        "delta(AQ_5 - F_R) >= 2 leaves a 31-cycle."
    )
    print("Hence the n=6 part of Case 3.1.3 is verified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
