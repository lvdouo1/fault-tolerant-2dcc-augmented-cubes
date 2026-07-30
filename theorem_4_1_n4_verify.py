#!/usr/bin/env python3
"""
Computer verification for Theorem 4.1 in the base case n = 4.

Theorem 4.1 (n = 4 case):
For every fault set F ⊆ E(AQ_4) with |F| <= 5, and for every integer l with 3 <= l <= 8,
there exist two vertex-disjoint cycles C1 and C2 in AQ_4 - F such that
    |C1| = l,
    |C2| = 16 - l,
and V(C1) ∪ V(C2) = V(AQ_4).

This script verifies that statement by exhaustive search, reduced by the automorphism
group of AQ_4.

Method overview
---------------
1. Build AQ_4 by its standard Cayley description on GF(2)^4.
2. Compute the affine automorphism group of AQ_4:
       x |-> A x + b
   where A preserves the generator set of AQ_4.
3. Enumerate orbit representatives of edge-fault sets F with |F| <= 5 under that
   automorphism group.
4. For each representative graph G = AQ_4 - F, compute subset-DP information:
       cycle_exists[mask] = whether G[mask] has a Hamiltonian cycle.
5. For each l = 3,4,...,8, search a vertex subset S of size l with
       cycle_exists[S] = True and cycle_exists[V \\ S] = True.
   Such a pair gives the desired two disjoint cycles covering all vertices.

Because the property is invariant under graph automorphisms, checking one
representative from each orbit of fault sets is sufficient.

Notes
-----
- This script is intended as a verification/certification tool for the paper.
- It prints the number of orbit representatives at each fault size and confirms
  whether every representative passes.
- Optional witness reconstruction is included.
- Exit status 0 means verified, 1 means a counterexample was found, and
  2 means the computation stopped without a conclusive result.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime
from itertools import combinations, product
import argparse
import json
from pathlib import Path
import time
from typing import Dict, Iterable, List, Optional, Sequence, TextIO, Tuple

N = 16
ALL_MASK = (1 << N) - 1


def bit(v: int) -> int:
    return 1 << v


def bits(mask: int) -> Iterable[int]:
    while mask:
        lb = mask & -mask
        yield lb.bit_length() - 1
        mask ^= lb


def fmt_vertex(v: int) -> str:
    return format(v, "04b")


def fmt_edge(edge: Tuple[int, int]) -> str:
    u, v = sorted(edge)
    return f"({fmt_vertex(u)}, {fmt_vertex(v)})"


def fmt_edge_set(edge_ids: Sequence[int], edge_list: Sequence[Tuple[int, int]]) -> str:
    return "{" + ", ".join(fmt_edge(edge_list[i]) for i in edge_ids) + "}"


def fmt_cycle(vertices: Sequence[int]) -> str:
    return "[" + ", ".join(fmt_vertex(v) for v in vertices) + "]"


# ---------------------------------------------------------------------------
# AQ_4 construction
# ---------------------------------------------------------------------------

def build_aq4() -> Tuple[List[int], List[Tuple[int, int]], Dict[Tuple[int, int], int]]:
    """
    Build AQ_4 as a Cayley graph on GF(2)^4.

    Generators:
        e1, e2, e3, e4, eps2, eps3, eps4
    where eps_i = e1 + ... + ei over GF(2).

    Vertex encoding:
        4-bit integer x in {0,...,15}.
    """
    adjacency = [0] * N

    for x in range(N):
        # e1, e2, e3, e4
        for i in range(4):
            adjacency[x] |= bit(x ^ (1 << i))
        # eps2, eps3, eps4
        for i in range(2, 5):
            adjacency[x] |= bit(x ^ ((1 << i) - 1))

    edges: List[Tuple[int, int]] = []
    edge_index: Dict[Tuple[int, int], int] = {}
    for u in range(N):
        for v in bits(adjacency[u]):
            if u < v:
                edge_index[(u, v)] = len(edges)
                edges.append((u, v))

    return adjacency, edges, edge_index


BASE_ADJ, EDGE_LIST, EDGE_INDEX = build_aq4()
M = len(EDGE_LIST)


def validate_base_graph() -> None:
    """Validate the basic structural data used throughout the verification."""
    if N != 16:
        raise RuntimeError(f"AQ_4 must have 16 vertices, found N={N}.")
    if len(BASE_ADJ) != N:
        raise RuntimeError("The adjacency table has the wrong size.")
    if M != 56:
        raise RuntimeError(f"AQ_4 must have 56 edges, found {M}.")
    if len(EDGE_INDEX) != M:
        raise RuntimeError("The edge-index map is inconsistent with the edge list.")

    for vertex, adjacency in enumerate(BASE_ADJ):
        if adjacency & bit(vertex):
            raise RuntimeError(f"A self-loop was found at vertex {vertex}.")
        if adjacency.bit_count() != 7:
            raise RuntimeError(
                f"AQ_4 must be 7-regular, but vertex {vertex} has "
                f"degree {adjacency.bit_count()}."
            )
        for neighbor in bits(adjacency):
            if not (BASE_ADJ[neighbor] & bit(vertex)):
                raise RuntimeError(
                    f"Adjacency is not symmetric for vertices {vertex}, {neighbor}."
                )

    for edge_id, (u, v) in enumerate(EDGE_LIST):
        if not (0 <= u < v < N):
            raise RuntimeError(f"Invalid canonical edge {EDGE_LIST[edge_id]}.")
        if EDGE_INDEX.get((u, v)) != edge_id:
            raise RuntimeError("The edge list and edge-index map disagree.")
        if not (BASE_ADJ[u] & bit(v)):
            raise RuntimeError(f"Listed edge {(u, v)} is absent from adjacency.")


validate_base_graph()


# ---------------------------------------------------------------------------
# Linear algebra over GF(2)
# ---------------------------------------------------------------------------

def gf2_mat_vec_mul(rows: Sequence[int], x: int) -> int:
    """
    rows is a list of 4 row bitmasks, each in {0,...,15}.
    Return A x over GF(2) as a 4-bit integer.
    """
    out = 0
    for i, row in enumerate(rows):
        parity = (row & x).bit_count() & 1
        out |= parity << i
    return out


def gf2_rank(rows: Sequence[int]) -> int:
    rows = [r for r in rows if r]
    rank = 0
    col = 3
    rows = rows[:]
    while col >= 0 and rows:
        pivot = None
        for i in range(rank, len(rows)):
            if (rows[i] >> col) & 1:
                pivot = i
                break
        if pivot is None:
            col -= 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and ((rows[i] >> col) & 1):
                rows[i] ^= rows[rank]
        rank += 1
        col -= 1
    return rank


def is_invertible(rows: Sequence[int]) -> bool:
    return gf2_rank(rows) == 4


# ---------------------------------------------------------------------------
# Automorphism group of AQ_4
# ---------------------------------------------------------------------------

def generator_set() -> List[int]:
    gens = []
    for i in range(4):
        gens.append(1 << i)          # e_i
    for i in range(2, 5):
        gens.append((1 << i) - 1)    # eps_i
    return gens


GENS = sorted(generator_set())
GENS_SET = set(GENS)


def compute_linear_stabilizer() -> List[Tuple[int, int, int, int]]:
    """
    Return all invertible 4x4 GF(2)-matrices A such that A(GENS) = GENS.
    Matrices are represented by row bitmasks.
    """
    mats: List[Tuple[int, int, int, int]] = []
    for rows in product(range(16), repeat=4):
        if not is_invertible(rows):
            continue
        image = {gf2_mat_vec_mul(rows, g) for g in GENS}
        if image == GENS_SET:
            mats.append(tuple(rows))
    return mats


def compute_automorphisms() -> List[List[int]]:
    """
    Compute the full affine automorphism group:
        x |-> A x + b
    where A preserves GENS and b is arbitrary.
    Return each automorphism as a length-16 list perm with perm[x] = image of x.
    """
    linear_parts = compute_linear_stabilizer()
    autos: List[List[int]] = []
    for A in linear_parts:
        for b in range(16):
            perm = [gf2_mat_vec_mul(A, x) ^ b for x in range(16)]
            autos.append(perm)

    # Deduplicate just in case
    uniq: Dict[Tuple[int, ...], List[int]] = {}
    for p in autos:
        uniq.setdefault(tuple(p), p)
    return list(uniq.values())


def edge_image_under_perm(edge_id: int, perm: Sequence[int]) -> int:
    u, v = EDGE_LIST[edge_id]
    a, b = perm[u], perm[v]
    if a > b:
        a, b = b, a
    return EDGE_INDEX[(a, b)]


def compute_edge_action_perms(vertex_perms: Sequence[Sequence[int]]) -> List[List[int]]:
    """
    Convert each vertex permutation into the induced permutation on edges.
    """
    out: List[List[int]] = []
    for vp in vertex_perms:
        ep = [0] * M
        for eid in range(M):
            ep[eid] = edge_image_under_perm(eid, vp)
        out.append(ep)
    return out


def validate_automorphism_actions(
    vertex_perms: Sequence[Sequence[int]],
    edge_perms: Sequence[Sequence[int]],
) -> None:
    """Check that all computed actions are genuine vertex/edge permutations."""
    if not vertex_perms:
        raise RuntimeError("No affine automorphisms were found.")
    if len(vertex_perms) != len(edge_perms):
        raise RuntimeError("Vertex and edge automorphism lists have different sizes.")

    expected_vertices = set(range(N))
    expected_edges = set(range(M))
    for automorphism_id, (vp, ep) in enumerate(zip(vertex_perms, edge_perms)):
        if len(vp) != N or set(vp) != expected_vertices:
            raise RuntimeError(
                f"Automorphism {automorphism_id} is not a permutation of vertices."
            )
        if len(ep) != M or set(ep) != expected_edges:
            raise RuntimeError(
                f"Automorphism {automorphism_id} is not a permutation of edges."
            )

        for edge_id, (u, v) in enumerate(EDGE_LIST):
            image = tuple(sorted((vp[u], vp[v])))
            if EDGE_LIST[ep[edge_id]] != image:
                raise RuntimeError(
                    f"Automorphism {automorphism_id} has an inconsistent edge image."
                )


# ---------------------------------------------------------------------------
# Canonical representatives of fault sets under automorphisms
# ---------------------------------------------------------------------------

def canonical_fault_set(edge_ids: Sequence[int], edge_autos: Sequence[Sequence[int]]) -> Tuple[int, ...]:
    s = tuple(sorted(edge_ids))
    best = None
    for ep in edge_autos:
        img = tuple(sorted(ep[e] for e in s))
        if best is None or img < best:
            best = img
    assert best is not None
    return best


def enumerate_fault_orbit_reps(max_faults: int, edge_autos: Sequence[Sequence[int]]) -> List[List[Tuple[int, ...]]]:
    """
    Enumerate orbit representatives of k-edge subsets for k = 0..max_faults.

    We build representatives incrementally:
      reps[k+1] = canonical forms of reps[k] with one added edge.
    """
    reps: List[List[Tuple[int, ...]]] = [[] for _ in range(max_faults + 1)]
    reps[0] = [tuple()]

    for k in range(max_faults):
        seen = set()
        nxt: List[Tuple[int, ...]] = []
        for rep in reps[k]:
            rep_set = set(rep)
            for e in range(M):
                if e in rep_set:
                    continue
                cand = tuple(sorted(rep + (e,)))
                canon = canonical_fault_set(cand, edge_autos)
                if canon not in seen:
                    seen.add(canon)
                    nxt.append(canon)
        reps[k + 1] = sorted(nxt)

    return reps


# ---------------------------------------------------------------------------
# DP for Hamiltonian paths / cycles on induced subgraphs
# ---------------------------------------------------------------------------

@dataclass
class InstanceDP:
    neighbors: List[int]
    dp: List[List[int]]
    cycle_exists: List[bool]


def precompute_instance(faulty_edge_ids: Sequence[int]) -> InstanceDP:
    """
    For G = AQ_4 - F, compute:
      dp[mask][v] : bitset of all starts u such that G[mask] has a Hamiltonian u-v path
      cycle_exists[mask] : whether G[mask] has a Hamiltonian cycle
    """
    neighbors = BASE_ADJ[:]
    for eid in faulty_edge_ids:
        a, b = EDGE_LIST[eid]
        neighbors[a] &= ~bit(b)
        neighbors[b] &= ~bit(a)

    size = 1 << N
    dp = [[0] * N for _ in range(size)]

    for v in range(N):
        dp[bit(v)][v] = bit(v)

    for mask in range(1, size):
        if mask & (mask - 1) == 0:
            continue
        m = mask
        while m:
            lb = m & -m
            v = lb.bit_length() - 1
            prev = mask ^ lb
            starts = 0
            for w in bits(neighbors[v] & prev):
                starts |= dp[prev][w]
            dp[mask][v] = starts
            m ^= lb

    cycle_exists = [False] * size
    for mask in range(size):
        if mask.bit_count() < 3:
            continue
        ok = False
        for v in bits(mask):
            if dp[mask][v] & neighbors[v] & mask:
                ok = True
                break
        cycle_exists[mask] = ok

    return InstanceDP(neighbors=neighbors, dp=dp, cycle_exists=cycle_exists)


def reconstruct_path(instance: InstanceDP, mask: int, start: int, end: int) -> List[int]:
    assert (instance.dp[mask][end] >> start) & 1
    if mask == bit(start) == bit(end):
        return [start]

    prev = mask ^ bit(end)
    for w in bits(instance.neighbors[end] & prev):
        if (instance.dp[prev][w] >> start) & 1:
            path = reconstruct_path(instance, prev, start, w)
            path.append(end)
            return path

    raise RuntimeError("Path reconstruction failed.")


def reconstruct_cycle(instance: InstanceDP, mask: int) -> List[int]:
    assert instance.cycle_exists[mask]
    for end in bits(mask):
        starts = instance.dp[mask][end] & instance.neighbors[end] & mask
        if starts:
            start = (starts & -starts).bit_length() - 1
            return reconstruct_path(instance, mask, start, end)
    raise RuntimeError("Cycle reconstruction failed.")


def find_cycle_cover_witness(instance: InstanceDP, l: int) -> Optional[Tuple[List[int], List[int]]]:
    """
    Find cycles C1, C2 with |C1| = l and |C2| = 16-l covering all vertices.
    Return one witness as two cyclic vertex lists, or None if no witness exists.
    """
    for subset in combinations(range(N), l):
        mask = 0
        for v in subset:
            mask |= bit(v)
        comp = ALL_MASK ^ mask
        if instance.cycle_exists[mask] and instance.cycle_exists[comp]:
            c1 = reconstruct_cycle(instance, mask)
            c2 = reconstruct_cycle(instance, comp)
            return c1, c2
    return None


def validate_cycle(instance: InstanceDP, cycle: Sequence[int]) -> None:
    """Independently check a reconstructed cycle against the remaining graph."""
    if len(cycle) < 3:
        raise RuntimeError("A reconstructed cycle has length less than three.")
    if len(set(cycle)) != len(cycle):
        raise RuntimeError("A reconstructed cycle repeats a vertex.")
    if any(vertex < 0 or vertex >= N for vertex in cycle):
        raise RuntimeError("A reconstructed cycle contains an invalid vertex.")

    for u, v in zip(cycle, cycle[1:] + cycle[:1]):
        if not (instance.neighbors[u] & bit(v)):
            raise RuntimeError(
                f"A reconstructed cycle uses the absent edge ({u}, {v})."
            )


def validate_cycle_cover_witness(
    instance: InstanceDP,
    l: int,
    witness: Tuple[List[int], List[int]],
) -> None:
    """Check lengths, disjointness, coverage, and all cycle edges."""
    c1, c2 = witness
    if len(c1) != l or len(c2) != N - l:
        raise RuntimeError(
            f"A witness has lengths {len(c1)} and {len(c2)}, expected "
            f"{l} and {N - l}."
        )
    validate_cycle(instance, c1)
    validate_cycle(instance, c2)

    vertices1 = set(c1)
    vertices2 = set(c2)
    if vertices1 & vertices2:
        raise RuntimeError("The reconstructed cycles are not vertex-disjoint.")
    if vertices1 | vertices2 != set(range(N)):
        raise RuntimeError("The reconstructed cycles do not cover all vertices.")


# ---------------------------------------------------------------------------
# Reporting and main verification
# ---------------------------------------------------------------------------


class VerificationReporter:
    """Write detailed records to a log while keeping console output compact."""

    def __init__(
        self,
        log_path: Path,
        progress_every: int = 100,
    ) -> None:
        if progress_every < 1:
            raise ValueError("progress_every must be positive.")
        log_path.parent.mkdir(parents=True, exist_ok=True)
        self.log_path = log_path
        self.progress_every = progress_every
        self._handle: TextIO = log_path.open("w", encoding="utf-8", buffering=1)

    @staticmethod
    def _timestamp() -> str:
        return datetime.now().astimezone().isoformat(timespec="seconds")

    def detail(self, message: str = "") -> None:
        """Write one detailed line to the log only."""
        self._handle.write(f"[{self._timestamp()}] {message}\n")

    def console(self, message: str = "") -> None:
        """Write a line both to the console and to the detailed log."""
        print(message, flush=True)
        self.detail(message)

    def progress(
        self,
        checked: int,
        total: int,
        fault_size: int,
    ) -> None:
        """Print only occasional progress updates, but retain all details in the log."""
        if checked == 1 or checked == total or checked % self.progress_every == 0:
            self.console(
                f"Progress: {checked}/{total} representatives checked "
                f"(current |F|={fault_size})."
            )

    def close(self) -> None:
        self._handle.close()


@dataclass(frozen=True)
class VerificationSummary:
    theorem: str
    case: str
    max_faults: int
    status: str
    verified: bool
    vertex_count: int
    edge_count: int
    automorphism_count: int
    representatives_by_fault_size: Dict[int, int]
    total_representatives: int
    checked_representatives: int
    elapsed_seconds: float
    counterexample: Optional[Dict[str, object]]


def run_verification_n4(
    max_faults: int,
    reporter: VerificationReporter,
    with_witness: bool = False,
) -> VerificationSummary:
    """Run the exact n=4 verification and return a machine-readable summary."""
    if max_faults < 0 or max_faults > M:
        raise ValueError(f"max_faults must lie between 0 and {M}.")

    start_time = time.time()
    reporter.console("Building automorphism group of AQ_4...")
    vertex_autos = compute_automorphisms()
    edge_autos = compute_edge_action_perms(vertex_autos)
    validate_automorphism_actions(vertex_autos, edge_autos)

    reporter.console(f"Number of vertex automorphisms found: {len(vertex_autos)}")
    reporter.console(f"Number of edges: {M}")
    reporter.console("Enumerating orbit representatives of fault sets...")

    reps_by_size = enumerate_fault_orbit_reps(max_faults, edge_autos)
    representatives_by_fault_size = {
        k: len(reps_by_size[k]) for k in range(max_faults + 1)
    }
    total_instances = sum(representatives_by_fault_size.values())

    for k, count in representatives_by_fault_size.items():
        reporter.console(f"  fault size {k}: {count} orbit representatives")
    reporter.console(f"Total representatives to check: {total_instances}")
    reporter.detail("")

    checked = 0
    for k in range(max_faults + 1):
        count_at_size = len(reps_by_size[k])
        reporter.console(
            f"Checking fault size {k}: {count_at_size} orbit representatives."
        )

        for local_index, fault_rep in enumerate(reps_by_size[k], start=1):
            checked += 1
            instance = precompute_instance(fault_rep)
            representative_text = fmt_edge_set(fault_rep, EDGE_LIST)

            for cycle_length in range(3, 9):
                witness = find_cycle_cover_witness(instance, cycle_length)
                if witness is None:
                    elapsed = time.time() - start_time
                    counterexample = {
                        "fault_size": k,
                        "fault_edge_ids": list(fault_rep),
                        "fault_edges": [list(EDGE_LIST[eid]) for eid in fault_rep],
                        "fault_edges_binary": representative_text,
                        "cycle_length": cycle_length,
                    }
                    reporter.console("[COUNTEREXAMPLE] Verification failed.")
                    reporter.console(f"Fault size: {k}")
                    reporter.console(f"Representative F = {representative_text}")
                    reporter.console(f"l = {cycle_length}")
                    return VerificationSummary(
                        theorem="Theorem 4.1",
                        case="n=4 base case",
                        max_faults=max_faults,
                        status="COUNTEREXAMPLE",
                        verified=False,
                        vertex_count=N,
                        edge_count=M,
                        automorphism_count=len(vertex_autos),
                        representatives_by_fault_size=representatives_by_fault_size,
                        total_representatives=total_instances,
                        checked_representatives=checked,
                        elapsed_seconds=elapsed,
                        counterexample=counterexample,
                    )

                validate_cycle_cover_witness(instance, cycle_length, witness)
                if with_witness:
                    cycle1, cycle2 = witness
                    reporter.detail(
                        f"PASS |F|={k} local={local_index}/{count_at_size} "
                        f"l={cycle_length} F={representative_text} "
                        f"C1={fmt_cycle(cycle1)} C2={fmt_cycle(cycle2)}"
                    )

            if not with_witness:
                reporter.detail(
                    f"PASS |F|={k} local={local_index}/{count_at_size} "
                    f"global={checked}/{total_instances} F={representative_text}"
                )
            reporter.progress(checked, total_instances, k)

        reporter.console(f"Completed fault size {k}.")

    elapsed = time.time() - start_time
    reporter.console("[VERIFIED] All orbit representatives passed.")
    reporter.console("Therefore Theorem 4.1 is verified for n = 4.")
    return VerificationSummary(
        theorem="Theorem 4.1",
        case="n=4 base case",
        max_faults=max_faults,
        status="VERIFIED",
        verified=True,
        vertex_count=N,
        edge_count=M,
        automorphism_count=len(vertex_autos),
        representatives_by_fault_size=representatives_by_fault_size,
        total_representatives=total_instances,
        checked_representatives=checked,
        elapsed_seconds=elapsed,
        counterexample=None,
    )


def verify_theorem_41_n4(
    max_faults: int = 5,
    verbose: bool = True,
    with_witness: bool = False,
) -> bool:
    """Compatibility wrapper returning only True or False."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = Path(f"theorem_4_1_n4_{timestamp}.log")
    reporter = VerificationReporter(
        log_path=log_path,
        progress_every=100 if verbose else 10**18,
    )
    try:
        summary = run_verification_n4(
            max_faults=max_faults,
            reporter=reporter,
            with_witness=with_witness,
        )
        return summary.verified
    finally:
        reporter.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify the n=4 base case of Theorem 4.1."
    )
    parser.add_argument(
        "--max-faults",
        type=int,
        default=5,
        help="Check all fault sizes from 0 through this value (default: 5).",
    )
    parser.add_argument(
        "--log-dir",
        default="verify_theorem_4_1_n4_logs",
        help="Directory for the detailed log and JSON summary.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help="Print one console progress update per this many representatives.",
    )
    parser.add_argument(
        "--with-witness",
        action="store_true",
        help="Store reconstructed cycle witnesses in the detailed log.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = Path(args.log_dir)
    log_path = log_dir / f"theorem_4_1_n4_{run_id}.log"
    summary_path = log_dir / f"theorem_4_1_n4_{run_id}_summary.json"
    reporter = VerificationReporter(
        log_path=log_path,
        progress_every=args.progress_every,
    )

    reporter.console(f"Detailed log: {log_path.resolve()}")
    reporter.console(f"Summary file: {summary_path.resolve()}")

    try:
        summary = run_verification_n4(
            max_faults=args.max_faults,
            reporter=reporter,
            with_witness=args.with_witness,
        )
        payload = asdict(summary)
        exit_code = 0 if summary.verified else 1
    except Exception as exc:
        reporter.console(f"[INCONCLUSIVE] Verification stopped: {exc}")
        payload = {
            "theorem": "Theorem 4.1",
            "case": "n=4 base case",
            "max_faults": args.max_faults,
            "status": "INCONCLUSIVE",
            "verified": False,
            "error": str(exc),
        }
        exit_code = 2
    finally:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with summary_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
        reporter.console(f"Detailed log saved to: {log_path.resolve()}")
        reporter.console(f"Summary saved to: {summary_path.resolve()}")
        reporter.close()

    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
