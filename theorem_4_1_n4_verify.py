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
       cycle_exists[S] = True and cycle_exists[V \ S] = True.
   Such a pair gives the desired two disjoint cycles covering all vertices.

Because the property is invariant under graph automorphisms, checking one
representative from each orbit of fault sets is sufficient.

Notes
-----
- This script is intended as a verification/certification tool for the paper.
- It prints the number of orbit representatives at each fault size and confirms
  whether every representative passes.
- Optional witness reconstruction is included.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

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
M = len(EDGE_LIST)  # should be 56


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


# ---------------------------------------------------------------------------
# Main verification
# ---------------------------------------------------------------------------

def verify_theorem_41_n4(
    max_faults: int = 5,
    verbose: bool = True,
    with_witness: bool = False,
) -> bool:
    """
    Verify the n=4 base case of Theorem 4.1.

    Because the property is monotone under deleting fewer faulty edges, we check all
    edge-fault sets of sizes 0..5. To keep the run feasible, we check one representative
    from each automorphism orbit of fault sets.
    """
    if verbose:
        print("Building automorphism group of AQ_4...")
    vertex_autos = compute_automorphisms()
    edge_autos = compute_edge_action_perms(vertex_autos)

    if verbose:
        print(f"Number of vertex automorphisms found: {len(vertex_autos)}")
        print(f"Number of edges: {M}")
        print("Enumerating orbit representatives of fault sets...")

    reps_by_size = enumerate_fault_orbit_reps(max_faults, edge_autos)

    total_instances = 0
    for k in range(max_faults + 1):
        if verbose:
            print(f"  fault size {k}: {len(reps_by_size[k])} orbit representatives")
        total_instances += len(reps_by_size[k])

    if verbose:
        print(f"Total representatives to check: {total_instances}\n")

    checked = 0
    for k in range(max_faults + 1):
        for F_rep in reps_by_size[k]:
            checked += 1
            instance = precompute_instance(F_rep)

            for l in range(3, 9):
                witness = find_cycle_cover_witness(instance, l)
                if witness is None:
                    print("FAILED")
                    print(f"Fault size: {k}")
                    print(f"Representative F = {fmt_edge_set(F_rep, EDGE_LIST)}")
                    print(f"l = {l}")
                    return False

                if verbose and with_witness:
                    c1, c2 = witness
                    print(
                        f"PASS  |F|={k}  l={l}  "
                        f"F={fmt_edge_set(F_rep, EDGE_LIST)}  "
                        f"C1={fmt_cycle(c1)}  C2={fmt_cycle(c2)}"
                    )

            if verbose and not with_witness:
                print(f"[{checked}/{total_instances}] passed: |F|={k}, F={fmt_edge_set(F_rep, EDGE_LIST)}")

    if verbose:
        print("\nAll representatives passed.")
        print("Therefore Theorem 4.1 is verified for n = 4.")
    return True


if __name__ == "__main__":
    verify_theorem_41_n4(verbose=True, with_witness=False)
