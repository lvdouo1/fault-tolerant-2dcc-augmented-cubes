# Computer verification for edge-fault-tolerant 2-DCC pancyclicity of augmented cubes

This repository contains the Python programs used for the computer-assisted verification in the paper

> **Edge-fault-tolerant two-disjoint-cycle-cover pancyclicity of augmented cubes**.

The programs verify the finite base and exceptional cases that arise in the proof of Theorem 4.1. All searches are exact: exhaustive dynamic programming, graph-automorphism reduction, and mixed-integer linear programming (MILP) are used to rule out counterexamples.

## Verified statements

The repository verifies the following cases.

### Theorem 4.1, base case `n = 4`

For every fault set $F \subseteq E(AQ_4)$ with $|F| \le 5$, and every integer $3 \le \ell \le 8$, the graph $AQ_4-F$ contains two vertex-disjoint cycles $C_1$ and $C_2$ such that

$$
|C_1|=\ell, \qquad |C_2|=16-\ell,
$$

and

$$
V(C_1)\cup V(C_2)=V(AQ_4).
$$

### Theorem 4.1, Case 2.3.1, exceptional case `n = 5`

For every six-edge fault set $F_R \subseteq E(AQ_4)$ satisfying

$$
\delta(AQ_4-F_R)\ge 2,
$$

the graph $AQ_4-F_R$ contains a cycle of length 15.

### Theorem 4.1, Case 3.1.3, exceptional case `n = 5`

For every seven-edge fault set $F_R \subseteq E(AQ_4)$ satisfying

$$
\delta(AQ_4-F_R)\ge 2,
$$

the graph $AQ_4-F_R$ contains a cycle of length 15.

### Theorem 4.1, Case 3.1.3, exceptional case `n = 6`

For every nine-edge fault set $F_R \subseteq E(AQ_5)$ satisfying

$$
\delta(AQ_5-F_R)\ge 2,
$$

the graph $AQ_5-F_R$ contains a cycle of length 31.

## Repository contents

| File | Purpose | Verification method |
|---|---|---|
| `theorem_4_1_n4_verify.py` | Verifies the complete base case `n = 4` for all fault sizes from 0 through 5 and all $\ell=3,\ldots,8$. | Affine-automorphism orbit reduction and exact subset dynamic programming. |
| `case_2_3_1_n5_verify.py` | Verifies Case 2.3.1 for `n = 5`. | Counterexample MILP, exact Held--Karp dynamic programming, and automorphism-orbit cycle cuts. |
| `case_3_1_3_n5_verify.py` | Verifies Case 3.1.3 for `n = 5`. | Counterexample MILP, exact Held--Karp dynamic programming, and automorphism-orbit cycle cuts. |
| `case_3_1_3_n6_verify.py` | Verifies Case 3.1.3 for `n = 6`. | Nine translation-symmetry subproblems, a master counterexample MILP, and exact Hamiltonian-cycle MILPs with subtour elimination. |
| `requirements.txt` | Lists the required Python packages. | — |

## Requirements

Python 3.10 or later is recommended. The required packages are:

- `networkx>=2.8`
- `numpy>=1.23`
- `scipy>=1.9`

SciPy 1.9 or later is required because the MILP-based programs use `scipy.optimize.milp` and the HiGHS solver.

Install the dependencies with:

```bash
python -m pip install -r requirements.txt
```

## Clone the repository

```bash
git clone https://github.com/lvdouo1/edge-fault-tolerant-2dcc-augmented-cubes.git
cd edge-fault-tolerant-2dcc-augmented-cubes
```

## Running the verification programs

Run all four programs from the repository directory:

```bash
python theorem_4_1_n4_verify.py
python case_2_3_1_n5_verify.py
python case_3_1_3_n5_verify.py
python case_3_1_3_n6_verify.py --jobs 9 --cycles-per-candidate 16
```

The `n = 6` computation is substantially more expensive than the other three computations. The number of parallel processes can be changed with `--jobs`.

## Options for the `n = 4` program

The default command

```bash
python theorem_4_1_n4_verify.py
```

checks every fault size from 0 through 5. It automatically creates timestamped files in

```text
verify_theorem_4_1_n4_logs/
```

including:

```text
theorem_4_1_n4_<timestamp>.log
theorem_4_1_n4_<timestamp>_summary.json
```

Useful optional arguments are:

```bash
# Change how often progress is printed to the console.
python theorem_4_1_n4_verify.py --progress-every 500

# Store reconstructed cycle witnesses in the detailed log.
python theorem_4_1_n4_verify.py --with-witness

# Choose another output directory.
python theorem_4_1_n4_verify.py --log-dir my_n4_logs

# Test a smaller range during development only.
python theorem_4_1_n4_verify.py --max-faults 3
```

For the full verification used in the paper, keep the default `--max-faults 5`.

## Options and output files for the `n = 6` program

A typical full run is:

```bash
python case_3_1_3_n6_verify.py --jobs 9 --cycles-per-candidate 16
```

By default, the program writes per-subproblem JSONL logs to

```text
verify_case_3_1_3_logs/
```

and writes the final machine-readable summary to

```text
case_3_1_3_n6_summary.json
```

Important options include:

```text
--jobs N                    number of parallel symmetry subproblems
--cycles-per-candidate N    maximum cycle cuts generated per candidate
--master-time-limit SEC     optional limit for each master MILP
--oracle-time-limit SEC     optional limit for each Hamiltonian-cycle MILP
--max-iterations N          safety cap for each symmetry subproblem
--log-dir PATH              directory for JSONL progress logs
--summary-json PATH         path for the final JSON summary
```

A solver time limit is never interpreted as a successful verification. If a limit is reached, the run is reported as inconclusive.

## Verification methods

### Base case `n = 4`

The program constructs $AQ_4$ as a Cayley graph on $\mathrm{GF}(2)^4$. Fault sets are grouped into orbits under affine graph automorphisms of the form

$$
x\longmapsto Ax+b,
$$

where $A$ preserves the generating set of $AQ_4$. Only one representative from each orbit is checked.

For every representative fault set, an exact subset dynamic program determines whether each induced vertex subset contains a Hamiltonian cycle. For each $\ell=3,\ldots,8$, the program searches for a subset $S$ with $|S|=\ell$ such that both $(AQ_4-F)[S]$ and $(AQ_4-F)[V\setminus S]$ contain Hamiltonian cycles. The reconstructed witnesses are independently validated before the representative is accepted.

### Exceptional cases `n = 5`

A binary master MILP selects a candidate fault set satisfying the required fault count and minimum-degree condition. Every previously discovered 15-cycle contributes a valid cut requiring at least one of its edges to be faulty.

For each candidate, an exact Held--Karp dynamic program checks all possible omitted vertices and searches for a Hamiltonian cycle on the remaining 15 vertices. Whenever a fault-free 15-cycle is found, its cut and all automorphic images are added to the master MILP. Verification is complete only when the master MILP is proved infeasible.

### Exceptional case `n = 6`

Translation symmetry reduces the search to nine subproblems, one for each generator edge $(0,s)$ of $AQ_5$. Each subproblem forces the corresponding edge to be faulty.

The master MILP selects nine faulty edges satisfying the degree condition and all previously generated 31-cycle cuts. For each candidate, the program removes each possible omitted vertex and solves Hamiltonicity exactly by a second MILP. Degree-two constraints first produce a 2-factor; disconnected subtours are then eliminated iteratively until either a Hamiltonian cycle is found or infeasibility is proved.

All nine translation-symmetry subproblems must be proved infeasible before the program reports the `n = 6` case as verified.

## Solver statuses and exit codes

The MILP-based programs regard only the following SciPy/HiGHS statuses as conclusive:

- `status == 0`: an optimal feasible solution was found;
- `status == 2`: infeasibility was proved.

A time limit, iteration limit, numerical issue, unbounded status, or any other termination status is treated as inconclusive and is never reported as verification success.

For `theorem_4_1_n4_verify.py` and `case_3_1_3_n6_verify.py`, the exit codes are:

| Exit code | Meaning |
|---:|---|
| `0` | The requested verification completed successfully. |
| `1` | A genuine counterexample was found. |
| `2` | The computation ended inconclusively or an error occurred. |

The two `n = 5` scripts exit successfully only after the counterexample MILP has been proved infeasible. A counterexample or inconclusive solver termination produces a nonzero exit.

## Recorded successful runs

All four programs have been run successfully for the full parameter ranges used in the paper.

| Program | Recorded result |
|---|---|
| `theorem_4_1_n4_verify.py` | Verified all 37,002 orbit representatives for fault sizes 0 through 5. |
| `case_2_3_1_n5_verify.py` | Verified after 5 master iterations and 640 cycle cuts. |
| `case_3_1_3_n5_verify.py` | Verified after 4 master iterations and 512 cycle cuts. |
| `case_3_1_3_n6_verify.py` | All nine translation-symmetry subproblems were proved infeasible. |

### Recorded `n = 4` base-case run

The successful full run constructed the 56-edge graph $AQ_4$, found 128 affine vertex automorphisms, and checked the following numbers of fault-set orbit representatives:

| Fault size $|F|$ | Orbit representatives |
|---:|---:|
| 0 | 1 |
| 1 | 3 |
| 2 | 36 |
| 3 | 343 |
| 4 | 3,593 |
| 5 | 33,026 |
| **Total** | **37,002** |

All 37,002 representatives passed for every $\ell=3,4,\ldots,8$. Thus the complete base case of Theorem 4.1 for $n=4$ was verified. The program also saved a timestamped detailed log and a machine-readable JSON summary.

### Recorded `n = 6` exceptional-case run

The recorded `n = 6` subproblem results were:

| Generator $s$ | Forced edge | Iterations | Cycle cuts | Wall time (s) |
|---:|---|---:|---:|---:|
| 1 | `(0, 1)` | 17 | 304 | 139128.70 |
| 2 | `(0, 2)` | 15 | 272 | 59322.33 |
| 3 | `(0, 3)` | 16 | 288 | 86435.30 |
| 4 | `(0, 4)` | 15 | 272 | 60493.50 |
| 7 | `(0, 7)` | 16 | 288 | 61589.97 |
| 8 | `(0, 8)` | 15 | 272 | 64653.21 |
| 15 | `(0, 15)` | 16 | 288 | 63622.21 |
| 16 | `(0, 16)` | 15 | 272 | 83096.10 |
| 31 | `(0, 31)` | 16 | 288 | 60205.28 |

Running times depend strongly on the processor, operating system, Python environment, and solver version. They are reported only as a reproducibility reference and are not part of the mathematical conclusion.

## Reproducibility

For an archival run, retain the following together:

1. the exact source-code revision or Git commit hash;
2. the output of `python --version`;
3. the installed package versions, for example from `python -m pip freeze`;
4. the complete console output;
5. the timestamped `n = 4` log and JSON summary;
6. the nine `n = 6` JSONL logs and final summary JSON.

The computations are exact searches at the model level, but the MILP programs use the floating-point HiGHS solver through SciPy. The scripts therefore validate returned binary solutions, set the MILP relative gap to zero, and reject every nonconclusive solver status.

## Citation

When referring to this repository, please cite the accompanying paper:

> *Edge-fault-tolerant two-disjoint-cycle-cover pancyclicity of augmented cubes*.

Complete bibliographic information can be added here after publication.
