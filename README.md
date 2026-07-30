Computer verification for edge-fault-tolerant 2-DCC pancyclicity of augmented cubes

This repository contains the computer-assisted verification programs accompanying the paper

Edge-fault-tolerant two-disjoint-cycle-cover pancyclicity of augmented cubes.

The programs verify the finite base and exceptional cases used in the proof of Theorem 4.1. The remaining cases are established by the mathematical arguments in the paper.

Repository: https://github.com/lvdouo1/edge-fault-tolerant-2dcc-augmented-cubes

Quick start

git clone https://github.com/lvdouo1/edge-fault-tolerant-2dcc-augmented-cubes.git
cd edge-fault-tolerant-2dcc-augmented-cubes
python -m pip install -r requirements.txt

Run the four verification programs with:

python theorem_4_1_n4_verify.py
python case_2_3_1_n5_verify.py
python case_3_1_3_n5_verify.py
python case_3_1_3_n6_verify.py --jobs 9 --cycles-per-candidate 16

The last computation is substantially more expensive than the other three.

Verification scope

Program

Case verified

Fault condition

Required cycle structure

theorem_4_1_n4_verify.py

Theorem 4.1, base case $n=4$

$F\subseteq E(AQ_4)$ and $

F

\le 5$

Two vertex-disjoint cycles of lengths $\ell$ and $16-\ell$, for every $3\le \ell\le 8$, covering all vertices

case_2_3_1_n5_verify.py

Theorem 4.1, Case 2.3.1, $n=5$

$

F_R

=6$ and $\delta(AQ_4-F_R)\ge 2$

A cycle of length $15$

case_3_1_3_n5_verify.py

Theorem 4.1, Case 3.1.3, $n=5$

$

F_R

=7$ and $\delta(AQ_4-F_R)\ge 2$

A cycle of length $15$

case_3_1_3_n6_verify.py

Theorem 4.1, Case 3.1.3, $n=6$

$

F_R

=9$ and $\delta(AQ_5-F_R)\ge 2$

A cycle of length $31$

For the base case, the required cycles $C_1$ and $C_2$ satisfy

$$|C_1|=\ell,\qquad |C_2|=16-\ell,$$

and

$$V(C_1)\cap V(C_2)=\varnothing,\qquadV(C_1)\cup V(C_2)=V(AQ_4).$$

Repository contents

File

Description

theorem_4_1_n4_verify.py

Exhaustive verification of the complete base case $n=4$

case_2_3_1_n5_verify.py

Verification of Case 2.3.1 for $n=5$

case_3_1_3_n5_verify.py

Verification of Case 3.1.3 for $n=5$

case_3_1_3_n6_verify.py

Parallel verification of Case 3.1.3 for $n=6$

requirements.txt

Python package requirements

Requirements

Python 3.10 or later is recommended. The programs use:

NetworkX;

NumPy;

SciPy, including scipy.optimize.milp and the HiGHS solver.

Install the dependencies with:

python -m pip install -r requirements.txt

For an archival run, also record the exact installed versions:

python --version
python -m pip freeze > requirements-lock.txt

Verification methods

Base case $n=4$

The program constructs $AQ_4$ as a Cayley graph on $\mathrm{GF}(2)^4$. Fault sets are grouped into orbits under affine automorphisms

$$x\longmapsto Ax+b,$$

where $A$ preserves the generating set of $AQ_4$. One representative from each orbit is checked.

For each representative fault set, an exact subset dynamic program determines whether every relevant induced vertex subset contains a Hamiltonian cycle. For each $\ell=3,\ldots,8$, the program searches for a subset $S$ with $|S|=\ell$ such that both

$$(AQ_4-F)[S]\quad\text{and}\quad(AQ_4-F)[V(AQ_4)\setminus S]$$

contain Hamiltonian cycles. Reconstructed witnesses are validated before the representative is accepted.

Exceptional cases $n=5$

A binary master MILP selects a candidate fault set satisfying the prescribed fault count and minimum-degree condition. Each previously discovered $15$-cycle yields a valid cut requiring at least one of its edges to be faulty.

For every candidate fault set, an exact Held--Karp dynamic program checks all possible omitted vertices. If a fault-free $15$-cycle is found, the corresponding cycle cut and all of its automorphic images are added to the master MILP. The verification ends successfully only when the master MILP is proved infeasible.

Exceptional case $n=6$

Translation symmetry reduces the search to nine subproblems, one for each generator edge $(0,s)$ of $AQ_5$. In each subproblem, the corresponding edge is forced to be faulty.

The master MILP selects nine faulty edges satisfying the degree condition and all previously generated $31$-cycle cuts. For each candidate fault set, the program removes each possible omitted vertex and checks Hamiltonicity by a second MILP. Degree-two constraints first generate a $2$-factor; disconnected subtours are eliminated iteratively until either a Hamiltonian cycle is found or infeasibility is proved.

All nine symmetry subproblems must be completed successfully.

Running the programs

Theorem 4.1, base case $n=4$

python theorem_4_1_n4_verify.py

The default run checks every fault size from $0$ through $5$. It creates timestamped output files in:

verify_theorem_4_1_n4_logs/

The generated files are:

theorem_4_1_n4_<timestamp>.log
theorem_4_1_n4_<timestamp>_summary.json

Useful options:

# Print progress every 500 representatives.
python theorem_4_1_n4_verify.py --progress-every 500

# Store reconstructed cycle witnesses in the detailed log.
python theorem_4_1_n4_verify.py --with-witness

# Choose another output directory.
python theorem_4_1_n4_verify.py --log-dir my_n4_logs

The option --max-faults is intended for development tests. The complete verification reported in the paper uses the default value 5.

Theorem 4.1, Case 2.3.1, $n=5$

python case_2_3_1_n5_verify.py

Theorem 4.1, Case 3.1.3, $n=5$

python case_3_1_3_n5_verify.py

Theorem 4.1, Case 3.1.3, $n=6$

python case_3_1_3_n6_verify.py --jobs 9 --cycles-per-candidate 16

The program writes one JSONL progress log per symmetry subproblem to

verify_case_3_1_3_logs/

and writes the final machine-readable summary to

case_3_1_3_n6_summary.json

Important options include:

--jobs N                    number of parallel worker processes
--cycles-per-candidate N    maximum number of cycle cuts per candidate
--master-time-limit SEC     optional time limit for each master MILP
--oracle-time-limit SEC     optional time limit for each Hamiltonian-cycle MILP
--max-iterations N          safety cap for each symmetry subproblem
--log-dir PATH              directory for JSONL logs
--summary-json PATH         path for the final JSON summary

A time limit is never interpreted as a successful verification.

Interpreting the results

A computation is accepted as successful only when it reaches its final verification message without an exception.

For the MILP-based programs, only the following SciPy/HiGHS statuses are treated as conclusive:

status == 0: an optimal feasible solution was found;

status == 2: infeasibility was proved.

Time limits, iteration limits, numerical problems, unbounded statuses, and all other solver terminations are treated as inconclusive.

For theorem_4_1_n4_verify.py and case_3_1_3_n6_verify.py, the exit codes are:

Exit code

Meaning

0

Verification completed successfully

1

A counterexample was found

2

The computation was inconclusive or an error occurred

The two $n=5$ programs return successfully only after the counterexample MILP has been proved infeasible. A counterexample or an inconclusive solver termination produces a nonzero process exit.

Recorded successful runs

All four programs have completed successfully for the full parameter ranges used in the paper.

Program

Recorded result

theorem_4_1_n4_verify.py

All $37{,}002$ fault-set orbit representatives passed

case_2_3_1_n5_verify.py

Verified after 5 master iterations and 640 cycle cuts

case_3_1_3_n5_verify.py

Verified after 4 master iterations and 512 cycle cuts

case_3_1_3_n6_verify.py

All nine translation-symmetry subproblems were proved infeasible

Base case $n=4$

The complete run constructed the 56-edge graph $AQ_4$, found 128 affine vertex automorphisms, and checked the following orbit representatives:

| Fault size $|F|$ | Orbit representatives ||---:|---:|| 0 | 1 || 1 | 3 || 2 | 36 || 3 | 343 || 4 | 3,593 || 5 | 33,026 || Total | 37,002 |

Every representative passed for every $\ell\in{3,4,5,6,7,8}$. The run also produced a detailed timestamped log and a JSON summary.

Exceptional case $n=6$

<details>
<summary>Recorded results for the nine symmetry subproblems</summary>

Generator $s$

Forced edge

Iterations

Cycle cuts

Wall time (s)

1

(0, 1)

17

304

139128.70

2

(0, 2)

15

272

59322.33

3

(0, 3)

16

288

86435.30

4

(0, 4)

15

272

60493.50

7

(0, 7)

16

288

61589.97

8

(0, 8)

15

272

64653.21

15

(0, 15)

16

288

63622.21

16

(0, 16)

15

272

83096.10

31

(0, 31)

16

288

60205.28

</details>

Wall-clock times depend strongly on the processor, operating system, Python environment, solver version, and number of parallel workers. They are included only as a reproducibility reference and are not part of the mathematical conclusion.

Reproducibility and archival records

For a publication or archival run, retain the following together:

the Git commit hash of the exact code used;

the Python version;

the exact dependency versions;

the operating system and CPU information;

the command-line arguments;

the complete console output;

the timestamped $n=4$ log and JSON summary;

the nine $n=6$ JSONL logs and final summary JSON.

A convenient command for recording the source revision is:

git rev-parse HEAD

The search procedures are exhaustive at the combinatorial-model level. The MILP programs use the floating-point HiGHS solver through SciPy; accordingly, the scripts validate returned binary solutions, set the relative MILP gap to zero, and reject every nonconclusive solver status.

Relation to the paper

These programs verify only the finite base and exceptional cases explicitly identified in the proof of Theorem 4.1. They are intended to be read together with the mathematical proof in the paper, which establishes all remaining cases and explains why the computationally verified statements are sufficient.

Citation

Please cite the accompanying paper when referring to this repository:

Edge-fault-tolerant two-disjoint-cycle-cover pancyclicity of augmented cubes.

Complete bibliographic information can be added after publication.
