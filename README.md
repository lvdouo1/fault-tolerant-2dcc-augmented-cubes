[README.md](https://github.com/user-attachments/files/30231947/README.md)
# Computer verification for fault-tolerant 2-DCC pancyclicity of augmented cubes

This repository contains the Python programs used for the computer-assisted verification in the paper

**Fault-tolerant two-disjoint-cycle-cover pancyclicity of augmented cubes**.

## Files

- `theorem_4_1_n4_verify.py`  
  Verifies the base case `n = 4` of Theorem 4.1. This file is unchanged.

- `case_2_3_1_n5_verify.py`  
  Verifies the exceptional `n = 5` case in Case 2.3.1.

- `case_3_1_3_n5_verify.py`  
  Verifies the exceptional `n = 5` case in Case 3.1.3.

- `case_3_1_3_n6_verify.py`  
  Verifies the exceptional `n = 6` case in Case 3.1.3.

## Requirements

Python 3.10 or later is recommended. Install the required packages with

```bash
python -m pip install -r requirements.txt
```

or

```bash
python -m pip install networkx numpy "scipy>=1.9"
```

SciPy 1.9 or later is required because the programs use `scipy.optimize.milp`.

## How to run

```bash
python theorem_4_1_n4_verify.py
python case_2_3_1_n5_verify.py
python case_3_1_3_n5_verify.py
python case_3_1_3_n6_verify.py --jobs 9 --cycles-per-candidate 16
```

The `n = 6` program is substantially more expensive than the other programs. It writes progress logs to `verify_case_3_1_3_logs/` and a final summary to `case_3_1_3_n6_summary.json` by default.

## Meaning of solver statuses

The MILP-based programs accept only the following SciPy/HiGHS statuses as conclusive:

- `status == 0`: an optimal solution was found;
- `status == 2`: the MILP was proved infeasible.

A time limit, iteration limit, unbounded status, numerical problem, or any other termination status is treated as inconclusive. In that situation the program raises an error and does **not** claim that the corresponding case has been verified.

## Successful output

The two `n = 5` scripts exit with status code `0` only after the counterexample MILP has been proved infeasible.

The `n = 6` script exits with status code `0` only after all nine translation-symmetry subproblems have been proved infeasible. It exits with status code `1` if a genuine counterexample is found and with status code `2` if the computation terminates inconclusively.
