# Computer verification for fault-tolerant 2-DCC pancyclicity of augmented cubes

This repository contains the Python programs used for the computer-assisted verification in the paper

**Fault-tolerant two-disjoint-cycle-cover pancyclicity of augmented cubes**.

## Files

- `theorem_4_1_n4_verify.py`  
  Verifies the base case \(n=4\) of Theorem 4.1.

- `subcase_2_4_i_n5_verify.py`  
  Verifies the exceptional case in Subcase 2.4(i) for \(n=5\).

- `subcase_3_4_i_n5_verify.py`  
  Verifies the exceptional case in Subcase 3.4(i) for \(n=5\).

- `subcase_3_4_i_n6_verify.py`  
  Verifies the exceptional case in Subcase 3.4(i) for \(n=6\).

## Requirements

The programs require Python 3 and the following packages:

```bash
pip install networkx scipy numpy
```

## How to run

Run the scripts from the command line:

```bash
python theorem_4_1_n4_verify.py
python subcase_2_4_i_n5_verify.py
python subcase_3_4_i_n5_verify.py
python subcase_3_4_i_n6_verify.py --jobs 9 --cycles-per-candidate 16
```

The last script may take longer than the others.

## Expected result

Each script should finish without reporting a counterexample.
