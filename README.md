# SPL Optimization - QAP Solver

This repository contains Python scripts for solving a small Quadratic Assignment
Problem (QAP) related to SPL optimization.

The problem is to assign 14 fixed 2D points to 14 named roles while minimizing a
weighted sum of Manhattan distances between selected role pairs.

## Problem Overview

The input data is embedded directly in the Python scripts:

- 14 fixed 2D points.
- 14 roles:
  - `a1` to `a6`
  - `b1` to `b3`
  - `c1` to `c3`
  - `d1` to `d2`
- 22 weighted objective terms.
- Manhattan distance, also called L1 distance, is used between assigned points.

Coordinates have one decimal digit. The solvers scale coordinates by 10 so that
all internal distance calculations can be done with integers.

## Objective

For every objective term `(weight, role_i, role_j)`, the solver adds:

```text
weight * ManhattanDistance(point_assigned_to_role_i, point_assigned_to_role_j)
```

The goal is to minimize the sum of all these terms.

## Files

| File | Description |
| --- | --- |
| `bb_solver.py` | Standalone pure-Python Branch and Bound solver. It does not require OR-Tools. |
| `exact_solver.py` | OR-Tools CP-SAT model using integer assignment variables and element constraints. |
| `exact_solver_v2.py` | Alternative CP-SAT model using a Boolean assignment encoding, with a Branch and Bound fallback if CP-SAT does not prove optimality. |
| `requirements.txt` | Python dependency list for the OR-Tools based solvers. |

## Requirements

- Python 3.8 or newer.
- `ortools`, only required for `exact_solver.py` and `exact_solver_v2.py`.

The standalone Branch and Bound solver only uses the Python standard library.

## Installation

Create and activate a virtual environment if desired, then install the
dependencies:

```bash
pip install -r requirements.txt
```

## Usage

Run the standalone Branch and Bound solver:

```bash
python bb_solver.py
```

Run the OR-Tools CP-SAT solvers:

```bash
python exact_solver.py
python exact_solver_v2.py
```

## Verified Result

On the current problem data, `bb_solver.py` proved the global optimum:

```text
f* = 141.8
```

The verified optimal assignment printed by `bb_solver.py` is:

| Role | Point |
| --- | --- |
| `a1` | `(10.5, 10.5)` |
| `a2` | `(4.5, 4.0)` |
| `a3` | `(13.5, 13.0)` |
| `a4` | `(3.0, 1.5)` |
| `a5` | `(4.7, 12.3)` |
| `a6` | `(4.7, 14.7)` |
| `b1` | `(10.5, 6.7)` |
| `b2` | `(1.3, 12.3)` |
| `b3` | `(13.5, 6.3)` |
| `c1` | `(10.5, 13.0)` |
| `c2` | `(1.5, 4.0)` |
| `c3` | `(1.3, 9.7)` |
| `d1` | `(4.7, 9.7)` |
| `d2` | `(1.3, 14.3)` |

This result was checked by running:

```bash
python bb_solver.py
```

During that run, the solver explored `3,725,985` nodes and verified the final
objective value as `141.8`. Runtime depends on the machine.

## Solver Notes

### Branch and Bound

`bb_solver.py` uses:

- A fixed role branching order.
- A current best upper bound.
- Lower-bound pruning based on minimum possible distances to remaining points.
- Exhaustive search with pruning to prove that no better assignment exists.

### CP-SAT

The CP-SAT scripts use Google OR-Tools to model the assignment problem exactly.
They are intended as alternative exact formulations, not as a packaged library
or command-line application.

## Limitations

- The problem data is hard-coded in the scripts.
- There is no separate configuration file or CLI argument parser.
- There is no automated test suite in this repository.
- The repository is a script-based research/optimization workspace, not a
  Python package.

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
