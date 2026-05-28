# SPL Optimization - QAP Solver

This project solves a 14-point Quadratic Assignment Problem (QAP): assign fixed
2D points to roles `a1` to `d2` while minimizing a weighted sum of Manhattan
distances between selected role pairs.

## Contents

- `bb_solver.py`: standalone pure-Python Branch and Bound solver.
- `exact_solver.py`: exact OR-Tools CP-SAT model using integer assignment
  variables and element constraints.
- `exact_solver_v2.py`: alternative CP-SAT implementation with a Boolean
  encoding and an optional Branch and Bound fallback.

## Requirements

- Python 3.8+
- OR-Tools, required for the CP-SAT solvers

Install dependencies with:

```bash
pip install -r requirements.txt
```

## Usage

Run the pure-Python solver:

```bash
python bb_solver.py
```

Run the CP-SAT solvers:

```bash
python exact_solver.py
python exact_solver_v2.py
```

## Problem Data

- 14 fixed 2D points.
- 14 roles: `a1-a6`, `b1-b3`, `c1-c3`, `d1-d2`.
- 22 weighted objective terms.
- Coordinates are scaled by 10 internally for integer calculations.

## License

MIT
