# LP Solver (Simplex Method)

A linear programming solver that reads LP problems from standard input and outputs optimal solutions, detecting infeasible and unbounded cases.

## Features

- **Dual Simplex Method** for handling initially infeasible dictionaries
- **Steepest Edge Pivot Rule** with Bland's rule fallback after 3 degenerate pivots
- **Primal-Dual Two-Phase Method** for initially infeasible LPs
- Dictionary-form pivot operations

## Usage

```bash
# Pipe input
cat input.txt | python simplex_alg.py

# Interactive
python simplex_alg.py
# Then paste LP into stdin
```

## Dependencies

- Python 3
- NumPy
