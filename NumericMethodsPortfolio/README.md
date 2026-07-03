# Numerical Methods Laboratory

From-scratch implementations of classical numerical methods in **C**, built around custom `Matrix.h`, `Vector`, and `Solve_by.h` libraries that evolve across assignments. Each module is self-contained, compiles with a single `gcc` call, and reads input data from plain text files.

> **Course:** Numerical Methods — CIMAT, M.Sc. in Computer Science.

---

## Overview

The repository spans the core curriculum of a graduate-level numerical methods course, organized into six thematic modules:

- **Linear Systems** — Direct solvers (Gaussian elimination, Crout LU, Cholesky) and iterative solvers (Jacobi, Gauss-Seidel), plus a 1D heat equation application with band-diagonal Cholesky.
- **Eigenvalue Problems** — Power method, inverse power method, subspace iteration (direct and inverse), Jacobi rotation for symmetric matrices, and conjugate gradient.
- **Interpolation** — Taylor, Lagrange, Neville, Newton divided differences, Hermite, natural cubic spline, and clamped cubic spline.
- **Least Squares & Integration** — Least-squares fitting with polynomial, cosine, and RBF bases; Newton-Cotes (open and closed), and Gaussian quadrature.
- **Nonlinear Systems** — Fixed-point iteration, Newton's method, and Broyden's method for systems of nonlinear equations.
- **Ordinary Differential Equations** — Euler, Heun (predictor-corrector), Taylor second-order, and Runge-Kutta 4th-order, for scalar ODEs and 2×2 systems.

---

## Repository Structure

```
.
├── linear-systems/
│   ├── AguilarKalebMetodosTarea4.c
│   ├── Matrix.h
│   ├── Solve_by.h
│   ├── A.txt, b.txt, SPD.txt, b_spd.txt, ...
│   └── README.txt
│
├── eigenvalue-problems/
│   ├── AguilarKalebMetodosTarea5.c
│   ├── AguilarKalebMetodosTarea6.c
│   ├── Matrix.h
│   ├── Solve_by.h
│   ├── Eigen_3x3.txt, Eigen_5x5.txt, ...
│   └── README.txt
│
├── interpolation/
│   ├── AguilarKalebMetodosTarea8.c
│   ├── AguilarKalebMetodosTarea9.c
│   ├── Taylor.txt, ValoresExp.txt, IV.txt
│   ├── P1a.txt, P1b.txt, Pato.txt, Snoopy.txt
│   └── README.txt
│
├── least-squares-and-integration/
│   ├── AguilarKalebMetodosTarea10.c
│   ├── Matrix.h
│   ├── Solve_by.h
│   ├── IntegracionNumerica.h
│   ├── ValoresSutherland1.txt, ValoresSutherland2.txt
│   └── README.txt
│
├── nonlinear-systems/
│   ├── AguilarKalebMetodosTarea11.c
│   ├── Matrix.h
│   ├── Solve_by.h
│   └── README.txt
│
├── ode-solvers/
│   ├── AguilarKalebMetodosTarea12.c
│   ├── Matrix.h
│   ├── Solve_by.h
│   ├── SolveEDO.h
│   └── README.txt
│
└── README.md
```

---

## Custom Libraries

All modules share a set of hand-written C libraries that grow across assignments:

**`Matrix.h`** — Defines `Matrix` and `Vector` structs with dynamic memory allocation. Provides creation from files, arithmetic (add, subtract, multiply, transpose), norms (Euclidean, Frobenius), and I/O utilities. All computations use `long double` precision.

**`Solve_by.h`** — Cascade solver that attempts Cholesky → Crout LU → Gauss-Seidel → Jacobi automatically. Also exposes individual solvers, Gaussian elimination with partial pivoting, QR factorization, and conjugate gradient.

**`IntegracionNumerica.h`** — Newton-Cotes quadrature (open n=0–3, closed n=1–4) and Gaussian quadrature (n=1–5) with hardcoded nodes and weights.

**`SolveEDO.h`** — ODE stepping functions (Euler, Heun, Taylor-2, RK4) for scalar and 2×2 systems, with a generic driver that iterates over [a,b] with N steps.

---

## Methods Summary

| Module | Methods |
|--------|---------|
| Linear Systems | Gaussian elimination, Crout LU, Cholesky, Jacobi, Gauss-Seidel, band-diagonal Cholesky |
| Eigenvalues | Power, inverse power, subspace iteration, inverse subspace iteration, Jacobi rotation, conjugate gradient |
| Interpolation | Taylor, Lagrange, Neville, Newton, Hermite, cubic spline (natural), cubic spline (clamped) |
| Least Squares & Integration | Polynomial / cosine / RBF bases, Newton-Cotes (open & closed), Gaussian quadrature |
| Nonlinear Systems | Fixed-point, Newton, Broyden |
| ODE Solvers | Euler, Heun, Taylor 2nd-order, RK4, 2×2 system extensions |

---

## Compilation

Every module compiles with a single `gcc` call. No external dependencies beyond the standard C library and `math.h`:

```bash
gcc AguilarKalebMetodosTarea4.c -o programa -lm
```

Each program supports a default mode (runs all problems) and individual problem selection via command-line arguments. See the per-module `README.txt` for usage examples.

---

## Design Choices

- **No external dependencies.** All linear algebra, interpolation, and integration routines are implemented from scratch — no LAPACK, no GSL.
- **`long double` precision** throughout for numerical stability in iterative methods and high-order interpolation.
- **File-driven I/O.** Input matrices and vectors are read from plain text files; results are printed to terminal and optionally saved to CSV/TXT for plotting.
- **Modular headers.** `Matrix.h` and `Solve_by.h` are reused and extended across assignments, simulating a growing numerical library.

---

## Notes

- Each module has its own `README.txt` with detailed function signatures, CLI usage, and file format specifications.
- The libraries evolve across assignments (e.g., `Matrix.h` gains `long double` support, in-place operations, and column extraction in later modules).
- Output CSV files are designed for direct plotting with gnuplot, matplotlib, or similar tools.

---

## Author

**Kaleb Aguilar**

Physics Engineer
M.Sc. Student in Computer Science at CIMAT

LinkedIn:
https://linkedin.com/in/kaleb-aguilar-70971129a
