# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

High order layer potential evaluation library in R^2 (chunkmatc). Fortran 77 numerical library for solving Laplace and Helmholtz boundary integral equations using chunk-based discretization with specialized singular quadratures.

## Build Commands

Each test/driver has its own `.make` file. The executable is always named `int2`. Building automatically runs the executable.

```bash
# Build and run a test (e.g., test8)
make -f test8.make

# Clean object files
make -f test8.make clean

# Build the quadrature driver
make -f hqsuppquad.make
```

**Compiler selection:** Edit the `HOST` variable in the `.make` file (uncomment one):
- `linux-gfortran` / `linux-gfortran-openmp`
- `macos-gfortran` / `macos-gfortran-openmp`

Compiler flags are set in `make.inc`. All builds use `gfortran -std=legacy -O3`. On macOS, the linker gets `-Wl,-stack_size,0x20000000` for a 512 MB stack.

## Architecture

### Core Library

- **chunkmatc.f** — Central discretizer and geometry engine. User-defined geometry via callback subroutines (`chunkpnt`). Builds interaction matrices chunk-by-chunk (`chunkmatc`, `chunkmatc_od` for off-diagonal, `chunkmatc_dd` for diagonal blocks). Geometry refinement via `genrefineinfo` variants (uniform, dyadic).

- **inter2dn.f** — Interaction kernel library. Provides Laplace (`lfinter1`–`lfinter4`) and Helmholtz (`hfinter1`–`hfinter4`) kernels: single layer (S), double layer (D), and hypersingular variants. Kernels encode geometry as 6-element arrays `[x, y, nx, ny, tx, ty]`.

- **hqsuppquad.f** — Specialized quadrature tables for smooth, logarithmic, Hilbert, and hypersingular kernels, up to 10th order. This is the largest source file (~12K lines of quadrature data).

### Supporting Routines

- **adapgauss.f** — Adaptive Gaussian quadrature on an interval (scalar functions)
- **cadavect.f** — Adaptive Gaussian integration on segments (vector/smooth functions)
- **cgmres6-rel.f** — Complex GMRES iterative solver with cyclic buffer
- **cqrsolve.f** — Complex QR decomposition solver (`cqrsolv`, `cqrdecom`, `cqrsolve`, `cqrinv`)
- **hank103.f** — Hankel function H_0^(1) and H_1^(1) computation
- **legeexps.f** — Legendre polynomial expansions and Gauss-Legendre nodes/weights
- **prini.f** — Formatted printing/debugging utilities

### Driver Programs

- **hqsuppquad_dr.f** — Driver for testing/generating quadrature tables. **Reads `npols` from stdin** (requires interactive input).

### Test Programs

All tests solve Helmholtz BVPs on a circle with a unit charge, comparing direct potential against solved density potential:

- **test8.f** — Dirichlet, first kind integral equation
- **test8dn.f** — Same + Dirichlet-Neumann map check
- **test9.f** — Dirichlet, second kind integral equation
- **test9dn.f** — Same + D-N map (requires hypersingular evaluation)

### Typical Solver Flow

1. Set parameters (wavenumber `rk`, polynomial order `norder`)
2. Get Legendre nodes via `legeexps()`
3. Define geometry callback, generate chunks via `chunkgeo()`
4. Build discretization points/normals/weights via `chunkallpnts()`, `chunkallwhts()`
5. Assemble interaction matrix via `chunkmatc()`
6. Solve with GMRES (`cgmres`) or QR (`cqrsolve`)
7. Evaluate/compare potentials

### Conventions

- Error codes propagated via integer `ier` (first argument)
- Extensible parameter passing via `par1`–`par4` arguments
- Geometry callbacks receive chunk index and return position, derivatives, normals
- Output files use Fortran unit numbers (16–18 for geometry, 117–118 for normals/tangents)
- **Naming note:** The comments at the top of `chunkmatc.f` list `fpatchpnt`, `spatchpnt`, etc., but the actual subroutine definitions (and test files) use `fchunkpnt`, `schunkpnt`, `echunkpnt`, `rchunkpnt`. The comments are stale.

## MATLAB

`matlab/hqsuppquad*.m` — MATLAB versions of the auxiliary quadrature tables (orders 1–10). These provide nodes and weights for smooth, log-singular, Hilbert, and hypersingular kernels.
