# test7 timing results

## Machine

- Apple M4 Max, 16 cores, 128 GB RAM
- macOS 26.4, arm64

## Compilers

- gfortran 15.2.0 (Homebrew GCC 15.2.0_1)
- gcc-14 14.3.0 (Homebrew GCC 14.3.0)
- f2c 20240504-mpfr-2.2.0

## Problem

Laplace equation in R^2, interior Dirichlet BVP via double-layer potential.
Chunked boundary discretization, norder=10, noversamp=40, npts=400.

## Results (best of 3 runs)

### 4 threads

| Configuration          | chunkmatc (s) | cqrdecom (s) | error      |
|------------------------|---------------|---------------|------------|
| f2c non-OMP            | 0.0614        | 0.0292        | 0.147E-13  |
| f2c OMP (4 thr)        | 0.0171        | 0.0290        | 0.147E-13  |
| f2c quad (no OMP)      | 2.470         | 4.498         | 0.200E-31  |
| f2c OMP quad (4 thr)   | 0.704         | 4.627         | 0.200E-31  |
| gfortran non-OMP       | 0.0573        | 0.0381        | 0.153E-13  |
| gfortran OMP (4 thr)   | 0.0163        | 0.0376        | 0.153E-13  |
| gfortran quad (no OMP) | 1.900         | 4.794         | 0.200E-31  |
| gfortran OMP quad (4 thr) | 0.540      | 4.922         | 0.200E-31  |

### 8 threads

| Configuration              | chunkmatc (s) | cqrdecom (s) | error      |
|----------------------------|---------------|---------------|------------|
| f2c OMP (8 thr)            | 0.0088        | 0.029         | 0.147E-13  |
| f2c OMP quad (8 thr)       | 0.356         | 4.575         | 0.200E-31  |
| gfortran OMP (8 thr)       | 0.0084        | 0.038         | 0.153E-13  |
| gfortran OMP quad (8 thr)  | 0.272         | 4.885         | 0.200E-31  |

## Notes

- chunkmatc is OpenMP-parallelized; cqrdecom is serial (no speedup from threads).
- Quad precision uses software __float128 (libquadmath); ~40x slower on chunkmatc, ~155x on cqrdecom.
- Quad accuracy: ~17 extra digits (error ~2E-32 vs ~1.5E-14).
- chunkmatc scales ~2x from 4 to 8 threads.
- gfortran quad OMP timers required omp_wrappers.c (single-underscore shim) to override
  libgomp's omp_get_wtime_, which returns double even under -freal-8-real-16.
- f2c quad OMP timers required omp_wrappers.c (double-underscore shim) compiled with -DF2C_FLOAT128.
