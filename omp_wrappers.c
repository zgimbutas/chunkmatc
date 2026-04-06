/* omp_wrappers.c — Name-mangling shims for OpenMP 2.5 runtime routines.

   f2c appends __ to Fortran subroutine/function names, so Fortran calls
   to omp_get_wtime() become calls to omp_get_wtime__() in C. These
   wrappers bridge the gap.

   Compile with: gcc -fopenmp -c omp_wrappers.c
   Link with f2c-generated code that uses OpenMP runtime routines.

   Covers all 16 runtime routines from OpenMP Fortran 2.0/2.5 spec. */

#ifdef _OPENMP
#include <omp.h>

/* ---- Execution environment routines ---- */

void omp_set_num_threads__(int *n) { omp_set_num_threads(*n); }
int omp_get_num_threads__(void) { return omp_get_num_threads(); }
int omp_get_max_threads__(void) { return omp_get_max_threads(); }
int omp_get_thread_num__(void) { return omp_get_thread_num(); }
int omp_get_num_procs__(void) { return omp_get_num_procs(); }
int omp_in_parallel__(void) { return omp_in_parallel(); }
void omp_set_dynamic__(int *flag) { omp_set_dynamic(*flag); }
int omp_get_dynamic__(void) { return omp_get_dynamic(); }
void omp_set_nested__(int *flag) { omp_set_nested(*flag); }
int omp_get_nested__(void) { return omp_get_nested(); }

/* ---- Lock routines ---- */

void omp_init_lock__(omp_lock_t *lock) { omp_init_lock(lock); }
void omp_destroy_lock__(omp_lock_t *lock) { omp_destroy_lock(lock); }
void omp_set_lock__(omp_lock_t *lock) { omp_set_lock(lock); }
void omp_unset_lock__(omp_lock_t *lock) { omp_unset_lock(lock); }
int omp_test_lock__(omp_lock_t *lock) { return omp_test_lock(lock); }

void omp_init_nest_lock__(omp_nest_lock_t *lock) { omp_init_nest_lock(lock); }
void omp_destroy_nest_lock__(omp_nest_lock_t *lock) { omp_destroy_nest_lock(lock); }
void omp_set_nest_lock__(omp_nest_lock_t *lock) { omp_set_nest_lock(lock); }
void omp_unset_nest_lock__(omp_nest_lock_t *lock) { omp_unset_nest_lock(lock); }
int omp_test_nest_lock__(omp_nest_lock_t *lock) { return omp_test_nest_lock(lock); }

/* ---- OpenMP 3.0 routines ---- */

void omp_set_max_active_levels__(int *n) { omp_set_max_active_levels(*n); }
int omp_get_max_active_levels__(void) { return omp_get_max_active_levels(); }
int omp_get_level__(void) { return omp_get_level(); }
int omp_get_ancestor_thread_num__(int *level) { return omp_get_ancestor_thread_num(*level); }
int omp_get_team_size__(int *level) { return omp_get_team_size(*level); }
int omp_get_active_level__(void) { return omp_get_active_level(); }
int omp_get_thread_limit__(void) { return omp_get_thread_limit(); }

/* ---- OpenMP 3.1 routines ---- */

int omp_in_final__(void) { return omp_in_final(); }

/* ---- Timing routines (f2c: double underscore) ---- */

#ifdef F2C_FLOAT128
/* Under -freal-8-real-16, f2c declares these as returning quadreal
   (__float128).  The wrapper must match the promoted ABI.
   Use a volatile intermediate to force the double->__float128 promotion
   through memory, preventing gcc from optimizing away the widening. */
__float128 omp_get_wtime__(void) { volatile double t = omp_get_wtime(); return (__float128)t; }
__float128 omp_get_wtick__(void) { volatile double t = omp_get_wtick(); return (__float128)t; }
#else
double omp_get_wtime__(void) { return omp_get_wtime(); }
double omp_get_wtick__(void) { return omp_get_wtick(); }
#endif

/* ---- Timing routines (gfortran: single underscore) ---- */
/* gfortran calls omp_get_wtime_ (Fortran name mangling).  Under
   -freal-8-real-16, gfortran expects __float128 return but libgomp's
   omp_get_wtime_ returns double.  Override it here; we call
   omp_get_wtime (C symbol, no underscore) to avoid recursion. */

#ifdef F2C_FLOAT128
__float128 omp_get_wtime_(void) { volatile double t = omp_get_wtime(); return (__float128)t; }
__float128 omp_get_wtick_(void) { volatile double t = omp_get_wtick(); return (__float128)t; }
#else
double omp_get_wtime_(void) { return omp_get_wtime(); }
double omp_get_wtick_(void) { return omp_get_wtick(); }
#endif

#else
/* Stub implementations for non-OpenMP builds */

void omp_set_num_threads__(int *n) { (void)n; }
int omp_get_num_threads__(void) { return 1; }
int omp_get_max_threads__(void) { return 1; }
int omp_get_thread_num__(void) { return 0; }
int omp_get_num_procs__(void) { return 1; }
int omp_in_parallel__(void) { return 0; }
void omp_set_dynamic__(int *flag) { (void)flag; }
int omp_get_dynamic__(void) { return 0; }
void omp_set_nested__(int *flag) { (void)flag; }
int omp_get_nested__(void) { return 0; }

void omp_init_lock__(void *lock) { (void)lock; }
void omp_destroy_lock__(void *lock) { (void)lock; }
void omp_set_lock__(void *lock) { (void)lock; }
void omp_unset_lock__(void *lock) { (void)lock; }
int omp_test_lock__(void *lock) { (void)lock; return 1; }

void omp_init_nest_lock__(void *lock) { (void)lock; }
void omp_destroy_nest_lock__(void *lock) { (void)lock; }
void omp_set_nest_lock__(void *lock) { (void)lock; }
void omp_unset_nest_lock__(void *lock) { (void)lock; }
int omp_test_nest_lock__(void *lock) { (void)lock; return 1; }

void omp_set_max_active_levels__(int *n) { (void)n; }
int omp_get_max_active_levels__(void) { return 1; }
int omp_get_level__(void) { return 0; }
int omp_get_ancestor_thread_num__(int *level) { (void)level; return 0; }
int omp_get_team_size__(int *level) { (void)level; return 1; }
int omp_get_active_level__(void) { return 0; }
int omp_get_thread_limit__(void) { return 1; }
int omp_in_final__(void) { return 1; }

#ifdef F2C_FLOAT128
__float128 omp_get_wtime__(void) { return 0.0Q; }
__float128 omp_get_wtick__(void) { return 0.0Q; }
__float128 omp_get_wtime_(void) { return 0.0Q; }
__float128 omp_get_wtick_(void) { return 0.0Q; }
#else
double omp_get_wtime__(void) { return 0.0; }
double omp_get_wtick__(void) { return 0.0; }
double omp_get_wtime_(void) { return 0.0; }
double omp_get_wtick_(void) { return 0.0; }
#endif

#endif
