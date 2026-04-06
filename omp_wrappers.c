/* Wrappers for OpenMP functions called from f2c-translated code.
   f2c mangles names containing underscores by appending __, so
   omp_get_wtime becomes omp_get_wtime__. */
#ifdef _OPENMP
#include <omp.h>
double omp_get_wtime__(void) { return omp_get_wtime(); }
#else
double omp_get_wtime__(void) { return 0.0; }
#endif
