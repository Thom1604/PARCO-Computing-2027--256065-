#ifndef HEAT_H
#define HEAT_H

#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>


#ifdef _OPENMP
  #include <omp.h>
#endif

// Macro per indicizzazione 1D di una griglia (N+2) x (N+2) 
#define INDEX(i, j, N) ((i) * ((N) + 2) + (j))

// ================= Predefined values ================= 
#define DEFAULT_N 256
#define DEFAULT_MAX_ITER 10000
#define DEFAULT_TOL 1e-4
#define DEFAULT_BOUNDARY_TEMP 100.0

// ================= Timing ============================
static inline double now_sec(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#endif
}

#endif // HEAT_H