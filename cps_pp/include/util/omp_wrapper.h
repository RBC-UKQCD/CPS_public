#ifndef OMP_WRAPPER_H
#define OMP_WRAPPER_H
#include <config.h>
#ifdef USE_OMP
#include <omp.h>
CPS_START_NAMESPACE
const int MAX_THREADS=64;
CPS_END_NAMESPACE
#else
inline int omp_get_num_threads(void) {return 1;}
inline int omp_get_thread_num(void) {return 0;}
inline void omp_set_num_threads(int n) {}
CPS_START_NAMESPACE
const int MAX_THREADS=1;
CPS_END_NAMESPACE
#endif
#endif
