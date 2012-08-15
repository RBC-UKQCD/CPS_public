#ifndef FAKE_OMP_H
#define FAKE_OMP_H
#ifdef USE_OMP
#include <omp.h>
#else
inline int omp_get_num_threads(void) {return 1;}
inline int omp_get_thread_num(void) {return 0;}
inline void omp_set_num_threads(int n) {}
#endif
#endif
