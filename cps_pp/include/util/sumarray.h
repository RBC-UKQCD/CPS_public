#ifndef SUM_ARRAY
#define SUM_ARRAY

#include "timer.h"
#include <complex>
#include <comms/sysfunc_cps.h>
//#include "crc32-reference.h"

#if 0
extern MPI_Comm QMP_COMM_WORLD;
#else
#define QMP_COMM_WORLD MPI_COMM_WORLD
#endif

inline int sumArray(long* recv, const long* send, const long n_elem) {
#ifdef USE_QMP
	return MPI_Allreduce((long*)send, recv, n_elem, MPI_LONG, MPI_SUM, QMP_COMM_WORLD);
#else
	memmove(recv, send, n_elem * sizeof(long));
	return 0;
#endif
}

inline int sumArray(uint32_t* recv, const uint32_t* send, const long n_elem) {
#ifdef USE_QMP
	return MPI_Allreduce((uint32_t*)send, recv, n_elem, MPI_UNSIGNED, MPI_SUM, QMP_COMM_WORLD);
#else
	memmove(recv, send, n_elem * sizeof(uint32_t));
	return 0;
#endif
}

inline int sumArray(int* recv, const int* send, const long n_elem) {
#ifdef USE_QMP
	return MPI_Allreduce((int*)send, recv, n_elem, MPI_INT, MPI_SUM, QMP_COMM_WORLD);
#else
	memmove(recv, send, n_elem * sizeof(int));
	return 0;
#endif
}

inline int sumArray(double* recv, const double* send, const long n_elem) {
#ifdef USE_QMP
	return MPI_Allreduce((double*)send, recv, n_elem, MPI_DOUBLE, MPI_SUM, QMP_COMM_WORLD);
#else
	memmove(recv, send, n_elem * sizeof(double));
	return 0;
#endif
}

template<class M>
int sumArray(M* vs, const long n_elem) {
	// M can be double or long
	int status = 0;
#ifdef USE_QMP
	M tmp[n_elem];
	status = sumArray(tmp, vs, n_elem);
	memcpy(vs, tmp, n_elem * sizeof(M));
#endif
	return status;
}

#endif
