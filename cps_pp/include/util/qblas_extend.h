#ifndef __QBLAS_EXTEND__CD
#define __QBLAS_EXTEND__CD
#include "qblas.h"

/*
  some convenience functions for 
  the case where all vectors are 
  continguous in memory
*/
inline void cblas_daxpy(const int N, 
                        const double alpha, 
                        const double *X,
                        double *Y)
{ cblas_daxpy(N,alpha,X,1,Y,1); }

inline double cblas_dnrm2(const int N, 
                          const double *X)
{ return cblas_dnrm2(N,X,1); }

inline void cblas_dcopy(const int N, 
                        const double *X, 
                        double *Y )
{ cblas_dcopy(N,X,1,Y,1); }

inline void cblas_dscal(const int N, 
                        const double alpha, 
                        double *X)
{ cblas_dscal(N,alpha,X,1); }

inline double cblas_ddot(const int N, 
			 const double *X, 
			 const double *Y )
{ return cblas_ddot(N,X,1,Y,1); }

inline double cblas_ddot(const int N, const double *X)
{ return cblas_ddot(N,X,1,X,1); }


#endif 
