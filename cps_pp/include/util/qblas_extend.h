#ifndef __QBLAS_EXTEND__CD
#define __QBLAS_EXTEND__CD
#ifdef HAVE_CBLAS_H
#include <cblas.h>
#elif (defined HAVE_MKL_CBLAS_H)
#include <mkl_cblas.h>
#elif (defined HAVE_GSL_CBLAS_H)
#include <gsl_cblas.h>
#endif


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

inline void cblas_zdotc_sub(const int N, const double *X, const double *Y,
			  double* dot)
{ return cblas_zdotc_sub(N,X,1,Y,1, dot); }

//inline void cblas_cdotc_sub(const int N, const float *X, const float *Y,
//		  float* dot)
//{ return cblas_cdotc_sub(N,X,Y,dot); }



#endif 
