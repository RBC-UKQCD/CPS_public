#ifndef BLAS_SUBS_H
#define BLAS_SUBS_H

#define USE_BLAS

#ifdef USE_BLAS
#include <util/qblas_extend.h>
#endif

#include<util/time_cps.h>

//#define DEBUG_DWF_DSLASH(msg,a ...) do \
//     printf("[%05d] %s:%d:QMP/%s(): " msg, UniqueID(), \
//                             __FILE__,__LINE__,__FUNCTION__,##a); \
//  while(0);

//#define DEBUG_DWF_DSLASH
#ifdef  DEBUG_DWF_DSLASH
#undef DEBUG_DWF_DSLASH
#define DEBUG_DWF_DSLASH(msg,a ...) do \
    if( UniqueID()%100==0 )			\
      printf("[%05d] " msg, UniqueID() ,##a);	\
  while(0);

static double time_now, time_prev;
static double time_elapse(){
  time_now = cps::dclock();
  double elp  = time_now - time_prev;
  time_prev=time_now;
  return elp;
}

#else
#define time_elapse() 0
#define DEBUG_DWF_DSLASH(msg,a ...) {}
#endif


#ifndef USE_BLAS
#define MOVE_FLOAT( pa, pb, n )  moveFloat(pa, pb, n)
#define VEC_TIMESEQU_FLOAT(py, fact, n ) vecTimesEquFloat( py, fact, n)
#define AXPY(n, fact, px, py)  fTimesV1PlusV2(py, fact, px, py, n)
#else
#define MOVE_FLOAT( pa, pb, n )  cblas_dcopy(n, pb, 1, pa, 1)
#define VEC_TIMESEQU_FLOAT(py, fact, n ) cblas_dscal( n,  fact, py,1 )
#define AXPY(n, fact, px, py)  cblas_daxpy(n, fact, px,1,py,1)
#endif



#endif
