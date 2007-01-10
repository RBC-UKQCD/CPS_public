#ifndef INCLUDED_VECTOR_ASM_H
#define INCLUDED_VECTOR_ASM_H               //!< Prevent multiple inclusion
#include<config.h>
/*!\file
  \brief  Declaration/definition of Assembly routines used in Vector class


  $Id: vector_asm.h,v 1.1 2007-01-10 18:30:36 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2007-01-10 18:30:36 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/vector_asm.h,v 1.1 2007-01-10 18:30:36 chulwoo Exp $
//  $Id: vector_asm.h,v 1.1 2007-01-10 18:30:36 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: vector_asm.h,v $
//  $Revision: 1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/vector_asm.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


CPS_START_NAMESPACE


class Vector; 	// forward declaration
class Matrix;


#if TARGET == QCDOC ||  TARGET == BGL || ( defined USE_QMP )

extern "C" {
  void xaxpy (Float *scalep,Float *InOutScale, Float *Add, int len);
  void xaxpy_norm (Float *scalep, Float *InOutScale, Float *Add, int len, Float *res);
  
  void vaxpy3(Vector *res,Float *scale,Vector *mult,Vector *add, int ncvec);
  inline void vaxpy3_m(Matrix *res,Float *scale,Matrix *mult,const Matrix *add, int ncvec){
    vaxpy3((Vector *)res, scale, (Vector *)mult,(Vector *)add,ncvec);
  }

  void vaxpy3_norm(Vector *res,Float *scale,Vector *mult,Vector *add, int ncvec,Float *norm);
  inline void vaxpy3_norm_m(Vector *res,Float *scale,Vector *mult,Vector *add, int ncvec,Float *norm){
    vaxpy3_norm((Vector *)res, scale, (Vector *)mult,(Vector *)add,ncvec,norm);
  }

  // SU(3) * SU(3) streaming operations
  void m1m2(Matrix *res, const Matrix *m1, const Matrix *m2, int *length);
  void m1m2dag(Matrix *res, const Matrix *m1, const Matrix *m2, int *length);
  void m1dagm2dag(Matrix *res, const Matrix *m1, const Matrix *m2, int *length);
  void gdagmdag(Matrix *res, const Matrix *g, const Matrix *m, int *length);

  // f*SU(3)*SU(3) + SU(3) streaming operations
  void fm1m2pm3(Matrix *res, Float *scale, const Matrix *m1, 
		const Matrix *m2, Matrix *m3, int *length);
  void fm1m2dagpm3(Matrix *res, Float *scale, const Matrix *m1, 
		   const Matrix *m2, Matrix *m3, int *length);
  void fm1dagm2dagpm3(Matrix *res, Float *scale, const Matrix *m1, 
		      const Matrix *m2, Matrix *m3, int *length);
  void fgdagm1dagpm2(Matrix *res, Float *scale, const Matrix *g, 
		     const Matrix *m1, Matrix *m2, int *length);

}

#endif
CPS_END_NAMESPACE
#endif
