#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:43:07 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/dense_matrix.h,v 1.3 2004-08-18 11:43:07 zs Exp $
//  $Id: dense_matrix.h,v 1.3 2004-08-18 11:43:07 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: dense_matrix.h,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/dense_matrix.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//
// lapack.h:
//     implement linear algebra package
//
// A general matrix is allocated as IFloat[row][column][real, imag], 
// the so-called canonical order in the QCDSP physics software system.
//
// Version 1.0:    
//


#ifndef INCLUDED_LAPACK_H
#define INCLUDED_LAPACK_H
 
CPS_END_NAMESPACE
#include<util/data_types.h>
CPS_START_NAMESPACE

enum MAT_INV_ALG { MAT_INV_ALG_LDL,     
		   // by first decomposing hermitian A into L D L^dag 
		   MAT_INV_ALG_LDL_CMPR  
		   // same as MAT_INV_ALG_LDL except that the hermitian 
		   // matrices are stored in compact format.
};

extern "C" {  

  void mat_inv(IFloat *out, const IFloat *in, int n, MAT_INV_ALG alg, 
	       IFloat *error_p);  
  //--------------------------------------------------------------------------
  // Purpose:
  //         B = 1/A.    
  // Arguments:
  //     A,B:   n x n complex matrix          (allocated as Float[n][n][2]
  //     alg:   algorithm specified to solve the inverse
  //     err:   sum { |A * B - I| } 
  //            (the absolute value of all the matrix elements)
  // Warning:
  //   in == out OK but never overlapped or embedded otherwise
  //--------------------------------------------------------------------------



  void mat_hrm_cmpr(IFloat *mat_out, const IFloat *mat_in, int mat_n);
  //--------------------------------------------------------------------------
  // Purpose:
  //   compress a hermitian matrix by storing only the lower triangular
  //   matrix elements in the same order.
  // Arguments:
  //   in:   n x n hermitian matrix allocated as IFloat[n*n*2]
  //   out:  compressed matrix allocated as IFloat[n*n]
  // Warning:
  //   in == out OK but never overlapped or embedded otherwise
  //--------------------------------------------------------------------------




  void mat_hrm_decm(IFloat *mat_out, const IFloat *mat_in, int mat_n);
  //--------------------------------------------------------------------------
  // Purpose:
  //   restore a hermitian matrix from its lower triangular half.
  // Arguments:
  //   out:   n x n hermitian matrix allocated as IFloat[n*n*2]
  //   in:    compressed matrix allocated as IFloat[n*n]
  // Warning:
  //   in == out OK, but never overlapped or embeded otherwise
  //--------------------------------------------------------------------------



  void mat_hrm_ldl(IFloat *L, IFloat *D, const IFloat *A, int n);
  //--------------------------------------------------------------------------
  // Purpose:
  //     to decompose a hermitian n by n matrix A into L D L^dag 
  // Attention:
  //    1. All matrices are stored in compact format for time efficiency.
  //    2. A, D and L must be different from each other.
  // Arguments:
  //    A:  lower triangular half of a hermitian matrix. allocated as Float[n*n]
  //    D:  real diagonal matrix. allocated as Float[n]
  //    L:  unit lower-triangular matrix. allocated as Float[n*(n-1)]
  //--------------------------------------------------------------------------

}     // extern "C"

#endif            /*  #ifndef INCLUDED_LAPACK_H   */




CPS_END_NAMESPACE
