#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/dense_matrix.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: dense_matrix.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:30  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:17  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: dense_matrix.h,v $
//  $Revision: 1.1.1.1 $
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
