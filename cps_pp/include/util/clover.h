#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/clover.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: clover.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:29  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:16  anj
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
//  $RCSfile: clover.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/clover.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// clover.h
//
// C header file for the clover fermion library
//
//--------------------------------------------------------------------------

#ifndef INCLUDED_CLOVER_H
#define INCLUDED_CLOVER_H


CPS_END_NAMESPACE
#include<util/wilson.h>
CPS_START_NAMESPACE

//--------------------------------------------------------------------------
// Type definitions
//--------------------------------------------------------------------------
struct Clover {
  IFloat clover_coef;         // The clover coefficient 
  IFloat *frm_buf0;           // buffers of the size of checkerboarded 
  IFloat *frm_buf1;           //    fermionic fields for DiracOpClover.
  int nsites[4];             // number of local sites.  
  Wilson *wilson_p;          // Pointer to the wilson structure
};


//--------------------------------------------------------------------------
// initialization routine:
//   performs all initializations needed before clover library
// funcs are called. It sets the addressing related arrays and reserves
// memory for the needed temporary buffers. It only needs to be called only
// once at the begining of the program (or after a clover_end call)
// before any number of calls to clover functions are made.
// clover_p is a  pointer to Clover struct.
//--------------------------------------------------------------------------
void clover_init(Clover *clover_p);

//--------------------------------------------------------------------------
// This routine frees any memory reserved by clover_init.
// clover_p is a  pointer to Clover struct.
//--------------------------------------------------------------------------
void clover_end(Clover *clover_p);

//--------------------------------------------------------------------------
// calculate Y = A X   where  Y, X:   spinors (checkboarded)
//                               A:   compressed Clover matrices (chbd'ed)
//                               n:   number of local sites (chbd'ed)
// assume:  Y != X
//--------------------------------------------------------------------------
inline void clover_mat_mlt(IFloat *Y, const IFloat *A, const IFloat *X, int n);





//--------------------------------------------------------------------------
// inline function definition:  NO attention for user to pay here!
//--------------------------------------------------------------------------
extern "C" void clover_mat_mlt_asm(IFloat*, const IFloat*, const IFloat*, int n);
extern "C" void clover_mat_mlt_C(IFloat*, const IFloat*, const IFloat*, int n);

void clover_mat_mlt(IFloat *Y, const IFloat *A, const IFloat *X, int n) 
{
#ifdef _TARTAN
  clover_mat_mlt_asm(Y, A, X, n);
#else
  clover_mat_mlt_C(Y, A, X, n);
#endif
}





#endif  // #ifndef INCLUDED_CLOVER_H
CPS_END_NAMESPACE
