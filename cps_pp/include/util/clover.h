#include<config.h>
CPS_START_NAMESPACE
/*! \file
  \brief Declarations of routines used internally in the DiracOpClover class.

  $Id: clover.h,v 1.5 2004-09-02 16:59:27 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-09-02 16:59:27 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/clover.h,v 1.5 2004-09-02 16:59:27 zs Exp $
//  $Id: clover.h,v 1.5 2004-09-02 16:59:27 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: clover.h,v $
//  $Revision: 1.5 $
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
#include <util/wilson.h>
CPS_START_NAMESPACE

//--------------------------------------------------------------------------
// Type definitions
//--------------------------------------------------------------------------
//! A container of data relevant to the clover matrix multiplication.
/*!
  This contains parameters describing the  local lattice geometry and
  pointers to workspace arrays needed by the Wilson-clover matrix
  multiplication routines.
 */

struct Clover {
  IFloat clover_coef;         //!< The clover coefficient 
  IFloat *frm_buf0;           //!< Workspace array 
  IFloat *frm_buf1;           //!< Workspace array 
  int nsites[4];             //!< The dimensions of the local lattice.
  Wilson *wilson_p;          //!< Wilson hopping term data
};


//--------------------------------------------------------------------------
//! Initialisation of parameters and memory used in the clover matrix multiplication.
//--------------------------------------------------------------------------
void clover_init(Clover *clover_p);

//--------------------------------------------------------------------------
//! Frees memory reserved by clover_init    
//--------------------------------------------------------------------------
void clover_end(Clover *clover_p);

//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
inline void clover_mat_mlt(IFloat *Y, const IFloat *A, const IFloat *X, int n);





//--------------------------------------------------------------------------
// inline function definition:  NO attention for user to pay here!
//--------------------------------------------------------------------------

extern "C" void clover_mat_mlt_asm(IFloat*, const IFloat*, const IFloat*, int n);
extern "C" void clover_mat_mlt_C(IFloat*, const IFloat*, const IFloat*, int n);

//--------------------------------------------------------------------------
//! Multiplication by the clover term.
/*!
  A spinor vector on a single parity is multiplied by the clover matrix.

  \param Y The resulting vector
  \param A The (lower triangular part of the) clover matrix 
  \param X The vector to be multiplied
  \param n ?
*/
//--------------------------------------------------------------------------
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
