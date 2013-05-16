#include<config.h>
CPS_START_NAMESPACE
/*! \file
  \brief Declarations of routine used internally in the DiracOpWilson class.

  $Id: naive.h,v 1.1 2013-05-16 04:14:50 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-05-16 04:14:50 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/naive.h,v 1.1 2013-05-16 04:14:50 chulwoo Exp $
//  $Id: naive.h,v 1.1 2013-05-16 04:14:50 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: naive.h,v $
//  $Revision: 1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/naive.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/****************************************************************************/
/* 10/16/97                                                                 */
/*                                                                          */
/* f_naive.h                                                                 */
/*                                                                          */
/* C header file for the fermion wilson library wilson.lib                  */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_NAIVE_H
#define INCLUDED_NAIVE_H

CPS_END_NAMESPACE
#ifdef USE_QMP
#include <qmp.h>
#endif
#include <util/wilson.h>
#include <util/data_types.h>
CPS_START_NAMESPACE

/*==========================================================================*/
/* The Wilson fermion matrix * vector routines                              */
/*==========================================================================*/

extern "C"{

//! Multiplication by the odd-even preconditioned Wilson matrix    
/*!
  Performs the multiplication \f$ \chi = M^\dagger M \psi \f$.
  
  \param chi The result vector
  \param u The gauge field
  \param psi The vector to be multiplied
  \param m_sq_p The value of \f$ |M \psi|^2 \f$
  \param kappa The Wilson matrix mass parameter
  \param wilson_p pointer to a Wilson struct.
*/
void naive_mdagm(IFloat  *chi,    
		  IFloat  *u,      
		  IFloat  *psi,    
		  IFloat  *mp_sq_p,
		  IFloat  Kappa,   
		  Wilson *wilson_p);

//--------------------------------------------------------------------

//! Multiplication by the Wilson matrix hopping term
/*!
  The multiplication is performed on all sites of a single parity.
  The resulting vector is then defined on sites of the opposite parity.
  
  \pre The gauge field and spinor vector should be in odd-even order.  

  \param chi The resulting vector
  \param u The gauge field
  \param psi The vector bing multiplied
  \param cb Multiplies a vector on odd parity sites if this is 0, even
            parity sites if this is 1.
  \param dag Multiply by the hermitian conjugate hopping term if this is 1,
             or not if this is 0.
  \param wilson_p The Wilson multiplication data structure.
 */
    
void naive_dslash(IFloat *chi, 
		   IFloat *u, 
		   IFloat *psi, 
		   int cb,
		   int dag,
		   Wilson *wilson_p);

void naive_dslash_two(Float *chi0, Float *chi1,
		   Float *u, 
		   Float *psi0, Float *psi1, 
		   int cb0, int cb1,
		   int dag,
		   Wilson *wp);

//--------------------------------------------------------------------

//! Multiplication by the odd-even preconditioned Wilson matrix
/*!
  The matrix-vector multiplication is performed on all sites with odd
  parity.

  \pre The gauge field and spinor vector should be in odd-even order.
  
  \param chi The matrix-vector product
  \param u The gauge field
  \param psi The vector to be mutliplied
  \param kappa The Wilson matrix hopping parameter
  \param wilson_p The Wilson multiplication data structure.
 */
    
void naive_m(IFloat *chi, 
	      IFloat *u, 
	      IFloat *psi, 
	      IFloat kappa,
	      Wilson *wilson_p);

//----------------------------------------------------------------------

//! Multiplication by the odd-even preconditioned Wilson matrix
/*!
  The matrix-vector multiplication is performed on all sites of odd
   parity.

  \pre The gauge field and spinor vector should be in odd-even order.
  
  \param chi The matrix-vector product
  \param u The gauge field
  \param psi The vector to be mutliplied
  \param kappa The Wilson matrix hopping parameter
  \param wilson_p The Wilson multiplication data structure.
 */
    
void naive_mdag(IFloat *chi, 
		 IFloat *u, 
		 IFloat *psi, 
		 IFloat kappa,
		 Wilson *wilson_p);

}

#endif

CPS_END_NAMESPACE
