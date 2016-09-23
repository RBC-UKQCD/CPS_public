#include <config.h>
#ifdef USE_SSE
#include "../sse/wilson_mdagm.C"
#else
CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.

  $Id: wilson_mdagm.C,v 1.5 2011/03/04 11:25:28 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2011/03/04 11:25:28 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_wilson/noarch/wilson_mdagm.C,v 1.5 2011/03/04 11:25:28 chulwoo Exp $
//  $Id: wilson_mdagm.C,v 1.5 2011/03/04 11:25:28 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_wilson/noarch/wilson_mdagm.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/***************************************************************************/
/*                                                                         */
/* wilson_mdagm: It calculates chi = M^dag * M * psi, where M is the       */
/* Wilson fermion matrix. M is a function of the gauge fields u.           */
/* The sum |M*psi|^2 is also calulated and stored in mp_sq_p.              */
/* Kappa is the usual Wilson parameter and lx,ly,lz,lt is the lattice      */
/* size.                                                                   */
/*                                                                         */
/* This routine is to be used with scalar machines.                        */
/*                                                                         */
/* WARNING:                                                                */
/*                                                                         */
/* This set of routines will work only if the node sublattices have        */
/* even number of sites in each direction.                                 */
/*                                                                         */
/***************************************************************************/

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/wilson.h>
CPS_START_NAMESPACE


//! Access to the elements of the \e SU(3) matrix
/*!
  Gets the element of the \e SU(3) matrix \e u with row \e row,
  column \e col and complex component \e d
*/
#define U(r,row,col,d,n,cb) *(u+(r+2*(row+3*(col+3*(d+4*(n+vol[0]*(cb)))))))
//! Access to the elements of a spinor vector.
/*!
  Gets the element of the spinor \e psi with spin \e s,
  colour \e c and complex component \e r
*/
#define PSI(r,c,s,n)     *(psi+(r+2*(c+3*(s+4*(n)))))
//! As above, but the vector is called chi
#define CHI(r,c,s,n)     *(chi+(r+2*(c+3*(s+4*(n)))))
#define TMP1(r,c,s,n)     *(tmp1+(r+2*(c+3*(s+4*(n)))))
#define TMP2(r,c,s,n)     *(tmp2+(r+2*(c+3*(s+4*(n)))))


/*--------------------------------------------------------------------------*/

/*!
  The matrix-vector product   \f$ M^\dagger M\psi \f$ is computed
  on all sites with odd parity.

  \pre The gauge field and spinor vector should be in odd-even order.
  
  \param chi_f The matrix-vector product
  \param u_f The gauge field
  \param psi_f The vector to be mutliplied
  \param mp_sq_p_f If this is initially non-NULL, then \f$ |M\psi|^2 \f$ is
  computed and stored at this address.
  \param kappa_f The Wilson matrix hopping parameter
  \param wilson_p The Wilson multiplication data structure.
 */

void wilson_mdagm(IFloat *chi_f, 
		  IFloat *u_f, 
		  IFloat *psi_f, 
		  IFloat *mp_sq_p_f,
		  IFloat kappa_f,
		  Wilson *wilson_p)
{
  IFloat *tmp1_f;
  IFloat *tmp2_f;
  Float sum;
  int vol;
  int r, c, s, n;


/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  vol =  wilson_p->vol[0];
  tmp1_f = wilson_p->af[0];
  tmp2_f = wilson_p->af[1];

  Float *chi = (Float *) chi_f;
  Float *psi = (Float *) psi_f;
  Float *mp_sq_p = (Float *) mp_sq_p_f;
  Float kappa = Float(kappa_f);
  Float *tmp2 = (Float *) tmp2_f;


/*--------------------------------------------------------------------------*/
/* Dslash_E0                                                                */
/*--------------------------------------------------------------------------*/
  wilson_dslash(tmp1_f, u_f, psi_f, 1, 0, wilson_p);

/*--------------------------------------------------------------------------*/
/* Dslash_0E                                                                */
/*--------------------------------------------------------------------------*/
  wilson_dslash(tmp2_f, u_f, tmp1_f, 0, 0, wilson_p);

/*--------------------------------------------------------------------------*/
/* 1_OO - kappa^2 * Dslash_0E * Dslash_E0                                   */
/*--------------------------------------------------------------------------*/
  sum = 0.0;
  for(n=0;n<vol;n++){
    for(s=0;s<4;s++){
      for(c=0;c<3;c++){
	for(r=0;r<2;r++){
	  TMP2(r,c,s,n) = ( PSI(r,c,s,n) - kappa*kappa * TMP2(r,c,s,n));
	  if(mp_sq_p != 0)
	    sum = sum + TMP2(r,c,s,n)*TMP2(r,c,s,n);
	}
      }
    }
  }
  if(mp_sq_p != 0)
    *mp_sq_p = sum;

/*--------------------------------------------------------------------------*/
/* DslashDag_E0                                                             */
/*--------------------------------------------------------------------------*/
  wilson_dslash(tmp1_f, u_f, tmp2_f, 1, 1, wilson_p);

/*--------------------------------------------------------------------------*/
/* DslashDag_0E                                                             */
/*--------------------------------------------------------------------------*/
  wilson_dslash(chi_f, u_f, tmp1_f, 0, 1, wilson_p);

/*--------------------------------------------------------------------------*/
/* [1_OO - kappa * DslashDag_0E * DslashDag_E0] *                           */
/*                                 [1_OO - kappa^2 * Dslash_0E * Dslash_E0] */
/*--------------------------------------------------------------------------*/
  for(n=0;n<vol;n++){
    for(s=0;s<4;s++){
      for(c=0;c<3;c++){
	for(r=0;r<2;r++){
	  CHI(r,c,s,n) = ( TMP2(r,c,s,n) - kappa*kappa * CHI(r,c,s,n));
	}
      }
    }
  }


}







CPS_END_NAMESPACE
#endif
