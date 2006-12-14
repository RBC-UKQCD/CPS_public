#include <config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
*/
/*--------------------------------------------------------------------*/

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


void wilson_m(IFloat *chi_f, 
	      IFloat *u_f, 
	      IFloat *psi_f, 
	      IFloat kappa_f,
	      Wilson *wilson_p)
{
  IFloat *tmp1_f;
  int vol;
  int r, c, s, n;

/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  vol =  wilson_p->vol[0];
  tmp1_f = wilson_p->spinor_tmp;

  Float *chi = (Float *) chi_f;
  Float *psi = (Float *) psi_f;
  Float kappa = Float(kappa_f);
  

/*--------------------------------------------------------------------------*/
/* Dslash_E0                                                                */
/*--------------------------------------------------------------------------*/
  wilson_dslash(tmp1_f, u_f, psi_f, 1, 0, wilson_p);

/*--------------------------------------------------------------------------*/
/* Dslash_0E                                                                */
/*--------------------------------------------------------------------------*/
  wilson_dslash(chi_f, u_f, tmp1_f, 0, 0, wilson_p);

/*--------------------------------------------------------------------------*/
/* 1_OO - kappa^2 * Dslash_0E * Dslash_E0                                   */
/*--------------------------------------------------------------------------*/
  for(n=0;n<vol;n++){
    for(s=0;s<4;s++){
      for(c=0;c<3;c++){
	for(r=0;r<2;r++){
	  CHI(r,c,s,n) = ( PSI(r,c,s,n) - kappa*kappa * CHI(r,c,s,n));
	}
      }
    }
  }

}







CPS_END_NAMESPACE
