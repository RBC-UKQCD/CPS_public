#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_scl/wilson_mdag.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/***************************************************************************/
/*                                                                         */
/* wilson_mdag: It calculates chi = M^dag * psi, where M is the            */
/* Wilson fermion matrix. M is a function of the gauge fields u.           */
/* Kappa is the usual Wilson parameter and lx,ly,lz,lt is the lattice      */
/* size.                                                                   */
/*                                                                         */
/* This routine is to be used with scalar machines.                        */
/*                                                                         */
/* WARNING:                                                                 */
/*                                                                          */
/* This set of routines will work only if the node sublattices have         */
/* even number of sites in each direction.                                  */
/*                                                                          */
/***************************************************************************/

CPS_END_NAMESPACE
#include<util/data_types.h>
#include<util/wilson.h>
CPS_START_NAMESPACE


#define U(r,row,col,d,n,cb) *(u+(r+2*(row+3*(col+3*(d+4*(n+vol[0]*(cb)))))))
#define PSI(r,c,s,n)     *(psi+(r+2*(c+3*(s+4*(n)))))
#define CHI(r,c,s,n)     *(chi+(r+2*(c+3*(s+4*(n)))))
#define TMP1(r,c,s,n)     *(tmp1+(r+2*(c+3*(s+4*(n)))))

/*--------------------------------------------------------------------------*/
/* Function declarations                                                    */
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* wilson_mdag:                                                             */
/*--------------------------------------------------------------------------*/
void wilson_mdag(IFloat *chi_f, 
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
  tmp1_f = wilson_p->af[0];

  Float *chi = (Float *) chi_f;
  Float *psi = (Float *) psi_f;
  Float kappa = Float(kappa_f);

/*--------------------------------------------------------------------------*/
/* DslashDag_E0                                                             */
/*--------------------------------------------------------------------------*/
  wilson_dslash(tmp1_f, u_f, psi_f, 1, 1, wilson_p);

/*--------------------------------------------------------------------------*/
/* DslashDag_0E                                                             */
/*--------------------------------------------------------------------------*/
  wilson_dslash(chi_f, u_f, tmp1_f, 0, 1, wilson_p);

/*--------------------------------------------------------------------------*/
/* [1_OO - kappa * DslashDag_0E * DslashDag_E0] ]                           */
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
