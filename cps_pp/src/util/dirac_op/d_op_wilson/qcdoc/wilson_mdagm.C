#include<config.h>
#include<stdio.h>
CPS_START_NAMESPACE
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

#define U(r,row,col,d,n,cb) *(u+(r+2*(row+3*(col+3*(d+4*(n+vol[0]*(cb)))))))
#define PSI(r,c,s,n)     *(psi+(r+2*(c+3*(s+4*(n)))))
#define CHI(r,c,s,n)     *(chi+(r+2*(c+3*(s+4*(n)))))
#define TMP1(r,c,s,n)     *(tmp1+(r+2*(c+3*(s+4*(n)))))
#define TMP2(r,c,s,n)     *(tmp2+(r+2*(c+3*(s+4*(n)))))


/*--------------------------------------------------------------------------*/
/* Function declarations                                                    */
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* wilson_mdagm:                                                            */
/*--------------------------------------------------------------------------*/
#define WHPAD 512
void wilson_mdagm(Float *chi_f, 
		  Float *u_f, 
		  Float *psi_f, 
		  Float *mp_sq_p_f,
		  Float kappa_f,
		  Wilson *wilson_p)
{
  Float *tmp1_f;
  Float *tmp2_f;
  Float sum;
  int vol;
  int r, c, s, n;
  printf("wilson_mdagm: mp_sq_p_f = %p\n",mp_sq_p_f);

/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  vol =  wilson_p->vol[0];
  printf("wilson_p=%p vol=%d\n",wilson_p,vol);
  tmp1_f = wilson_p->spinor_tmp;
  tmp2_f  = tmp1_f + SPINOR_SIZE * vol + (WHPAD/sizeof(Float)) ;

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
  if(mp_sq_p != 0){
   for(n=0;n<vol;n++){
    for(s=0;s<4;s++){
     for(c=0;c<3;c++){
      for(r=0;r<2;r++){
  	  TMP2(r,c,s,n) = ( PSI(r,c,s,n) - kappa*kappa * TMP2(r,c,s,n));
	    sum = sum + TMP2(r,c,s,n)*TMP2(r,c,s,n);
            printf("sum(%d,%d,%d,%d) = %e\n",r,c,s,n,sum);
      }
     }
    }
   }
    *mp_sq_p = sum;
    printf("mp_sq_p=%e\n",sum);
  }

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
