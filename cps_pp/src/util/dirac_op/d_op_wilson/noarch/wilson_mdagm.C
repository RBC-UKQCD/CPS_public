#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/noarch/wilson_mdagm.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: wilson_mdagm.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:22  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:48  anj
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
//  Revision 1.2  2001/05/25 06:16:06  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: wilson_mdagm.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/noarch/wilson_mdagm.C,v $
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
#include<util/data_types.h>
#include<util/wilson.h>
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
