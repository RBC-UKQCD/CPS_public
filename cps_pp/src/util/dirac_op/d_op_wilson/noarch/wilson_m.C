#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.

  $Id: wilson_m.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/noarch/wilson_m.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: wilson_m.C,v 1.2 2003-07-24 16:53:54 zs Exp $
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
//  Revision 1.2  2001/06/19 18:12:47  anj
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
//  $RCSfile: wilson_m.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/noarch/wilson_m.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

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
  tmp1_f = wilson_p->af[0];

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
