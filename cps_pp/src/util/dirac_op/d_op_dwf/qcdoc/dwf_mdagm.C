#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-07-01 21:22:37 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdoc/dwf_mdagm.C,v 1.2 2004-07-01 21:22:37 chulwoo Exp $
//  $Id: dwf_mdagm.C,v 1.2 2004-07-01 21:22:37 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.4.1  2004/06/09 04:31:54  chulwoo
//  *** empty log message ***
//
//  Revision 1.1.2.1  2004/05/20 14:36:32  pab
//  Files for three optimised dirac operators.
//  Added some patches to the ReadLattice too.
//
//  Revision 1.1.1.1  2003/06/22 13:34:46  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.2  2001/06/19 18:12:39  anj
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
//  Revision 1.2  2001/05/25 06:16:05  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: dwf_mdagm.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdoc/dwf_mdagm.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// dwf_mdagm.C
//
// dwf_mdagm M^dag M where M is the fermion matrix.
// The in, out fields are defined on the checkerboard lattice.
// <out, in> = <dwf_mdagm*in, in> is returned in dot_prd.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/dirac_op.h>
#include "dwf_internal.h"
CPS_START_NAMESPACE


void dwf_mdagm(Vector *out, 
	       Matrix *gauge_field, 
	       Vector *in, 
	       Float *dot_prd,
	       Float mass, 
	       Dwf *dwf_lib_arg)
{
  Vector *frm_tmp1 = (Vector *) dwf_lib_arg->frm_tmp1;
  int f_size = 24 * dwf_lib_arg->vol_4d * dwf_lib_arg->ls / 2;
  Float minus_kappa_sq = -dwf_lib_arg->dwf_kappa * dwf_lib_arg->dwf_kappa;
  Vector *frm_tmp2 = (Vector *) dwf_lib_arg->frm_tmp2;
  

//------------------------------------------------------------------
// Apply M frm_tmp1 <- in
//------------------------------------------------------------------

  //  dwf_m(frm_tmp1, gauge_field, in, mass, dwf_lib_arg);

  //------------------------------------------------------------------
  // Apply Dslash E <- O
  //------------------------------------------------------------------
  dwf_dslash(frm_tmp2, gauge_field, in, mass, 1, 0, dwf_lib_arg);

  //------------------------------------------------------------------
  // Apply Dslash O <- E
  //------------------------------------------------------------------
  dwf_dslash(frm_tmp1, gauge_field, frm_tmp2, mass, 0, 0, dwf_lib_arg);

//------------------------------------------------------------------
// Calculate the dot product <out, in> = <M in, M in>
//------------------------------------------------------------------
  if(dot_prd != 0){

    vaxpy3_norm((Float *)frm_tmp1,
		&minus_kappa_sq,
	       (Float *)frm_tmp1,
	       (Float *)in,f_size/6,dot_prd);
    DiracOp::CGflops+=4*f_size;

  } else { 

    vaxpy3((Float *)frm_tmp1,
	   &minus_kappa_sq,
	  (Float *)frm_tmp1,
	  (Float *)in,f_size/6);
    DiracOp::CGflops+=2*f_size;

  }

//------------------------------------------------------------------
// Apply M^dag out <- frm_tmp1
//------------------------------------------------------------------
  dwf_mdag(out, gauge_field, frm_tmp1, mass, dwf_lib_arg);

}


CPS_END_NAMESPACE
