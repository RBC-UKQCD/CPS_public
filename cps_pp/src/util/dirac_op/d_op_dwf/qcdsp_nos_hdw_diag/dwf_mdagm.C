#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos_hdw_diag/dwf_mdagm.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: dwf_mdagm.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:12:43  anj
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
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos_hdw_diag/dwf_mdagm.C,v $
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
CPS_START_NAMESPACE


void dwf_mdagm(Vector *out, 
	       Matrix *gauge_field, 
	       Vector *in, 
	       Float *dot_prd,
	       Float mass, 
	       Dwf *dwf_lib_arg)
{
  Vector *frm_tmp1 = (Vector *) dwf_lib_arg->frm_tmp1;

//------------------------------------------------------------------
// Apply M
//------------------------------------------------------------------
  dwf_m(frm_tmp1, gauge_field, in, mass, dwf_lib_arg);

//------------------------------------------------------------------
// Calculate the dot product <out, in> = <M in, M in>
//------------------------------------------------------------------
  if(dot_prd != 0){
    int f_size = 24 * dwf_lib_arg->vol_4d * dwf_lib_arg->ls / 2;
    *dot_prd = frm_tmp1->NormSqNode(f_size);
  }

//------------------------------------------------------------------
// Apply M^dag
//------------------------------------------------------------------
  dwf_mdag(out, gauge_field, frm_tmp1, mass, dwf_lib_arg);

}


CPS_END_NAMESPACE
