#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:06 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdsp_nos_hdw_diag/dwf_mdagm.C,v 1.2 2004-06-04 21:14:06 chulwoo Exp $
//  $Id: dwf_mdagm.C,v 1.2 2004-06-04 21:14:06 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.2 $
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
