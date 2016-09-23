#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:57:50 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/qcdoc/dwf_mdagm.C,v 1.6 2004/08/18 11:57:50 zs Exp $
//  $Id: dwf_mdagm.C,v 1.6 2004/08/18 11:57:50 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: dwf_mdagm.C,v $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/qcdoc/dwf_mdagm.C,v $
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
//#include "dwf_internal.h"
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

    vaxpy3_norm(frm_tmp1, &minus_kappa_sq, frm_tmp1, in,f_size/6,dot_prd);
    DiracOp::CGflops+=4*f_size;

  } else { 

    vaxpy3(frm_tmp1, &minus_kappa_sq, frm_tmp1, in,f_size/6);
    DiracOp::CGflops+=2*f_size;

  }


//------------------------------------------------------------------
// Apply M^dag out <- frm_tmp1
//------------------------------------------------------------------
  dwf_mdag(out, gauge_field, frm_tmp1, mass, dwf_lib_arg);

}


CPS_END_NAMESPACE
