#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:57:50 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/qcdoc/dwf_mdag.C,v 1.7 2004/08/18 11:57:50 zs Exp $
//  $Id: dwf_mdag.C,v 1.7 2004/08/18 11:57:50 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: dwf_mdag.C,v $
//  $Revision: 1.7 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/qcdoc/dwf_mdag.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// dwf_mdag.C
//
// dwf_mdag is the dagger of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
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


void  dwf_mdag(Vector *out, 
	       Matrix *gauge_field, 
	       Vector *in, 
	       Float mass, 
	       Dwf *dwf_lib_arg)
{

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------
  int f_size = 24 * dwf_lib_arg->vol_4d * dwf_lib_arg->ls / 2;
  Float minus_kappa_sq = - dwf_lib_arg->dwf_kappa * dwf_lib_arg->dwf_kappa;
  Vector *frm_tmp2 = (Vector *) dwf_lib_arg->frm_tmp2;

//------------------------------------------------------------------
// Apply Dslash E <- O
//------------------------------------------------------------------
  dwf_dslash(frm_tmp2, gauge_field, in, mass, 1, 1, dwf_lib_arg);


//------------------------------------------------------------------
// Apply Dslash O <- E
//------------------------------------------------------------------
  dwf_dslash(out, gauge_field, frm_tmp2, mass, 0, 1, dwf_lib_arg);

//------------------------------------------------------------------
// out = in - dwf_kappa_sq * out
//------------------------------------------------------------------
  //out->FTimesV1PlusV2(minus_kappa_sq, out, in, f_size); 
  vaxpy3(out,&minus_kappa_sq,out,in,f_size/6);
  DiracOp::CGflops+=2*f_size;
}
CPS_END_NAMESPACE
