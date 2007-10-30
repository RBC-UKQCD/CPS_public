#include<config.h>
#include<stdlib.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/bgl/dwf_mdag.C,v $
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
#include<util/dirac_op.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
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
//  out->FTimesV1PlusV2(minus_kappa_sq, out, in, f_size); 
  Float *out_f = (Float *)out;
  Float *in_f = (Float *)in;
  xaxpy(&minus_kappa_sq,out_f,in_f,f_size/6);
  DiracOp::CGflops+=2*f_size;

}
CPS_END_NAMESPACE
