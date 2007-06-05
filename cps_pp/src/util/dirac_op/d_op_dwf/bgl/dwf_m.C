#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/bgl/dwf_m.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// dwf_m.C
//
// dwf_m is the fermion matrix.  
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
CPS_START_NAMESPACE


void  dwf_m(Vector *out, 
	    Matrix *gauge_field, 
	    Vector *in, 
	    Float mass, 
	    Dwf *dwf_lib_arg)
{

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------
  int f_size = 24 * dwf_lib_arg->vol_4d * dwf_lib_arg->ls / 2;
  Float minus_kappa_sq = -dwf_lib_arg->dwf_kappa * dwf_lib_arg->dwf_kappa;
  Vector *frm_tmp2 = (Vector *) dwf_lib_arg->frm_tmp2;
  
//------------------------------------------------------------------
// Apply Dslash E <- O
//------------------------------------------------------------------
  sync();
  VRB.Func("dwf_m","dwf_dslash_4");
//  dwf_dslash(frm_tmp2, gauge_field, in, mass, 1, 0, dwf_lib_arg);
  dwf_dslash_4(frm_tmp2, gauge_field, in, 1, 0, dwf_lib_arg);
  sync();
  VRB.Func("dwf_m","dwf_dslash_5_plus");
  dwf_dslash_5_plus(frm_tmp2, in, mass, 0, dwf_lib_arg);

//------------------------------------------------------------------
// Apply Dslash O <- E
//------------------------------------------------------------------
  sync();
  VRB.Func("dwf_m","dwf_dslash_4");
//  dwf_dslash(out, gauge_field, frm_tmp2, mass, 0, 0, dwf_lib_arg);
  dwf_dslash_4(out, gauge_field, frm_tmp2, 0, 0, dwf_lib_arg);
  sync();
  VRB.Func("dwf_m","dwf_dslash_5_plus");
  dwf_dslash_5_plus(out, frm_tmp2, mass, 0, dwf_lib_arg);

  sync();
  
//------------------------------------------------------------------
// out = in - dwf_kappa_sq * out
//------------------------------------------------------------------
  out->FTimesV1PlusV2(minus_kappa_sq, out, in, f_size); 
  DiracOp::CGflops+=2*f_size;

}






CPS_END_NAMESPACE
