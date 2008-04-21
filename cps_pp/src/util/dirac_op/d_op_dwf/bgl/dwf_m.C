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
#include<util/time_cps.h>
#include<util/dirac_op.h>
CPS_START_NAMESPACE


void  dwf_m(Vector *out, 
	    Matrix *gauge_field, 
	    Vector *in, 
	    Float mass, 
	    Dwf *dwf_lib_arg)
{

  static int called=0;
  called++;
  Float dslash_4_time=0.;
  Float dslash_5_time=0.;

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------
  int f_size = 24 * dwf_lib_arg->vol_4d * dwf_lib_arg->ls / 2;
  Float minus_kappa_sq = -dwf_lib_arg->dwf_kappa * dwf_lib_arg->dwf_kappa;
  Vector *frm_tmp2 = (Vector *) dwf_lib_arg->frm_tmp2;
  
//------------------------------------------------------------------
// Apply Dslash E <- O
//------------------------------------------------------------------
//  sync();
  dwf_dslash(frm_tmp2, gauge_field, in, mass, 1, 0, dwf_lib_arg);
#if 0
  dslash_4_time -=dclock();
  dwf_dslash_4(frm_tmp2, gauge_field, in, 1, 0, dwf_lib_arg);
  dslash_4_time +=dclock();
//  sync();
  dslash_5_time -=dclock();
  dwf_dslash_5_plus(frm_tmp2, in, mass, 0, dwf_lib_arg);
  dslash_5_time +=dclock();
#endif

//------------------------------------------------------------------
// Apply Dslash O <- E
//------------------------------------------------------------------
//  sync();
  dwf_dslash(out, gauge_field, frm_tmp2, mass, 0, 0, dwf_lib_arg);
#if 0
  dslash_4_time -=dclock();
  dwf_dslash_4(out, gauge_field, frm_tmp2, 0, 0, dwf_lib_arg);
  dslash_4_time +=dclock();
//  sync();
  dslash_5_time -=dclock();
  dwf_dslash_5_plus(out, frm_tmp2, mass, 0, dwf_lib_arg);
  dslash_5_time +=dclock();
#endif

//  sync();
  
//------------------------------------------------------------------
// out = in - dwf_kappa_sq * out
//------------------------------------------------------------------
//  out->FTimesV1PlusV2(minus_kappa_sq, out, in, f_size); 
  Float *out_f = (Float *)out;
  Float *in_f = (Float *)in;
  xaxpy(&minus_kappa_sq,out_f,in_f,f_size/6);
  DiracOp::CGflops+=2*f_size;
#if 0
  if(called%1000==0){
    print_time("dwf_m","dslash_4_time",dslash_4_time/1000.);
    print_time("dwf_m","dslash_5_time",dslash_5_time/1000.);
    dslash_4_time=0.;
    dslash_5_time=0.;
  }
#endif

}






CPS_END_NAMESPACE
