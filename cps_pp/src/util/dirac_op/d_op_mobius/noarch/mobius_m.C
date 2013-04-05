#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_m.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// mobius_m.C
//
// mobius_m is the fermion matrix.  
// The in, out fields are defined on the checkerboard lattice
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/mobius.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/dirac_op.h>
#include<util/time_cps.h>

#include "blas-subs.h"

CPS_START_NAMESPACE


//4d precond. mobius Dirac op:
// M_5 - kappa_b^2 M4eo M_5^-1 M4oe
void  mobius_m(Vector *out, 
	       Matrix *gauge_field, 
	       Vector *in, 
	       Float mass, 
	       Dwf *mobius_lib_arg)
{

  //------------------------------------------------------------------
  // Initializations
  //------------------------------------------------------------------
  const int f_size = 24 * mobius_lib_arg->vol_4d * mobius_lib_arg->ls / 2;
  const Float kappa_ratio = mobius_lib_arg->mobius_kappa_b/mobius_lib_arg->mobius_kappa_c;
  const Float minus_kappa_b_sq = -mobius_lib_arg->mobius_kappa_b * mobius_lib_arg->mobius_kappa_b;
  Vector  *frm_tmp2 = (Vector *) mobius_lib_arg->frm_tmp2;
  //Vector *temp = (Vector *) smalloc(f_size * sizeof(Float));
  Float norm;

  
  //  out = [ 1 + kappa_b/kappa_c 1/2 dslash_5  - kappa_b^2 Meo M5inv Moe] in
  // (dslash_5 uses (1+-g5), not P_R,L, i.e. no factor of 1/2 which is here out front)
  //    1. ftmp2 = Meo M5inv Moe in
  //    2. out <-  in
  //    3. out += -kappa_b^2 ftmp2
  //    4. out +=  -kappa_b/kappa_c /2 dslash_5 in
  //         (done by the dslash_5 with a5_inv = -kappa_b/kappa_c/2 *GJP.MobiusA5Inv() ) 


  //--------------------------------------------------------------
  //    1. ftmp2 = Meo M5inv Moe in
  //--------------------------------------------------------------
  // Apply Dslash O <- E
  //------------------------------------------------------------------
  time_elapse();
  mobius_dslash_4(out, gauge_field, in, 0, 0, mobius_lib_arg, mass);
  DEBUG_MOBIUS_DSLASH("mobius_dslash_4 %e\n", time_elapse());

  //------------------------------------------------------------------
  // Apply M_5^-1 (hopping in 5th dir + diagonal)
  //------------------------------------------------------------------
  mobius_m5inv(out, mass, 0, mobius_lib_arg);
  DEBUG_MOBIUS_DSLASH("mobius_m5inv %e\n", time_elapse());
  
  //------------------------------------------------------------------
  // Apply Dslash E <- O
  //------------------------------------------------------------------
  mobius_dslash_4(frm_tmp2, gauge_field, out, 1, 0, mobius_lib_arg, mass);
  DEBUG_MOBIUS_DSLASH("mobius_dslash_4 %e\n", time_elapse());
  
  //------------------------------------------------------------------
  //    2. out <-  in
  //------------------------------------------------------------------
#ifndef USE_BLAS
  moveFloat((IFloat*)out, (IFloat*)in, f_size);
#else
  cblas_dcopy(f_size, (IFloat*)in, 1, (IFloat*)out, 1);
#endif
  DEBUG_MOBIUS_DSLASH("out <- in %e\n", time_elapse());
  
  //------------------------------------------------------------------
  //    3. out += -kap2 ftmp2
  //------------------------------------------------------------------
#ifndef USE_BLAS
  fTimesV1PlusV2((IFloat*)out, minus_kappa_b_sq, (IFloat*)frm_tmp2,
		 (IFloat *)out, f_size);
#else

  cblas_daxpy(f_size, minus_kappa_b_sq, (IFloat*)frm_tmp2,1, 
            (IFloat *)out,  1);
#endif
  DEBUG_MOBIUS_DSLASH("mobius out+= kap2 %e\n", time_elapse());
  
  //------------------------------------------------------------------
  //    4. out +=  kappa_b/kappa_c dslash_5 in
  //------------------------------------------------------------------
  mobius_kappa_dslash_5_plus(out, in, mass, 0, mobius_lib_arg, kappa_ratio);
  DEBUG_MOBIUS_DSLASH("mobius_kappa_dslash_5_plus %e\n", time_elapse());

  // Flops count in this function is two AXPY = 4 flops per vector elements
  //DiracOp::CGflops +=  3*f_size; 

}






CPS_END_NAMESPACE
