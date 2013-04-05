#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_mdag.C,v $
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
// ( M_5 - kappa_b^2 M4eo M_5^-1 M4oe )^ dag
void  mobius_mdag(Vector *out, 
		  Matrix *gauge_field, 
		  Vector *in, 
		  Float mass, 
		  Dwf *mobius_lib_arg)
{
  

  //------------------------------------------------------------------
  // Initializations
  //------------------------------------------------------------------
  IFloat *g_field;
  Wilson *wilson_p;
  int size_cb[2];

  const int ls = mobius_lib_arg->ls;
  const int f_size = 24 * mobius_lib_arg->vol_4d * ls / 2;
  const Float minus_kappa_b_sq = -mobius_lib_arg->mobius_kappa_b * mobius_lib_arg->mobius_kappa_b;
  const Float kappa_ratio = mobius_lib_arg->mobius_kappa_b/mobius_lib_arg->mobius_kappa_c;
  const Float b_coeff =  GJP.Mobius_b();
  const Float c_coeff =  GJP.Mobius_c();
  Vector  *frm_tmp2 = (Vector *) mobius_lib_arg->frm_tmp2;
  Vector  *frm_tmp3 = (Vector *) mobius_lib_arg->frm_tmp3;
  Float norm;
  IFloat* frm_in;
  IFloat* frm_out;
  g_field = (IFloat *) gauge_field;
  wilson_p = mobius_lib_arg->wilson_p;
  size_cb[0] = 24*wilson_p->vol[0];
  size_cb[1] = 24*wilson_p->vol[1];
  
  //  out = [ 1 + kappa_b/kappa_c 1/2 dslash_5  - kappa_b^2 Meo M5inv Moe] in
  //    1. ftmp2 = Meo M5inv Moe in
  //    2. out <-  in
  //    3. out += -kappa_b^2 ftmp2
  //    4. out +=  -kappa_b/kappa_c /2 dslash_5 in
  //         (done by the dslash_5 with a5_inv = -kappa_b/kappa_c/2 *GJP.MobiusA5Inv() ) 
  time_elapse();
  // O <- E
  int parity = 0;
  // out = D^dag_W * frm_in
  frm_in = (IFloat*)in;
  frm_out = (IFloat*)out;
  for(int i=0; i<ls; i++){

    // Apply on 4-dim "parity" checkerboard part
    //------------------------------------------------------------
    wilson_dslash(frm_out, g_field, frm_in, parity, 1, wilson_p);
    
    frm_in = frm_in + size_cb[parity];
    frm_out = frm_out + size_cb[parity];
  }
  DEBUG_MOBIUS_DSLASH("mobius_dslash_4_dag %e\n", time_elapse());

  //TIZB  "b_coeff" part
  moveFloat((IFloat*)frm_tmp2, (IFloat*)out, f_size);
  vecTimesEquFloat( (IFloat*)out, b_coeff, f_size);
  DEBUG_MOBIUS_DSLASH("mobius_dslash_5_dag_b %e\n", time_elapse());
  mobius_kappa_dslash_5_plus(out, frm_tmp2, mass, 1, mobius_lib_arg, c_coeff);
  DEBUG_MOBIUS_DSLASH("mobius_dslash_5_dag_c %e\n", time_elapse());
  mobius_m5inv(out, mass, 1, mobius_lib_arg);
  DEBUG_MOBIUS_DSLASH("mobius_m5inv_dag %e\n", time_elapse());

  // E <- O
  parity = 1;
  frm_in = (IFloat*)out;
  frm_out = (IFloat*)frm_tmp2;
  for(int i=0; i<ls; i++){

    // Apply on 4-dim "parity" checkerboard part
    //------------------------------------------------------------
    wilson_dslash(frm_out, g_field, frm_in, parity, 1, wilson_p);
    
    frm_in = frm_in + size_cb[parity];
    frm_out = frm_out + size_cb[parity];
  }
  DEBUG_MOBIUS_DSLASH("mobius_dslash_4_dag %e\n", time_elapse());

  moveFloat((IFloat*)frm_tmp3, (IFloat*)frm_tmp2, f_size);
  DEBUG_MOBIUS_DSLASH("mobius dag moveFloat %e\n", time_elapse());
  vecTimesEquFloat( (IFloat*)frm_tmp3, b_coeff, f_size);
  DEBUG_MOBIUS_DSLASH("mobius dag times b %e\n", time_elapse());

  mobius_kappa_dslash_5_plus(frm_tmp3, frm_tmp2, mass, 1, mobius_lib_arg, c_coeff);
  DEBUG_MOBIUS_DSLASH("mobius_dslash_5_dag_c_plus %e\n", time_elapse());
  

  // M_5 part, add
  //------------------------------------------------------------------
  //    2. out <-  in
  //------------------------------------------------------------------
#ifndef USE_BLAS
  moveFloat((IFloat*)out, (IFloat*)in, f_size);
#else
  cblas_dcopy(f_size, (IFloat*)in, 1, (IFloat*)out, 1);
#endif
  DEBUG_MOBIUS_DSLASH("mobius dag out<-in %e\n", time_elapse());
  
  //------------------------------------------------------------------
  //    3. out += -kap2 ftmp3
  //------------------------------------------------------------------
#ifndef USE_BLAS
  fTimesV1PlusV2((IFloat*)out, minus_kappa_b_sq, (IFloat*)frm_tmp3,
		 (IFloat *)out, f_size);
#else
  cblas_daxpy(f_size, minus_kappa_b_sq, (IFloat*)frm_tmp3, 1, 
	      (IFloat *)out,  1);
#endif
  DEBUG_MOBIUS_DSLASH("mobius dag += -kap2 ftmp3 %e\n", time_elapse());
  
  //------------------------------------------------------------------
  //    4. out +=  kap  dslash_5 in
  //------------------------------------------------------------------
  mobius_kappa_dslash_5_plus(out, in, mass, 1, mobius_lib_arg, kappa_ratio);
  DEBUG_MOBIUS_DSLASH("mobius_kappa_dslash_5_dag_plus %e\n", time_elapse());

  // Flops count in this function is two AXPY = 4 flops per vector elements
  DiracOp::CGflops +=  3*f_size; 

}






CPS_END_NAMESPACE
