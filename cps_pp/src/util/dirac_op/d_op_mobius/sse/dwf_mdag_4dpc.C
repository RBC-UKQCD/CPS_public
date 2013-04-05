#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/sse/dwf_mdag_4dpc.C,v $
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
#include<util/time_cps.h>

#include "blas-subs.h"

CPS_START_NAMESPACE


// added  the last argument, a_five_inv to change the all over coeff
void dwf_dslash_5_plus_a_five(Vector *out, 
			      Vector *in,
			      Float mass,
			      int dag,
			      Dwf *dwf_lib_arg,
			      Float a_five_inv );



//4d precond. dwf Dirac op:
// ( M_5 - kappa^2 M4eo M_5^-1 M4oe )^ dag
void  dwf_mdag_4dpc(Vector *out, 
		    Matrix *gauge_field, 
		    Vector *in, 
		    Float mass, 
		    Dwf *dwf_lib_arg)
{


  //------------------------------------------------------------------
  // Initializations
  //------------------------------------------------------------------
  const int f_size = 24 * dwf_lib_arg->vol_4d * dwf_lib_arg->ls / 2;
  const Float minus_kappa = -dwf_lib_arg->dwf_kappa;
  const Float minus_kappa_sq = -dwf_lib_arg->dwf_kappa * dwf_lib_arg->dwf_kappa;
  Vector  *frm_tmp2 = (Vector *) dwf_lib_arg->frm_tmp2;
  //Vector *temp = (Vector *) smalloc(f_size * sizeof(Float));
  Float norm;

  
  //  out = [ 1 - kap dslash_5  - kap2 Meo M5inv Moe] in
  //
  //    1. ftmp2 = Meo M5inv Moe in
  //    2. out <-  in
  //    3. out += -kap2 ftmp2
  //    4. out +=  -kap  dslash_5 in
  //         (done by the dslash_5 with a5_inv = -kap *GJP.DwfA5Inv() ) 


  //--------------------------------------------------------------
  //    1. ftmp2 = Meo M5inv Moe in
  //--------------------------------------------------------------
  // Apply Dslash O <- E
  //------------------------------------------------------------------
  time_elapse();
  //dwf_dslash_4(frm_tmp2, gauge_field, in, 1, 1, dwf_lib_arg);
  dwf_dslash_4(out, gauge_field, in, 0, 1, dwf_lib_arg);
  DEBUG_DWF_DSLASH("dag1 dwf_dslash_4(o<-e) %e\n", time_elapse());

  //------------------------------------------------------------------
  // Apply M_5^-1 (hopping in 5th dir + diagonal)
  //------------------------------------------------------------------
  dwf_m5inv(out, mass, 1, dwf_lib_arg);
  DEBUG_DWF_DSLASH("dag1 dwf_m5inv %e\n", time_elapse());
  
  //------------------------------------------------------------------
  // Apply Dslash E <- O
  //------------------------------------------------------------------
  //dwf_dslash_4(out, gauge_field, frm_tmp2, 0, 0, dwf_lib_arg);
  dwf_dslash_4(frm_tmp2, gauge_field, out, 1, 1, dwf_lib_arg);
  DEBUG_DWF_DSLASH("dag1 dwf_dslash_4(e<-o) %e\n", time_elapse());
  
  //------------------------------------------------------------------
  //    2. out <-  in
  //------------------------------------------------------------------
#ifndef USE_BLAS
  moveFloat((IFloat*)out, (IFloat*)in, f_size);
#else
    cblas_dcopy(f_size, (IFloat*)in, 1, (IFloat*)out, 1);
#endif
  DEBUG_DWF_DSLASH("out <- in %e\n", time_elapse());
  
  //------------------------------------------------------------------
  //    3. out += -kap2 ftmp2
  //------------------------------------------------------------------
#ifndef USE_BLAS
  fTimesV1PlusV2((IFloat*)out, minus_kappa_sq, (IFloat*)frm_tmp2,
		 (IFloat *)out, f_size);
#else
  //    #define NN 24
  //    for(int x=0; x<f_size/NN; x++)
  //    cblas_daxpy(NN, minus_kappa_sq, (IFloat*)frm_tmp2+x*NN,1, 
  // 		(IFloat *)out+x*NN,  1);

  cblas_daxpy(f_size, minus_kappa_sq, (IFloat*)frm_tmp2,1, 
            (IFloat *)out,  1);
#endif
  
  DEBUG_DWF_DSLASH("dag1 out+= kap2 frm_tmp2 %e\n", time_elapse());

  //------------------------------------------------------------------
  //    4. out +=  -kap  dslash_5 in
  //         (done by the dslash_5 with a5_inv = -kap *GJP.DwfA5Inv() ) 
  //------------------------------------------------------------------
  //------------------------------------------------------------------
  // Apply Dslash_5 to in,  with  minus_kappa as a all over constants
  //------------------------------------------------------------------
  dwf_dslash_5_plus_a_five(out, in, mass, 1, dwf_lib_arg, GJP.DwfA5Inv()*minus_kappa );
  DEBUG_DWF_DSLASH("dag1 dwf_dslash_5_plus %e\n", time_elapse());

  // Flops count in this function is two AXPY = 4 flops per vector elements
  DiracOp::CGflops +=  3*f_size; 

}






CPS_END_NAMESPACE
