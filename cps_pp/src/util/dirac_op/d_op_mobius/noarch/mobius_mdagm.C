#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_mdagm.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// mobius_mdagm.C
//
// mobius_mdagm M^dag M where M is the fermion matrix.
// The in, out fields are defined on the checkerboard lattice.
// <out, in> = <mobius_mdagm*in, in> is returned in dot_prd.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/mobius.h>
#include<util/gjp.h>
#include<util/dirac_op.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE


void mobius_mdagm(Vector *out, 
		  Matrix *gauge_field, 
		  Vector *in, 
		  Float *dot_prd,
		  Float mass, 
		  Dwf *mobius_lib_arg)
{
  Vector *frm_tmp1 = (Vector *) mobius_lib_arg->frm_tmp1;
  
//------------------------------------------------------------------
// Apply M
//------------------------------------------------------------------
  mobius_m(frm_tmp1, gauge_field, in, mass, mobius_lib_arg);

//------------------------------------------------------------------
// Calculate the dot product <out, in> = <M in, M in>
//------------------------------------------------------------------
  if(dot_prd != 0){
    int f_size = 24 * mobius_lib_arg->vol_4d * mobius_lib_arg->ls / 2;
    *dot_prd = frm_tmp1->NormSqNode(f_size);
    DiracOp::CGflops+=2*f_size;
  }

//------------------------------------------------------------------
// Apply M^dag
//------------------------------------------------------------------
  mobius_mdag(out, gauge_field, frm_tmp1, mass, mobius_lib_arg);
}


//
//  H = Gm5 M
//  (H-mu)(H-mu) =   (Gm5 M-mu)(Gm5 M-mu) = (M^dag -mu Gm5)  (M -mu Gm5)
// = M^dag M - 2 mu Gm5 M  + mu^2
//
// perhaps bit faster in the latter form as only issue Gm5 onece
// although dot product would be slower, but dotproduct won't be used too much.
//
//  1. tmp <- M in
//  2. out <-  M^dag  tmp
//  3. tmp <-  Gm5 tmp
//  4. out -= 2 mu tmp
//  5. out += mu^2 in
//

// this might be broken!
void mobius_mdagm_shift(Vector *out, 
			Matrix *gauge_field, 
			Vector *in, 
			Float *dot_prd,
			Float mass, 
			Dwf *mobius_lib_arg,
			Float mu)
{
  Vector *frm_tmp1 = (Vector *) mobius_lib_arg->frm_tmp1;
  Vector *frm_tmp2 = (Vector *) mobius_lib_arg->frm_tmp2;

  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int ls = mobius_lib_arg->ls;
  const int f_size = vol_4d_cb * ls * 24;
  

  //  1. tmp1 <- M in
  //------------------------------------------------------------------
  // Apply M
  //------------------------------------------------------------------
  mobius_m(frm_tmp1, gauge_field, in, mass, mobius_lib_arg);

  //  2. out <-  M^dag  tmp1
  //------------------------------------------------------------------
  // Apply M^dag
  //------------------------------------------------------------------
  mobius_mdag(out, gauge_field, frm_tmp1, mass, mobius_lib_arg);
  
  //  3. tmp2 <-  Gm5 tmp1
  ReflectAndMultGamma5( frm_tmp2, frm_tmp1, vol_4d_cb, ls);
    
  //  4. out += -2 mu tmp2
#ifndef USE_BLAS
  fTimesV1PlusV2( (IFloat*)out, -2*mu, (IFloat*)frm_tmp2, (IFloat*)out, f_size);
#else
  cblas_daxpy(  (IFloat*)out, -2*mu, (IFloat*)frm_tmp2, (IFloat*)out, f_size);
#endif
  
//  5. out += mu^2 in

#ifndef USE_BLAS
  fTimesV1PlusV2( (IFloat*)out, mu*mu, (IFloat*)in, (IFloat*)out, f_size);
#else
  cblas_daxpy(  (IFloat*)out, mu*mu, (IFloat*)in, (IFloat*)out, f_size);
#endif

  
  //------------------------------------------------------------------
  // Calculate the dot product <out, in> 
  //------------------------------------------------------------------
  if(dot_prd != 0){
    int f_size = 24 * mobius_lib_arg->vol_4d * mobius_lib_arg->ls / 2;
    *dot_prd = frm_tmp1->NormSqNode(f_size);
    DiracOp::CGflops+=2*f_size;
  }

}

CPS_END_NAMESPACE
