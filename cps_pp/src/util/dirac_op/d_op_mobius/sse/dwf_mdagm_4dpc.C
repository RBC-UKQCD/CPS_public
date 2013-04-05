#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/sse/dwf_mdagm_4dpc.C,v $
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
#include<util/dirac_op.h>
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
  dwf_m_4dpc(frm_tmp1, gauge_field, in, mass, dwf_lib_arg);

//------------------------------------------------------------------
// Calculate the dot product <out, in> = <M in, M in>
//------------------------------------------------------------------
  if(dot_prd != 0){
    int f_size = 24 * dwf_lib_arg->vol_4d * dwf_lib_arg->ls / 2;
    *dot_prd = frm_tmp1->NormSqNode(f_size);
    DiracOp::CGflops+=2*f_size;
  }

//------------------------------------------------------------------
// Apply M^dag
//------------------------------------------------------------------
  dwf_mdag_4dpc(out, gauge_field, frm_tmp1, mass, dwf_lib_arg);

}


// FIXME: Using SSE, how to flip the sign ?

void ReflectAndMultGamma5( Vector *out, const Vector *in,  int nodevol, int ls)
{
  for(int s=0; s< ls; ++s) { 
    IFloat *p = (IFloat *)out + 24*nodevol*s;
    IFloat *q = (IFloat *)in + 24*nodevol*(ls-1-s);
    for(int n = 0; n < nodevol; ++n)
      {
	int i;
	for(i = 0; i < 12; ++i)
	  *p++ = *q++;
	for(i = 0; i < 12; ++i)
	  *p++ = - *q++;
      }
  }

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
void dwf_mdagm_shift(Vector *out, 
		     Matrix *gauge_field, 
		     Vector *in, 
		     Float *dot_prd,
		     Float mass, 
		     Dwf *dwf_lib_arg,
		     Float mu)
{
  Vector *frm_tmp1 = (Vector *) dwf_lib_arg->frm_tmp1;
  Vector *frm_tmp2 = (Vector *) dwf_lib_arg->frm_tmp2;

  const int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  const int ls = dwf_lib_arg->ls;
  const int f_size = vol_4d_cb * ls * 24;
  //const Float mu = dwf_lib_arg->eigen_shift;
  

  //printf("shift=%e\n",mu);exit(1);
  //  1. tmp1 <- M in
  //------------------------------------------------------------------
  // Apply M
  //------------------------------------------------------------------
  dwf_m_4dpc(frm_tmp1, gauge_field, in, mass, dwf_lib_arg);

  //  2. out <-  M^dag  tmp1
  //------------------------------------------------------------------
  // Apply M^dag
  //------------------------------------------------------------------
  dwf_mdag_4dpc(out, gauge_field, frm_tmp1, mass, dwf_lib_arg);
  
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
    int f_size = 24 * dwf_lib_arg->vol_4d * dwf_lib_arg->ls / 2;
    *dot_prd = frm_tmp1->NormSqNode(f_size);
    DiracOp::CGflops+=2*f_size;
  }


}

CPS_END_NAMESPACE
