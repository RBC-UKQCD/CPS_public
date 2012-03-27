#ifdef USE_OMP_DWF_FORCE
#include<omp.h>
#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/noarch/dwf_dslash_5_plus_omp.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf_dslash_5_plus.C
//
// dwf_dslash_5_plus is the derivative part of the 5th direction
// part of the fermion matrix. This routine accumulates the result
// on the out field 
// The in, out fields are defined on the checkerboard lattice.
// The action of this operator is the same for even/odd
// checkerboard fields because there is no gauge field along
// the 5th direction.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//
// Storage order for DWF fermions
//------------------------------------------------------------------
//  
//  |     |
//  | |r| |
//  | |i| |
//  |     |
//  | |r| | = |spin comp|
//  | |i| |
//  |     |
//  | |r| |
//  | |i| |
//  |     |
//  
//  
//  |             |
//  | |spin comp| |
//  |             |
//  |             |
//  | |spin comp| | 
//  |             | = |spinor|
//  |             |
//  | |spin comp| |
//  |             |
//  |             |
//  | |spin comp| |
//  |             |
//  
//  
//  |            |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |  |spinor|  |
//  |     .      | = |s-block|   The spinors are arranged in Wilson
//  |     .      |               order with odd - even 4d-checkerboard
//  |     .      |               storage.
//  |evn/odd vol |
//  |     .      |
//  |  |spinor|  |
//  |            |
//  
//  
//  |                |
//  | |s-block even| |  For even chckerboard
//  | |s-block odd|  |
//  | |s-block even| |
//  | |s-block odd|  |
//  |       .        |
//  |       .        |
//  |       .        |
//  |                |
//
//
//  |                |
//  | |s-block odd|  |  For odd chckerboard
//  | |s-block even| |
//  | |s-block odd|  |
//  | |s-block even| |
//  |       .        |
//  |       .        |
//  |       .        |
//  |                |
//
//------------------------------------------------------------------


CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/gjp.h>
#include<util/dirac_op.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/smalloc.h>
#include<comms/scu.h>
CPS_START_NAMESPACE


inline void fTimesV1PlusV2_single(IFloat *a, IFloat b, const IFloat *c,
const IFloat *d, int len)
{
for(int i = 0; i < len; ++i) {
*a++ = b * *c++ + *d++;
}
}


void dwf_dslash_5_plus(Vector *out, 
		       Vector *in, 
		       Float mass,
		       int dag, 
		       Dwf *dwf_lib_arg)
{
  int x;
  int s;

// Initializations
//------------------------------------------------------------------
  int local_ls = GJP.SnodeSites(); 
  int s_nodes = GJP.Snodes();
  int s_node_coor = GJP.SnodeCoor();
  int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  int ls_stride = 24 * vol_4d_cb;
  IFloat *f_in;
  IFloat *f_out;
//  IFloat *f_temp;
  IFloat *comm_buf = dwf_lib_arg->comm_buf;
  IFloat two_over_a5 = 2.0 * GJP.DwfA5Inv();
  IFloat neg_mass_two_over_a5 = -2.0 * mass * GJP.DwfA5Inv();

// [1 + gamma_5] term (if dag=1 [1 - gamma_5] term)
//
// out[s] = [1 + gamma_5] in[s-1]
//------------------------------------------------------------------
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  if(dag == 1){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  f_out = f_out + ls_stride; 
omp_set_num_threads(64);
#pragma omp parallel for default(shared)
  for(x=0; x<vol_4d_cb*(local_ls-1); x++){
  IFloat *f_in_t = f_in+24*x;
  IFloat *f_out_t = f_out+24*x;

      fTimesV1PlusV2_single(f_out_t, two_over_a5, f_in_t, f_out_t, 12);

  }


// [1 + gamma_5] for lower boundary term (if dag=1 [1 - gamma_5] term)
// If there's only one node along fifth direction, no communication
// is necessary; Otherwise data from adjacent node in minus direction
// will be needed.
// If the lower boundary is the s=0 term
// out[0] = - m_f * [1 + gamma_5] in[ls-1]
// else, out[s] = [1 + gamma_5] in[s-1]
//
//------------------------------------------------------------------

  f_in  = (IFloat *) in;  
  f_in = f_in + (local_ls-1)*ls_stride; 
  f_out = (IFloat *) out;
  
  if(dag == 1){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  
#pragma omp parallel for default(shared)
  for(x=0; x<vol_4d_cb; x++){
  IFloat *f_in_t = f_in+24*x;
  IFloat *f_out_t = f_out+24*x;
    
    IFloat *f_temp = f_in_t;
    
    if(s_node_coor == 0) { 
      fTimesV1PlusV2_single(f_out_t, neg_mass_two_over_a5, f_temp, f_out_t, 12);
    }
    else {
      fTimesV1PlusV2_single(f_out_t, two_over_a5, f_temp, f_out_t, 12);
    }
    
  }


// [1 - gamma_5] term (if dag=1 [1 + gamma_5] term)
// 
// out[s] = [1 - gamma_5] in[s+1]
//------------------------------------------------------------------
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  f_in = f_in + ls_stride;
#pragma omp parallel for default(shared)
  for(x=0; x<vol_4d_cb*(local_ls-1); x++){
  IFloat *f_in_t = f_in+24*x;
  IFloat *f_out_t = f_out+24*x;

      fTimesV1PlusV2_single(f_out_t, two_over_a5, f_in_t, f_out_t, 12);

  }


// [1 - gamma_5] for upper boundary term (if dag=1 [1 + gamma_5] term)
// If there's only one node along fifth direction, no communication
// is necessary; Otherwise data from adjacent node in minus direction
// will be needed.
// If the upper boundary is the s=ls term
// out[ls-1] = - m_f * [1 - gamma_5] in[0]
// else out[s] = [1 - gamma_5] in[s+1]
//
//------------------------------------------------------------------

  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;

  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }

  f_out = f_out + (local_ls-1)*ls_stride;
#pragma omp parallel for default(shared)
  for(x=0; x<vol_4d_cb; x++){
  IFloat *f_in_t = f_in+24*x;
  IFloat *f_out_t = f_out+24*x;

    IFloat *f_temp = f_in_t;
   
    if(s_node_coor == s_nodes - 1) { 
      fTimesV1PlusV2_single(f_out_t, neg_mass_two_over_a5, f_temp, f_out_t, 12);
    }
    else {
      fTimesV1PlusV2_single(f_out_t, two_over_a5, f_temp, f_out_t, 12);
    }
    
  }
  DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;


}



CPS_END_NAMESPACE
#endif
