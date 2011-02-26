#ifdef USE_SSE
//#include "dwf_dslash_5_plus-0411.C"
#include "dwf_dslash_5_plus-nonowait.C"
//#include "dwf_dslash_5_plus-orig.C"
#else

#include <unistd.h>
#include<config.h>
#include <pmmintrin.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/sse/dwf_dslash_5_plus.C,v $
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
#include<util/smalloc.h>
#include<util/qblas_extend.h>
#include<comms/scu.h>
CPS_START_NAMESPACE

void dwf_dslash_5_plus_dag0(Vector *out, 
		       Vector *in, 
		       Float mass,
		       Dwf *dwf_lib_arg)
{
#pragma omp parallel default(shared)
  {
    // Initializations
    //------------------------------------------------------------------
    
    int idx;
    IFloat *f_in;
    IFloat *f_out;
    int x;
    int s;
 
    const IFloat two_over_a5 = 2.0 * GJP.DwfA5Inv();
    const IFloat neg_mass_two_over_a5 = -2.0 * mass * GJP.DwfA5Inv();
    const int local_ls    = GJP.SnodeSites(); 
    const int s_nodes     = GJP.Snodes();
    const int s_node_coor = GJP.SnodeCoor();
    const int vol_4d_cb   = dwf_lib_arg->vol_4d / 2;
    const int max_dex((local_ls-1)*vol_4d_cb);
    const int ls_stride   = 24 * vol_4d_cb;
    
    IFloat *comm_buf   = dwf_lib_arg->comm_buf;
    
    // [1 + gamma_5] term (if dag=1 [1 - gamma_5] term)
    //
    // out[s] = [1 + gamma_5] in[s-1]
    //------------------------------------------------------------------        
    f_in  = (IFloat *) in;
    f_out = (IFloat *) out + ls_stride; 
       
    register __m128d al;
    al = _mm_loaddup_pd(&two_over_a5);
    register __m128d nal;
    nal = _mm_loaddup_pd(&neg_mass_two_over_a5);

    register __m128d x0, x1, x2, x3, x4, x5;		\
    register __m128d y0, y1, y2, y3, y4, y5;		\
    register __m128d z0, z1;				\


#define DAXPY( _A, _X, _Y, _x, _y )		\
    _x = _mm_load_pd( _X );			\
    _y = _mm_load_pd( _Y );			\
    _x = _mm_add_pd( _y, _mm_mul_pd(_A, _x) );	\
    _mm_store_pd( _Y, _x );			\

#define DIST_XY 6
#define DAXPY12( _A, _X, _Y )			\
    DAXPY(_A, _X+0,  _Y+0,  x0, y0);		\
    DAXPY(_A, _X+2,  _Y+2,  x1, y1);		\
    DAXPY(_A, _X+4,  _Y+4,  x2, y2);		\
    DAXPY(_A, _X+6,  _Y+6,  x3, y3);		\
    DAXPY(_A, _X+8,  _Y+8,  x4, y4);		\
    DAXPY(_A, _X+10, _Y+10, x5, y5);		\
    								\
    _mm_prefetch((char *)(_X + DIST_XY*24), _MM_HINT_T0);	\
    _mm_prefetch((char *)(_Y + DIST_XY*24), _MM_HINT_T0);	\


#pragma omp  for schedule(static)
    for (idx=0;idx<max_dex;idx++)
      {
	//	cblas_daxpy(12,two_over_a5,f_in+24*idx,f_out+24*idx);
	DAXPY12( al, f_in+24*idx, f_out+24*idx );
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
 
    f_in  = (IFloat *) in+ (local_ls-1)*ls_stride; 
    f_out = (IFloat *) out;
    
    if(s_node_coor == 0) {
#pragma omp for schedule(static)
    for(x=0; x<vol_4d_cb; x++)
      {
	const int shift(24*x);
	//IFloat *f_temp(f_in+shift);
	//cblas_daxpy(12,neg_mass_two_over_a5,f_temp,f_out+shift);
	DAXPY12( nal, f_in+shift, f_out+shift );

      }
    } else {
#pragma omp for schedule(static)
    for(x=0; x<vol_4d_cb; x++)
      {
	const int shift(24*x);
	//IFloat *f_temp(f_in+shift);
	//cblas_daxpy(12,two_over_a5,f_temp,f_out+shift);
	DAXPY12( al, f_in+shift, f_out+shift);
      }

    }

    
    // [1 - gamma_5] term (if dag=1 [1 + gamma_5] term)
    // 
    // out[s] = [1 - gamma_5] in[s+1]
    //------------------------------------------------------------------
    f_in  = (IFloat *) in + 12 + ls_stride;
    f_out = (IFloat *) out + 12;
    
#pragma omp for  schedule(static)
    for (idx=0;idx<max_dex;idx++)
      {
	int shift(24*idx);
	//cblas_daxpy(12,two_over_a5,f_in+shift,f_out+shift);
	DAXPY12( al, f_in+shift, f_out+shift); 
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

    f_in  = (IFloat *) in +12;
    f_out = (IFloat *) out+12+ (local_ls-1)*ls_stride;

    if(s_node_coor == s_nodes - 1) { 
#pragma omp for schedule(static)
      for(x=0; x<vol_4d_cb; x++){
	const int shift(24*x);
	//IFloat *f_temp (f_in+shift);
	//cblas_daxpy(12,neg_mass_two_over_a5,f_temp,f_out+shift);
	DAXPY12( nal, f_in+shift, f_out+shift );
      }
    } else {
#pragma omp for schedule(static)
      for(x=0; x<vol_4d_cb; x++){
	const int shift(24*x);
	//IFloat *f_temp (f_in+shift);
	//cblas_daxpy(12,two_over_a5,f_temp,f_out+shift);
	DAXPY12( al, f_in+shift, f_out+shift );
      }
    }
  }
    
  const int local_ls    = GJP.SnodeSites(); 
  const int vol_4d_cb   = dwf_lib_arg->vol_4d / 2;
  DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}

void dwf_dslash_5_plus_dag1(Vector *out, 
		       Vector *in, 
		       Float mass,
		       Dwf *dwf_lib_arg)
{

#pragma omp parallel default(shared)
  {
    // Initializations
    //------------------------------------------------------------------
    
    int idx;
    IFloat *f_in;
    IFloat *f_out;
    int x;
    int s;

    const IFloat two_over_a5 = 2.0 * GJP.DwfA5Inv();
    const IFloat neg_mass_two_over_a5 = -2.0 * mass * GJP.DwfA5Inv();
    const int local_ls    = GJP.SnodeSites(); 
    const int s_nodes     = GJP.Snodes();
    const int s_node_coor = GJP.SnodeCoor();
    const int vol_4d_cb   = dwf_lib_arg->vol_4d / 2;
    const int max_dex((local_ls-1)*vol_4d_cb);
    const int ls_stride   = 24 * vol_4d_cb;
  
    IFloat *comm_buf   = dwf_lib_arg->comm_buf;
    
    // [1 + gamma_5] term (if dag=1 [1 - gamma_5] term)
    //
    // out[s] = [1 + gamma_5] in[s-1]
    //------------------------------------------------------------------        
    f_in  = (IFloat *) in;
    f_out = (IFloat *) out;
    //if(dag == 1){
      f_in  =  f_in + 12;
      f_out = f_out + 12;
      //}
    f_out = f_out + ls_stride; 
       
#pragma omp  for schedule(static)
    for (idx=0;idx<max_dex;idx++)
      {
	cblas_daxpy(12,two_over_a5,f_in+24*idx,f_out+24*idx);
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
    
    //if(dag == 1){
      f_in  =  f_in + 12;
      f_out = f_out + 12;
      //}

#pragma omp for schedule(static)
    for(x=0; x<vol_4d_cb; x++)
      {
	int shift(24*x);
	IFloat *f_temp(f_in+shift);
	
	if (s_nodes != 1 ) {
	  f_temp = comm_buf;
	  getMinusData(f_temp, f_in+shift, 12, 4);
	}
	
	if(s_node_coor == 0) {
	  cblas_daxpy(12,neg_mass_two_over_a5,f_temp,f_out+shift);
	}
	else {
	  cblas_daxpy(12,two_over_a5,f_temp,f_out+shift);
	}
      }
    
    
    // [1 - gamma_5] term (if dag=1 [1 + gamma_5] term)
    // 
    // out[s] = [1 - gamma_5] in[s+1]
    //------------------------------------------------------------------
    f_in  = (IFloat *) in;
    f_out = (IFloat *) out;
#if 0
    if(dag == 0){
      f_in  =  f_in + 12;
      f_out = f_out + 12;
    }
#endif
    f_in = f_in + ls_stride;
    
#pragma omp for schedule(static)
    for (idx=0;idx<max_dex;idx++)
      {
	const int shift(24*idx);
	cblas_daxpy(12,two_over_a5,f_in+shift,f_out+shift);
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
    
#if 0
    if(dag == 0){
      f_in  =  f_in + 12;
      f_out = f_out + 12;
    }
#endif    
    f_out = f_out + (local_ls-1)*ls_stride;

#pragma omp for schedule(static)
    for(x=0; x<vol_4d_cb; x++){
      const int shift(24*x);
      IFloat *f_temp (f_in+shift);
      
      if (s_nodes != 1 ) {
	f_temp = comm_buf;
	getPlusData(f_temp, f_in+shift, 12, 4);
      }
      
      if(s_node_coor == s_nodes - 1) { 
	cblas_daxpy(12,neg_mass_two_over_a5,f_temp,f_out+shift);
      }
      else {
	cblas_daxpy(12,two_over_a5,f_temp,f_out+shift);
      }
    }
  } // omp parallel

  const int local_ls    = GJP.SnodeSites(); 
  const int vol_4d_cb   = dwf_lib_arg->vol_4d / 2;

    DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
}


void dwf_dslash_5_plus(Vector *out, 
		       Vector *in, 
		       Float mass,
		       int dag, 
		       Dwf *dwf_lib_arg)
{
  if (dag == 0) 
    dwf_dslash_5_plus_dag0(out, in, mass, dwf_lib_arg);
  else
    dwf_dslash_5_plus_dag1(out, in, mass, dwf_lib_arg);
}

CPS_END_NAMESPACE
#endif
