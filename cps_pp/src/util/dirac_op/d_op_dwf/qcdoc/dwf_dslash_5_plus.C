#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004/10/14 22:03:48 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/qcdoc/dwf_dslash_5_plus.C,v 1.7 2004/10/14 22:03:48 chulwoo Exp $
//  $Id: dwf_dslash_5_plus.C,v 1.7 2004/10/14 22:03:48 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: dwf_dslash_5_plus.C,v $
//  $Revision: 1.7 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/qcdoc/dwf_dslash_5_plus.C,v $
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
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/smalloc.h>
#include<comms/scu.h>
#include<util/dirac_op.h>
CPS_START_NAMESPACE

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>


extern "C" { 
  void FtV1pV2Skip(double*,double,const double*,const double*,int);
  void FtV1pV2Skip_asm(double* out,
		       const double *scale,
		       const double* V1,
		       const double* V2,
		       int ntwo_spin);
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
  IFloat *f_temp = dwf_lib_arg->comm_buf;
//  IFloat *f_temp= (IFloat *) fmalloc(ls_stride*sizeof(IFloat));
//  IFloat *comm_buf = dwf_lib_arg->comm_buf;
  IFloat two_over_a5 = 2.0 * GJP.DwfA5Inv();
  IFloat neg_mass_two_over_a5 = -2.0 * mass * GJP.DwfA5Inv();

//  VRB.Func("DiracOpDwf","dwf_dslash_5_plus()");
  f_in  = (IFloat *) in;

  if (s_nodes>1){
    sys_cacheflush(0);
    (dwf_lib_arg->PlusArg[0])->Addr(f_in+(local_ls-1)*ls_stride);
    (dwf_lib_arg->MinusArg[0])->Addr(f_in);
  }

// [1 + gamma_5] term (if dag=1 [1 - gamma_5] term)
//
// out[s] = [1 + gamma_5] in[s-1]
//------------------------------------------------------------------

  SCUDirArgMulti *Plus = dwf_lib_arg->Plus;
  if (s_nodes >1)  Plus->StartTrans();
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  if(dag == 1){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  f_out = f_out + ls_stride; 

  struct timeval start, stop;
  //  gettimeofday(&start,NULL);

  //  FtV1pV2Skip(f_out,two_over_a5,f_in,f_out,(local_ls-1)*vol_4d_cb);
  FtV1pV2Skip_asm(f_out,&two_over_a5,f_in,f_out,(local_ls-1)*vol_4d_cb);

// [1 + gamma_5] for lower boundary term (if dag=1 [1 - gamma_5] term)
// If there's only one node along fifth direction, no communication
// is necessary; Otherwise data from adjacent node in minus direction
// will be needed.
// If the lower boundary is the s=0 term
// out[0] = - m_f * [1 + gamma_5] in[ls-1]
// else, out[s] = [1 + gamma_5] in[s-1]
//
//------------------------------------------------------------------

  if( s_nodes>1) f_in = f_temp;
  else{
    f_in  = (IFloat *) in;  
    f_in = f_in + (local_ls-1)*ls_stride; 
  }
  f_out = (IFloat *) out;
  
  if(dag == 1){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
    if (s_nodes >1) Plus->TransComplete();
 
    if(s_node_coor == 0) { 
      FtV1pV2Skip_asm(f_out, &neg_mass_two_over_a5, f_in, f_out, vol_4d_cb);
    }
    else {
      FtV1pV2Skip_asm(f_out, &two_over_a5, f_in, f_out, vol_4d_cb);
    }



// [1 - gamma_5] term (if dag=1 [1 + gamma_5] term)
// 
// out[s] = [1 - gamma_5] in[s+1]
//------------------------------------------------------------------
  SCUDirArgMulti *Minus = dwf_lib_arg->Minus;
  if (s_nodes >1) Minus->StartTrans();
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  f_in = f_in + ls_stride;

  FtV1pV2Skip_asm(f_out,&two_over_a5,f_in,f_out,(local_ls-1)*vol_4d_cb);

// [1 - gamma_5] for upper boundary term (if dag=1 [1 + gamma_5] term)
// If there's only one node along fifth direction, no communication
// is necessary; Otherwise data from adjacent node in minus direction
// will be needed.
// If the upper boundary is the s=ls term
// out[ls-1] = - m_f * [1 - gamma_5] in[0]
// else out[s] = [1 - gamma_5] in[s+1]
//
//------------------------------------------------------------------


  if (s_nodes >1) f_in  = f_temp;
  else   f_in  = (IFloat *) in;
  f_out = (IFloat *) out;

  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }

  f_out = f_out + (local_ls-1)*ls_stride;

  if (s_nodes >1) Minus->TransComplete();
  if(s_node_coor == s_nodes - 1) { 
    FtV1pV2Skip_asm(f_out, &neg_mass_two_over_a5, f_in, f_out, vol_4d_cb);
  } else {
    FtV1pV2Skip_asm(f_out, &two_over_a5, f_in, f_out, vol_4d_cb);
  }
  DiracOp::CGflops+=local_ls*vol_4d_cb*24*2;

  //  gettimeofday(&stop,NULL);

  //  double flops = local_ls*vol_4d_cb*24*2;
  //  double micros = stop.tv_usec - start.tv_usec;
  //  micros +=  (stop.tv_sec - start.tv_sec)*1.E6;

  //  printf("Dslash_5: %f mflops\n",flops/micros );
}



CPS_END_NAMESPACE
