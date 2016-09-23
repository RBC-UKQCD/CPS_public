#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/bgl/dwf_dslash_5_plus_slice.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf_dslash_5_plus_slice.C
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

extern "C" {
  void FtV1pV2Skip(double*,double,const double*,const double*,int);
  void FtV1pV2Skip_asm(double* out,
               const double *scale,
               const double* V1,
               const double* V2,
               int ntwo_spin);
}

static  int local_ls;
static  int s_nodes;
static  int s_node_coor;
static  int vol_4d_cb;
static  int ls_stride;
static  Vector *in_p = NULL;
static  Float * rbuf_up = NULL;
static  Float * rbuf_down = NULL;
static  QMP_msgmem_t msgmem_up[2];
static  QMP_msgmem_t msgmem_down[2];
static  QMP_msghandle_t msghandle_up[2];
static  QMP_msghandle_t msghandle_down[2];

#undef USE_GETPLUS
void dwf_dslash_5_plus_start(Vector *out, 
		       Vector *in, 
		       Float mass,
		       int dag, 
		       Dwf *dwf_lib_arg)
{
  int x;
  int s;

//    getPlusData(comm_buf, f_in, 24*vol_4d_cb, 4);
//  f_in  = (IFloat *) in;  
//  f_in = f_in + (local_ls-1)*ls_stride; 
//    getMinusData(comm_buf, f_in, 24*vol_4d_cb, 4);
// Initializations
//------------------------------------------------------------------
  local_ls = GJP.SnodeSites();
  s_nodes = GJP.Snodes();
  s_node_coor = GJP.SnodeCoor();
  vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  ls_stride = 24 * vol_4d_cb;
  if(s_nodes<2) return;
  if(!rbuf_up){
    rbuf_up = (Float *)malloc(24*vol_4d_cb*sizeof(IFloat));
    msgmem_up[1] = QMP_declare_msgmem(rbuf_up,24*vol_4d_cb*sizeof(IFloat));
    msghandle_up[1] = QMP_declare_receive_relative(msgmem_up[1],4,+1,0);
//    fprintf(stderr,"rbuf_up=%p\n",rbuf_up);
//    fprintf(stderr,"msgmem_up[1]=%p\n",msgmem_up[1]);
//    fprintf(stderr,"msghandle_up[1]=%p\n",msghandle_up[1]);
  }
  if(!rbuf_down){
    rbuf_down = (Float *)malloc(24*vol_4d_cb*sizeof(IFloat));
    msgmem_down[1] = QMP_declare_msgmem(rbuf_down,24*vol_4d_cb*sizeof(IFloat));
    msghandle_down[1] = QMP_declare_receive_relative(msgmem_down[1],4,-1,0);
  }

  if (in_p != in){
    if (in_p!=NULL){
      QMP_free_msghandle(msghandle_up[0]); 
      QMP_free_msgmem(msgmem_up[0]);   
      QMP_free_msghandle(msghandle_down[0]); 
      QMP_free_msgmem(msgmem_down[0]);   
    }
    in_p = in;
    IFloat *f_in = (IFloat *)in;
    msgmem_up[0] = QMP_declare_msgmem(f_in,24*vol_4d_cb*sizeof(IFloat));
    msghandle_up[0] = QMP_declare_send_relative(msgmem_up[0],4,-1,0);
    f_in = f_in + (local_ls-1)*ls_stride;
    msgmem_down[0] = QMP_declare_msgmem(f_in,24*vol_4d_cb*sizeof(IFloat));
    msghandle_down[0] = QMP_declare_send_relative(msgmem_down[0],4,1,0);
  }
#ifndef USE_GETPLUS
  QMP_start(msghandle_up[0]);
  QMP_start(msghandle_up[1]);
  QMP_start(msghandle_down[0]);
  QMP_start(msghandle_down[1]);
#endif
}


void dwf_dslash_5_plus_slice(Vector *out, 
		       Vector *in, 
		       Float mass,
		       int dag, 
		       Dwf *dwf_lib_arg, 
		       int s_slice)
{
  int x;
  int s;

// Initializations
//------------------------------------------------------------------
#if 0
  int local_ls = GJP.SnodeSites(); 
  int s_nodes = GJP.Snodes();
  int s_node_coor = GJP.SnodeCoor();
  int vol_4d_cb = dwf_lib_arg->vol_4d / 2;
  int ls_stride = 24 * vol_4d_cb;
#endif

  IFloat *f_in;
  IFloat *f_out;
  IFloat *f_temp;
  IFloat *comm_buf = dwf_lib_arg->comm_buf;
  IFloat two_over_a5 = 2.0 * GJP.DwfA5Inv();
  IFloat neg_mass_two_over_a5 = -2.0 * mass * GJP.DwfA5Inv();

// [1 + gamma_5] term (if dag=1 [1 - gamma_5] term)
//
// out[s] = [1 + gamma_5] in[s-1]
//------------------------------------------------------------------

  if (s_slice<0 || s_slice >=local_ls)
  ERR.General("","dwf_dslash_5_plus_slice","s_slice=%d local_ls=%d!\n",s_slice,local_ls);
  
  if(s_slice>0 ){
    f_in  = (IFloat *) in;
    f_out = (IFloat *) out;
    f_in += (s_slice-1)*ls_stride;
    f_out += (s_slice)*ls_stride;
    if(dag == 1){
      f_in  =  f_in + 12;
      f_out = f_out + 12;
    }
    FtV1pV2Skip_asm(f_out,&two_over_a5,f_in,f_out,vol_4d_cb);
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

if (s_slice ==  0 ){
  f_in  = (IFloat *) in;  
  f_in = f_in + (local_ls-1)*ls_stride; 
  f_out = (IFloat *) out;
  
  if(dag == 1){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  
  f_temp = f_in;
  if (s_nodes > 1 ) {
#ifdef USE_GETPLUS
    getMinusData(comm_buf, f_in, 24*vol_4d_cb, 4);
    f_temp = comm_buf;
#else
      QMP_status_t send_status = QMP_wait(msghandle_down[0]);
      if (send_status != QMP_SUCCESS) 
      QMP_error("Send failed in dwf_dslash_5_plus_slice: %s\n", QMP_error_string(send_status));
      QMP_status_t recv_status = QMP_wait(msghandle_down[1]);
    if (recv_status != QMP_SUCCESS)
      QMP_error("Receive failed in dwf_dslash_5_plus_slice: %s\n", QMP_error_string(recv_status));
     f_temp = rbuf_down; 
     if(dag == 1) f_temp = f_temp + 12;
#endif
  }
    if(s_node_coor == 0) { 
      FtV1pV2Skip_asm(f_out,&neg_mass_two_over_a5,f_temp,f_out,vol_4d_cb);
    } else {
      FtV1pV2Skip_asm(f_out,&two_over_a5,f_temp,f_out,vol_4d_cb);
    }
}


// [1 - gamma_5] term (if dag=1 [1 + gamma_5] term)
// 
// out[s] = [1 - gamma_5] in[s+1]
//------------------------------------------------------------------
if(s_slice > 0 ){
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;
  f_in += (s_slice)*ls_stride;
  f_out += (s_slice-1)*ls_stride;
  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }
  FtV1pV2Skip_asm(f_out,&two_over_a5,f_in,f_out,vol_4d_cb);
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

if(s_slice == (local_ls-1) ){
  f_in  = (IFloat *) in;
  f_out = (IFloat *) out;

  if(dag == 0){
    f_in  =  f_in + 12;
    f_out = f_out + 12;
  }

  f_out = f_out + (local_ls-1)*ls_stride;
    f_temp = f_in;
  if (s_nodes > 1 ) {
#ifdef USE_GETPLUS
    getPlusData(comm_buf, f_in, 24*vol_4d_cb, 4);
    f_temp = comm_buf;
#else
      QMP_status_t send_status = QMP_wait(msghandle_up[0]);
      if (send_status != QMP_SUCCESS) 
      QMP_error("Send failed in dwf_dslash_5_plus_slice: %s\n", QMP_error_string(send_status));
      QMP_status_t recv_status = QMP_wait(msghandle_up[1]);
    if (recv_status != QMP_SUCCESS)
      QMP_error("Receive failed in dwf_dslash_5_plus_slice: %s\n", QMP_error_string(recv_status));
     f_temp = rbuf_up; 
     if(dag == 0) f_temp = f_temp + 12;
#endif
  }
    if(s_node_coor == s_nodes - 1) { 
      FtV1pV2Skip_asm(f_out,&neg_mass_two_over_a5,f_temp,f_out,vol_4d_cb);
    } else {
      FtV1pV2Skip_asm(f_out,&two_over_a5,f_temp,f_out,vol_4d_cb);
    }
}
//  DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;


}



CPS_END_NAMESPACE
