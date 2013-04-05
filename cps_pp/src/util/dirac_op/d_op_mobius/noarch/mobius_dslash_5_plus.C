#ifdef USE_SSE
#include "mobius_dslash_5_plus-nonowait.C"
#else
#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_dslash_5_plus.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// mobius_dslash_5_plus.C
//
// mobius_dslash_5_plus is the derivative part of the 5th direction
// part of the fermion matrix. This routine accumulates the result
// on the out field 
// The in, out fields are defined on the checkerboard lattice.
// The action of this operator is the same for even/odd
// checkerboard fields because there is no gauge field along
// the 5th direction.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//
// Storage order for MOBIUS fermions
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

//  Also works for 4d preconditioned operator,
//  |                |
//  | |s-block odd|  |  For odd chckerboard
//  | |s-block odd| |
//  | |s-block odd|  |
//  |       .        |
//  |       .        |
//  |       .        |
//  |                |

//
//------------------------------------------------------------------


CPS_END_NAMESPACE
#include<util/mobius.h>
#include<util/dwf.h>
#include<util/gjp.h>
#include<util/dirac_op.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/smalloc.h>
#include<comms/scu.h>
CPS_START_NAMESPACE


void mobius_kappa_dslash_5_plus(Vector *out, 
				Vector *in, 
				Float mass,
				int dag, 
				Dwf *mobius_lib_arg,
				Float fact)
{
  int x;
  int s;

// Initializations
//------------------------------------------------------------------
  const int local_ls = GJP.SnodeSites(); 
  const int s_nodes = GJP.Snodes();
  const int s_node_coor = GJP.SnodeCoor();
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  IFloat *f_in;
  IFloat *f_out;
  IFloat *f_temp;
  IFloat *comm_buf = mobius_lib_arg->comm_buf;
  const IFloat neg_mass = -mass*fact;

// here we are acting on psi_in with P_L,R = 1/2 (1-+g5)  as
// opposed to dwf where the result is multiplied by 2

//
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
  for(s=1; s < local_ls; s++){

    for(x=0; x<vol_4d_cb; x++){

      fTimesV1PlusV2(f_out, fact, f_in, f_out, 12);

      f_in  =  f_in + 24;
      f_out = f_out + 24;
    }
  }


//
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
  
  for(x=0; x<vol_4d_cb; x++){
    
    f_temp = f_in;
    
    if (s_nodes != 1 ) {
      f_temp = comm_buf;
      getMinusData(f_temp, f_in, 12, 4);
    }
    
    if(s_node_coor == 0) { 
      fTimesV1PlusV2(f_out, neg_mass, f_temp, f_out, 12);
    }
    else {
      fTimesV1PlusV2(f_out, fact, f_temp, f_out, 12);
    }
    
    f_in  =  f_in + 24;
    f_out = f_out + 24;
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
  for(s=0; s < local_ls - 1; s++){

    for(x=0; x<vol_4d_cb; x++){

      fTimesV1PlusV2(f_out, fact, f_in, f_out, 12);

      f_in  =  f_in + 24;
      f_out = f_out + 24;
    }
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
  for(x=0; x<vol_4d_cb; x++){

    f_temp = f_in;
   
    if (s_nodes != 1 ) {
      f_temp = comm_buf;
      getPlusData(f_temp, f_in, 12, 4);
    }

    if(s_node_coor == s_nodes - 1) { 
      fTimesV1PlusV2(f_out, neg_mass, f_temp, f_out, 12);
    }
    else {
      fTimesV1PlusV2(f_out, fact, f_temp, f_out, 12);
    }
    
    f_in  =  f_in + 24;
    f_out = f_out + 24;
  }
  DiracOp::CGflops+=2*2*vol_4d_cb*local_ls*12;
  

}

void mobius_dslash_5_plus(Vector *out, 
			  Vector *in, 
			  Float mass,
			  int dag, 
			  Dwf *mobius_lib_arg)
{
  mobius_kappa_dslash_5_plus(out, in, mass, dag, mobius_lib_arg, GJP.DwfA5Inv());
}



CPS_END_NAMESPACE

#endif
