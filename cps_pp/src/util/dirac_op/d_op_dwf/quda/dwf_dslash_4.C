#include<config.h>
#ifdef USE_SSE
#include "../sse/sse-dwf_dslash_4.C"
#else
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-12-05 19:55:38 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/quda/dwf_dslash_4.C,v 1.2 2012-12-05 19:55:38 chulwoo Exp $
//  $Id: dwf_dslash_4.C,v 1.2 2012-12-05 19:55:38 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: dwf_dslash_4.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/quda/dwf_dslash_4.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf_dslash_4.C
//
// dwf_dslash_4 is the derivative part of the space-time part of
// the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//------------------------------------------------------------------

#include<config.h>
#include<util/dwf.h>
#include<util/wilson.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
CPS_START_NAMESPACE

void wilson_dslash_vec(IFloat *chi_p_f,
                        IFloat *u_p_f,
                        IFloat *psi_p_f,
                        int cb,
                        int dag,
                        Wilson *wilson_p,
                        int vec_len,
                        unsigned long vec_offset);

void dwf_dslash_4(Vector *out, 
		  Matrix *gauge_field, 
		  Vector *in, 
		  int cb, 
		  int dag, 
		  Dwf *dwf_lib_arg)
{
  int i;
  int ls;
  IFloat *frm_in;
  IFloat *frm_out;
  IFloat *g_field;
  Wilson *wilson_p;
  int size_cb[2];
  int parity;

  //----------------------------------------------------------------
  // Initializations
  //----------------------------------------------------------------
  ls = dwf_lib_arg->ls;
  frm_in = (IFloat *) in;
  frm_out = (IFloat *) out;
  g_field = (IFloat *) gauge_field;
  wilson_p = dwf_lib_arg->wilson_p;
  size_cb[0] = 24*wilson_p->vol[0];
  size_cb[1] = 24*wilson_p->vol[1];
  
//#ifndef USE_TEST
#if 1
  //----------------------------------------------------------------
  // Apply 4-dimensional Dslash
  //----------------------------------------------------------------
  int vec_len=1;
  for(i=0; i<ls; i+= vec_len){

    // parity of 4-D checkerboard
    //------------------------------------------------------------
    parity = (i + cb) % 2;

    // Apply on 4-dim "parity" checkerboard part
    //------------------------------------------------------------
  if(vec_len==1)
    wilson_dslash(frm_out, g_field, frm_in, parity, dag, wilson_p);
  else{
#if TARGET == NOARCH
    ERR.NotImplemented("","dwf_dslash_4(..)","wilson_dslash_two() doesn't exists\n");
#else
    wilson_dslash_two(frm_out, frm_out+size_cb[parity], g_field, frm_in, frm_in+size_cb[parity], parity, 1-parity,dag, wilson_p);
#endif
    }
    frm_in = frm_in + vec_len*size_cb[parity];
    frm_out = frm_out + vec_len*size_cb[parity];
  }
#else
  //----------------------------------------------------------------
  // Apply vectorized 4-dimensional Dslash
  //----------------------------------------------------------------
  wilson_dslash_vec(frm_out, g_field, frm_in, cb, dag, wilson_p,ls/2,2*size_cb[parity]);
  frm_in = frm_in + size_cb[parity];
  frm_out = frm_out + size_cb[parity];
  wilson_dslash_vec(frm_out, g_field, frm_in, 1-cb, dag, wilson_p,ls/2,2*size_cb[parity]);
#endif


}

CPS_END_NAMESPACE
#endif
