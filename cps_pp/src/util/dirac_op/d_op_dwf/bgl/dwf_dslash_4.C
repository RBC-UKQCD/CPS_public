#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-02-08 18:35:07 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/bgl/dwf_dslash_4.C,v 1.5 2008-02-08 18:35:07 chulwoo Exp $
//  $Id: dwf_dslash_4.C,v 1.5 2008-02-08 18:35:07 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: dwf_dslash_4.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/bgl/dwf_dslash_4.C,v $
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

CPS_END_NAMESPACE
#include<config.h>
#include<util/dwf.h>
#include<util/wilson.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif
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
  
  //----------------------------------------------------------------
  // Apply 4-dimensional Dslash
  //----------------------------------------------------------------
#if TARGET == BGL
  int vec_len=1;
#else
  int vec_len=2;
#endif
  for(i=0; i<ls; i+= vec_len){

    // parity of 4-D checkerboard
    //------------------------------------------------------------
    parity = (i + cb) % 2;

    // Apply on 4-dim "parity" checkerboard part
    //------------------------------------------------------------
  if(vec_len==1)
    wilson_dslash(frm_out, g_field, frm_in, parity, dag, wilson_p);
  else
    wilson_dslash_two(frm_out, frm_out+size_cb[parity], g_field, frm_in, frm_in+size_cb[parity], parity, 1-parity,dag, wilson_p);
    frm_in = frm_in + vec_len*size_cb[parity];
    frm_out = frm_out + vec_len*size_cb[parity];
  }


}



CPS_END_NAMESPACE
