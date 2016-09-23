#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/02/08 18:35:07 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/qcdoc/dwf_dslash_4.C,v 1.5 2008/02/08 18:35:07 chulwoo Exp $
//  $Id: dwf_dslash_4.C,v 1.5 2008/02/08 18:35:07 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: dwf_dslash_4.C,v $
//  $Revision: 1.5 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/qcdoc/dwf_dslash_4.C,v $
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
#include<util/wfm.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif

#define MAX_LS 64

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
  /*
   * PAB: I don't want to modify the dwf_lib_arg if I can avoid.
   * Also a malloc/free per call is bad, so no choice but to put in a 
   * static limit;
   */
  Float *dwf_psis[MAX_LS];
  Float *dwf_chis[MAX_LS];
  int    dwf_cbs [MAX_LS];

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
  

  // I'm a bad programmer
  if ( size_cb[0] != size_cb[1] ) { 
    printf("dwf_dslash_4:size_cb[0](%d) != size_cb[1]\n(%d)\n",size_cb[0],size_cb[1]);
    exit(-1);
  }
  if ( ls > MAX_LS ) { 
    printf("dwf_dslash_4:ls(%d) is larger than MAX_LS(%d)\n",ls,MAX_LS);
    exit(-1);
  }

  for(i=0; i<ls; i++){
    dwf_cbs[i] = (i+cb)%2;
    dwf_chis[i]= frm_in  + i*size_cb[0];
    dwf_psis[i]= frm_out + i*size_cb[0];
  }

  //----------------------------------------------------------------
  // Apply optimised 4-dimensional Dslash
  //----------------------------------------------------------------

  wfm_dslash_vec(ls,dwf_psis,g_field,dwf_chis,dwf_cbs,dag);

}



CPS_END_NAMESPACE
