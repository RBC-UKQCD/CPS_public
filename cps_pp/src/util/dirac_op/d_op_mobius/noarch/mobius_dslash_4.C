#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:30 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_dslash_4.C,v 1.2 2013-04-05 17:46:30 chulwoo Exp $
//  $Id: mobius_dslash_4.C,v 1.2 2013-04-05 17:46:30 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: mobius_dslash_4.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_dslash_4.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// mobius_dslash_4.C
//
// mobius_dslash_4 is the derivative part of the space-time part of
// the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<config.h>
#include<util/mobius.h>
#include<util/dwf.h>
#include<util/wilson.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<comms/scu.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif


void mobius_dslash_4(Vector *out, 
		     Matrix *gauge_field, 
		     Vector *in, 
		     int cb, 
		     int dag, 
		     Dwf *mobius_lib_arg,
		     Float mass)
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
  ls = mobius_lib_arg->ls;
  const int f_size = 24 * mobius_lib_arg->vol_4d / 2;
  Float b_coeff = GJP.Mobius_b();
  Float c_coeff = GJP.Mobius_c();

  frm_in = (IFloat *) in;
  frm_out = (IFloat *) out;
  g_field = (IFloat *) gauge_field;
  wilson_p = mobius_lib_arg->wilson_p;
  size_cb[0] = 24*wilson_p->vol[0];
  size_cb[1] = 24*wilson_p->vol[1];
  
  IFloat* frm_;
  Vector  *frm_tmp3 = (Vector *) mobius_lib_arg->frm_tmp3;
  frm_ = (IFloat*)frm_tmp3;

  //----------------------------------------------------------------
  // Apply 4-dimensional Dslash
  //----------------------------------------------------------------

  // frm_ = b * Psi(s)
  vecEqualsVecTimesEquFloat(frm_, frm_in, b_coeff, f_size*ls);
  // frm_ += c * P_L * Psi(s+1) + c * P_R * Psi(s-1)
  mobius_kappa_dslash_5_plus((Vector*)frm_, in, mass, dag, mobius_lib_arg, c_coeff);

  // out = D_W * frm_
  for(i=0; i<ls; i++){

    // parity of 4-D checkerboard
    //------------------------------------------------------------
    parity = cb;//4d odd-even preconditioning

    // Apply on 4-dim "parity" checkerboard part
    //------------------------------------------------------------
    wilson_dslash(frm_out, g_field, frm_, parity, dag, wilson_p);
    
    frm_ = frm_ + size_cb[parity];
    frm_out = frm_out + size_cb[parity];
  }
  
}



CPS_END_NAMESPACE
