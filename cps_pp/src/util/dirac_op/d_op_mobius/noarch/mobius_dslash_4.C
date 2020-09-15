#include<config.h>
CPS_START_NAMESPACE
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
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE

void mobius_dslash_4(Vector *out, 
		     Matrix *gauge_field, 
		     Vector *in, 
		     int cb, 
		     int dag, 
		     Dwf *mobius_lib_arg,
		     Float mass)
{

    zmobius_dslash_4(out,gauge_field,in,cb,dag,mobius_lib_arg,mass);
}

void mobius_Booee(Vector *out, 
		     Vector *in, 
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
  const size_t f_size = 24 * mobius_lib_arg->vol_4d / 2;
  Float b_coeff = GJP.Mobius_b();
  Float c_coeff = GJP.Mobius_c();
  VRB.Debug("","mobius_Booee","b_coeff=%g c_coeff=%g\n",b_coeff,c_coeff);
  if(dag)
  ERR.General("","mobius_Booee", "Only implemented for dag=0\n");

  frm_in = (IFloat *) in;
  frm_out = (IFloat *) out;
  wilson_p = mobius_lib_arg->wilson_p;
  size_cb[0] = 24*wilson_p->vol[0];
  size_cb[1] = 24*wilson_p->vol[1];
  
  IFloat* frm_;
//  Vector  *frm_tmp3 = (Vector *) mobius_lib_arg->frm_tmp3;
  frm_ = (IFloat*)frm_out;

  //----------------------------------------------------------------
  // Apply 4-dimensional Dslash
  //----------------------------------------------------------------

  // frm_ = b * Psi(s)
  vecEqualsVecTimesEquFloat(frm_, frm_in, b_coeff, f_size*ls);
  // frm_ += c * P_L * Psi(s+1) + c * P_R * Psi(s-1)
  mobius_kappa_dslash_5_plus((Vector*)frm_, in, mass, dag, mobius_lib_arg, c_coeff);

}

CPS_END_NAMESPACE
