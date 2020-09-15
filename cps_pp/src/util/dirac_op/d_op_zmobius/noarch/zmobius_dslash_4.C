#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// zmobius_dslash_4.C
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
#include<util/zmobius.h>
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



void zmobius_dslash_4_dag0(Vector *out, 
		     Matrix *gauge_field, 
		     Vector *in, 
		     int cb, 
		     Zmobus *mobius_lib_arg,
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
  const char *fname="zmobius_dslash_4_dag0";


  //----------------------------------------------------------------
  // Initializations
  //----------------------------------------------------------------
  const int dag=0;
  ls = mobius_lib_arg->ls;
  const size_t f_size = 24 * mobius_lib_arg->vol_4d / 2;
#if 1
//  Float b_coeff = GJP.Mobius_b();
//  Float c_coeff = GJP.Mobius_c();
  Complex* b_coeff = mobius_lib_arg->zmobius_b.data();
  Complex* c_coeff = mobius_lib_arg->zmobius_c.data();
#else
  Complex* b_coeff = GJP.ZMobius_b();
  Complex* c_coeff = GJP.ZMobius_c();
#endif

  frm_in = (IFloat *) in;
  frm_out = (IFloat *) out;
  g_field = (IFloat *) gauge_field;
  wilson_p = mobius_lib_arg->wilson_p;
  size_cb[0] = 24*wilson_p->vol[0];
  size_cb[1] = 24*wilson_p->vol[1];
  
  IFloat* frm_;
  Vector  *frm_tmp3 = (Vector *) mobius_lib_arg->frm_tmp3;
//  VRB.Debug("",fname,"frm_tmp3=%p\n", mobius_lib_arg->frm_tmp3);

  frm_ = (IFloat*)frm_tmp3;

  //----------------------------------------------------------------
  // Apply 4-dimensional Dslash
  //----------------------------------------------------------------
#if 0
    // frm_ = b(s) * Psi(s)
    vecEqualsVecTimesEquFloat((Complex*)frm_, (Complex*)frm_in, b_coeff, f_size*ls);
    // frm_ += c(s) * P_L * Psi(s+1) + c(s) * P_R * Psi(s-1)
    zmobius_kappa_dslash_5_plus_cmplx((Vector*)frm_, in, mass, dag, mobius_lib_arg, c_coeff);
#else
    for(i=0; i<ls; i++){
      // frm_ = b(s) * Psi(s)
      // b or b^dagger

      Complex b;
      if(dag)  b = conj(b_coeff[i]);
      else   b = b_coeff[i]; 
      IFloat* _frm_ = frm_ + i* f_size;
      IFloat* _frm_in = frm_in + i* f_size;
//      printf("vecEqualsVecTimesEquComplex(%p %p (%g %g) %d)\n",
//	(Complex*)_frm_, (Complex*)_frm_in, b.real(), b.imag(), f_size);
      vecEqualsVecTimesEquComplex((Complex*)_frm_, (Complex*)_frm_in, b, f_size);
    }
    
    // frm_ += c(s) * P_L * Psi(s+1) + c(s) * P_R * Psi(s-1)
    zmobius_kappa_dslash_5_plus_cmplx((Vector*)frm_, in, mass, dag, mobius_lib_arg, c_coeff);
#endif

  // out = D_W * frm_
  for(i=0; i<ls; i++){
    // parity of 4-D checkerboard
    //------------------------------------------------------------
    parity = cb;//4d odd-even preconditioning

    // Apply on 4-dim "parity" checkerboard part
    //------------------------------------------------------------
//    VRB.Result("",fname,"wilson_dslash(%p %p %p %d %d %p)\n",frm_out, g_field, frm_, parity, dag, wilson_p);
    wilson_dslash(frm_out, g_field, frm_, parity, dag, wilson_p);
    
    frm_ = frm_ + size_cb[parity];
    frm_out = frm_out + size_cb[parity];
  }
  
}

void zmobius_dslash_4_dag1(Vector *out, 
		     Matrix *gauge_field, 
		     Vector *in, 
		     int cb, 
		     Zmobus *mobius_lib_arg,
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
  const int dag=1;
  ls = mobius_lib_arg->ls;
  const size_t f_size = 24 * mobius_lib_arg->vol_4d / 2;
#if 1
//  Float b_coeff = GJP.Mobius_b();
//  Float c_coeff = GJP.Mobius_c();
  Complex* b_coeff = mobius_lib_arg->zmobius_b.data();
  Complex* c_coeff = mobius_lib_arg->zmobius_c.data();
#else
  Complex* b_coeff = GJP.ZMobius_b();
  Complex* c_coeff = GJP.ZMobius_c();
#endif

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


  // frm_ = D_W * Psi(s)
  for(i=0; i<ls; i++){
    IFloat* _frm_ = frm_ + i* f_size;
    IFloat* _frm_in = frm_in + i* f_size;

    // parity of 4-D checkerboard
    //------------------------------------------------------------
    parity = cb;//4d odd-even preconditioning

    // Apply on 4-dim "parity" checkerboard part
    //------------------------------------------------------------
    wilson_dslash(_frm_, g_field, _frm_in, parity, dag, wilson_p);
  }
  
  // out = b(s) * Psi(s)
  // b or b^dagger
  for(i=0; i<ls; i++){
      Complex b;
      if(dag)  b = conj(b_coeff[i]);
      else   b = b_coeff[i]; 

      IFloat* _frm_ = frm_ + i* f_size;
      IFloat* _frm_out = frm_out + i* f_size;

      vecEqualsVecTimesEquComplex((Complex*)_frm_out, (Complex*)_frm_, b, f_size);
    }
    // out += c(s) * P_L * Psi(s+1) + c(s) * P_R * Psi(s-1)
  zmobius_kappa_dslash_5_plus_cmplx((Vector*)frm_out, (Vector*)frm_, mass, dag, mobius_lib_arg, c_coeff);


}



void zmobius_dslash_4(Vector *out, 
		     Matrix *gauge_field, 
		     Vector *in, 
		     int cb, 
		     int dag, 
		     Zmobus *mobius_lib_arg,
		      Float mass)
{
  if(!dag)
    zmobius_dslash_4_dag0(out, gauge_field, in, cb, mobius_lib_arg, mass);
  else
    zmobius_dslash_4_dag1(out, gauge_field, in, cb, mobius_lib_arg, mass);
}


CPS_END_NAMESPACE
