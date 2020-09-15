#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_dslash.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// mobius_dslash.C
//
// mobius_dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the full lattice
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/zmobius.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE


void zmobius_dslash(Vector *out, 
		   Matrix *gauge_field, 
		   Vector *in, 
		   Float mass,
		   int cb, 
		   int dag, 
		   Zmobus *mobius_lib_arg)
{
//------------------------------------------------------------------
// Apply 4-dimensional Dslash
//------------------------------------------------------------------
  zmobius_dslash_4(out, gauge_field, in, cb, dag, mobius_lib_arg, mass);


//------------------------------------------------------------------
// Apply 5th-direction Dslash
//------------------------------------------------------------------
  Complex* kappa_ratio = mobius_lib_arg->zmobius_kappa_ratio.data();
  zmobius_kappa_dslash_5_plus_cmplx(out, in, mass, dag, mobius_lib_arg, kappa_ratio);

}

CPS_END_NAMESPACE
