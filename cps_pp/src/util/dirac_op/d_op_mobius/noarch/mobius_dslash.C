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
#include<util/mobius.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE


void mobius_dslash(Vector *out, 
		   Matrix *gauge_field, 
		   Vector *in, 
		   Float mass,
		   int cb, 
		   int dag, 
		   Dwf *mobius_lib_arg)
{
//------------------------------------------------------------------
// Apply 4-dimensional Dslash
//------------------------------------------------------------------
  mobius_dslash_4(out, gauge_field, in, cb, dag, mobius_lib_arg, mass);


//------------------------------------------------------------------
// Apply 5th-direction Dslash
//------------------------------------------------------------------
  Float kappa_c_inv = 1.0/mobius_lib_arg->mobius_kappa_c;
  mobius_kappa_dslash_5_plus(out, in, mass, dag, mobius_lib_arg, kappa_c_inv);

}

CPS_END_NAMESPACE
