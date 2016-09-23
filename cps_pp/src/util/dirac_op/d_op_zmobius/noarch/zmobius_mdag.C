#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_m.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// mobius_m.C
//
// mobius_m is the fermion matrix.  
// The in, out fields are defined on the checkerboard lattice
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/zmobius.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/dirac_op.h>
#include<util/time_cps.h>

#include "blas-subs.h"

CPS_START_NAMESPACE


#include "zmobius_mdag-orig.h"
#include "zmobius_mdag-sym1.h"
#include "zmobius_mdag-sym1-MIT.h"
#include "zmobius_mdag-sym2.h"
#include "zmobius_mdag-sym2-MIT.h"
#include "zmobius_mdag-sym3.h"


void  zmobius_mdag(Vector *out, 
		   Matrix *gauge_field, 
		   Vector *in, 
		   Float mass, 
		   Zmobus *mobius_lib_arg)
{
  switch( mobius_lib_arg-> pc_type ){
  case   ZMOB_PC_ORIG:
    zmobius_mdag_orig(out, gauge_field, in, mass, mobius_lib_arg);
    break;
  case   ZMOB_PC_SYM1:
    zmobius_mdag_sym1(out, gauge_field, in, mass, mobius_lib_arg);
    break;
  case   ZMOB_PC_SYM2:
    zmobius_mdag_sym2(out, gauge_field, in, mass, mobius_lib_arg);
    break;
  case   ZMOB_PC_SYM1_MIT:
    zmobius_mdag_sym1_MIT(out, gauge_field, in, mass, mobius_lib_arg);
    break;
  case   ZMOB_PC_SYM2_MIT:
    zmobius_mdag_sym2_MIT(out, gauge_field, in, mass, mobius_lib_arg);
    break;
  case   ZMOB_PC_SYM3:
    zmobius_mdag_sym3(out, gauge_field, in, mass, mobius_lib_arg);
    break;
  default:
    ERR.NotImplemented("","zmobius_mdag(...)");
  }    
}



CPS_END_NAMESPACE
