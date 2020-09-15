#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
// mobius_m.C
//
// mobius_m is the fermion matrix.  
// The in, out fields are defined on the checkerboard lattice
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/mobius.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/dirac_op.h>
#include<util/time_cps.h>

#include "blas-subs.h"

CPS_START_NAMESPACE


//4d precond. mobius Dirac op:
// M_5 - kappa_b^2 M4eo M_5^-1 M4oe
void  mobius_m(Vector *out, 
	       Matrix *gauge_field, 
	       Vector *in, 
	       Float mass, 
	       Dwf *mobius_lib_arg)
{

  const char *fname="mobius_m()";
//currently Zmobus = Dwf
  VRB.Debug("",fname,"pc=%d \n",mobius_lib_arg ->pc_type);
//  if (mobius_lib_arg ->pc_type != ZMOB_PC_ORIG)
//  ERR.General("",fname,"Only ZMOB_PC_ORIG tested\n");
  zmobius_m(out,gauge_field,in,mass,mobius_lib_arg);

}

CPS_END_NAMESPACE
