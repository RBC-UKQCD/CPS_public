#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_mdag.C,v $
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
// ( M_5 - kappa_b^2 M4eo M_5^-1 M4oe )^ dag
void  mobius_mdag(Vector *out, 
		  Matrix *gauge_field, 
		  Vector *in, 
		  Float mass, 
		  Dwf *mobius_lib_arg)
{
  

  //------------------------------------------------------------------
  // Initializations
  //------------------------------------------------------------------
  const int ls = mobius_lib_arg->ls;
  const size_t f_size = 24 * mobius_lib_arg->vol_4d * ls / 2;

  zmobius_mdag(out,gauge_field,in,mass,mobius_lib_arg);


  // Flops count in this function is two AXPY = 4 flops per vector elements
  DiracOp::CGflops +=  3*f_size; 

}






CPS_END_NAMESPACE
