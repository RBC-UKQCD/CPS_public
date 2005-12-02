#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_int.C
//
// AlgInt is abstract base class from which all integrators 
// are derived.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<math.h>
#include<alg/alg_hmd.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_int.h>
CPS_START_NAMESPACE

AlgInt::AlgInt()
{
  cname = "AlgInt()";
  traj = -1;
}

AlgInt::~AlgInt()
{

}

/*
AlgInt operator+(AlgInt &A, AlgInt &B) {
  Lattice &lat = AlgLattice();
  AlgIntSum sum(&A, &B, lat, c_arg);
  return sum;
}
*/

CPS_END_NAMESPACE

