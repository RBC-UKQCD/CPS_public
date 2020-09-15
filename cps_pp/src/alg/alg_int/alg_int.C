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

std::vector<int> AlgInt::phi_veloc_all;
std::vector<int> AlgInt::md_veloc_all;
int AlgInt::veloc_id=0;
int AlgInt::traj_num=0;

AlgInt::AlgInt()
{
  cname = "AlgInt()";
//  traj = -1;
}

AlgInt::~AlgInt()
{

}

void AlgInt::copyConjLattice(){
  //do nothing for all derived classes bar AlgActionGauge 
}

// Calculate preliminary force for force gradient. This is just a trap
// for all derived classes that don't implement force gradient
// evolution.
void AlgInt::prepare_fg(Matrix * force, Float dt_ratio)
{
  const char fname[] = "prepare_fg(M*,F)";
  ERR.General(cname, fname, "Force gradient evolution not defined.\n");
}

/*
AlgInt operator+(AlgInt &A, AlgInt &B) {
  Lattice &lat = AlgLattice();
  AlgIntSum sum(&A, &B, lat, c_arg);
  return sum;
}
*/

CPS_END_NAMESPACE

