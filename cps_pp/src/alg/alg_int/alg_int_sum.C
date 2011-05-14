#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_int_sum.C
//
// AlgIntSum is derived from AlgInt, it is produces the composite
// integrator of A and B
// 
// If [A,B] = [B,A] then the order of the parameters to the
// constructor obviously does not matter (i.e., both operator
// corresponding to a momentum update.  If they do not commute, then
// the care must be taken when the constructor.
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
#include<util/checksum.h>
CPS_START_NAMESPACE


AlgIntSum::AlgIntSum(AlgInt &A, AlgInt &B, IntABArg &ab_arg)
  : AlgIntAB(A,B,ab_arg)
{
  A_calls = 1;
  B_calls = 1;
}

AlgIntSum::~AlgIntSum() {

}

void AlgIntSum::prepare_fg(Matrix * force, Float dt_ratio)
{
  A -> prepare_fg(force, dt_ratio);
  B -> prepare_fg(force, dt_ratio);
}

void AlgIntSum::evolve(Float dt, int steps) 
{
  char * fname = "evolve(Float, int)";
  step_cnt = 0;
  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(step_cnt);
  for (int i=0; i<steps; i++) {
    if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    A -> evolve(dt, A_steps);
    B -> evolve(dt, B_steps);
  }
  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);

}

CPS_END_NAMESPACE
