#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_int_leap.C
//
// AlgIntLeap is derived from AlgInt, it is an implementation
// of the leapfrog integrator using abstract operators, i.e.
//
// To construct a QPQ integrator, the update to the coordinate
// must be the first argument, the momentum the second.
//
// To construct a PQP integrator, the update to the momentum
// must be the second argument, the coordinate the second.
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


AlgIntLeap::AlgIntLeap(AlgInt &a, AlgInt &b, IntABArg &ab_arg) 
  : AlgIntAB(a,b,ab_arg)
{
  int_type = INT_LEAP;
  A_calls = 2;
  B_calls = 1;
}

AlgIntLeap::~AlgIntLeap() {

}

void AlgIntLeap::evolve(Float dt, int steps) 
{
  char * fname = "evolve(Float, int)";

  step_cnt = 0;
  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(step_cnt);

  A->evolve(0.5*dt/(Float)A_steps, A_steps);

  for (int i=0; i<steps; i++) {
    if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);

    B->evolve(dt/(Float)B_steps, B_steps);    
    if (i < steps-1) A->evolve(dt/(Float)A_steps, A_steps);
    else A->evolve(0.5*dt/(Float)A_steps, A_steps);
  }

  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);

}

CPS_END_NAMESPACE
