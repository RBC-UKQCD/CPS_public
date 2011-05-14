#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_int_omelyan.C
//
// AlgIntOmelyan is derived from AlgIntAB, it is an implementation of
// the omelyan integrator using abstract operators.
// 
// To construct a QPQPQ integrator, the update to the coordinate
// must be the first argument, the momentum the second.
//
// To construct a PQPQP integrator, the update to the momentum
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


AlgIntOmelyan::AlgIntOmelyan(AlgInt &A, AlgInt &B, 
			     IntABArg &arg_ab) 
  : AlgIntAB(A,B,arg_ab)
{

  int_type = INT_OMELYAN;
  A_calls = 3;
  B_calls = 2;
  lambda = ab_arg->lambda;
  
}

AlgIntOmelyan::~AlgIntOmelyan() {

}

void AlgIntOmelyan::evolve(Float dt, int steps) 
{
  char * fname = "evolve(Float, int)";

  step_cnt = 0;
  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(step_cnt);

  A->evolve(lambda*dt/(Float)A_steps, A_steps);

  for (int i=0; i<steps; i++) {
    if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    
    B->evolve(dt/(2.0*(Float)B_steps), B_steps);
    A->evolve((1-2*lambda)*dt/(Float)A_steps, A_steps);
    B->evolve(dt/(2.0*(Float)B_steps), B_steps);
    if (i < steps-1) A->evolve(2.0*lambda*dt/(Float)A_steps, A_steps);
    else A->evolve(lambda*dt/(Float)A_steps, A_steps);
  }

  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    
}

CPS_END_NAMESPACE
