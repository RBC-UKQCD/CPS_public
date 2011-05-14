#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_int_campostrini.C
//
// AlgIntCampostrini is derived from AlgIntAB, it is an implementation of
// the 4th order Campostrini integrator using abstract operators.
// 
// To construct a QPQPQPQ integrator, the update to the coordinate
// must be the first argument, the momentum the second.
//
// To construct a PQPQPQP integrator, the update to the momentum
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


AlgIntCampostrini::AlgIntCampostrini(AlgInt &A, AlgInt &B, 
			     IntABArg &arg_ab) 
  : AlgIntAB(A,B,arg_ab)
{

  int_type = INT_CAMPOSTRINI;
  A_calls = 4;
  B_calls = 3;
  sigma = pow(2.0, 1.0/3.0);
}

AlgIntCampostrini::~AlgIntCampostrini() {

}

void AlgIntCampostrini::evolve(Float dt, int steps) 
{
  char * fname = "evolve(Float, int)";
  epsilon = dt/(2.0 - sigma);

  step_cnt = 0;
  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(step_cnt);

  A->evolve(epsilon/(2.0*(Float)A_steps), A_steps);

  for (int i=0; i<steps; i++) {
    if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    
    B->evolve(epsilon/(Float)B_steps, B_steps);
    A->evolve((1-sigma)*epsilon/(2.0*(Float)A_steps), A_steps);
    B->evolve(-epsilon*sigma/(Float)B_steps, B_steps);
    A->evolve((1-sigma)*epsilon/(2.0*(Float)A_steps), A_steps);
    B->evolve(epsilon/(Float)B_steps, B_steps);
    if (i < steps-1) A->evolve(epsilon/(Float)A_steps, A_steps);
    else A->evolve(epsilon/(2.0*(Float)A_steps), A_steps);
  }

  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    
}

CPS_END_NAMESPACE
