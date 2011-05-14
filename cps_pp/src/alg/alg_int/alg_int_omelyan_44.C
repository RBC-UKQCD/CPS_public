#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_int_omelyan_44.C
//
// AlgIntOmelyan44 is derived from AlgIntAB, it is an implementation
// of the 4th order Omelyan integrator using abstract operators.  This
// integrator uses 4 force evaluations per step, see AlgIntOmelyan45
// for the 5 force version.
// 
// To construct a QPQPQPQPQ integrator, the update to the coordinate
// must be the first argument, the momentum the second.
//
// To construct a PQPQPQPQP integrator, the update to the momentum
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


AlgIntOmelyan44::AlgIntOmelyan44(AlgInt &A, AlgInt &B, 
				 IntABArg &arg_ab) 
  : AlgIntAB(A,B,arg_ab)
{

  int_type = INT_OMELYAN_44;
  A_calls = 5;
  B_calls = 4;
  rho = 0.1786178958448091;
  theta = -0.06626458266981843;
  lambda = 0.7123418310626056;
  
}

AlgIntOmelyan44::~AlgIntOmelyan44() {

}

void AlgIntOmelyan44::evolve(Float dt, int steps) 
{
  char * fname = "evolve(Float, int)";
  step_cnt = 0;
  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(step_cnt);

  A->evolve(rho*dt/(Float)A_steps, A_steps);

  for (int i=0; i<steps; i++) {
    if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    
    B->evolve(lambda*dt/(Float)B_steps, B_steps);
    A->evolve(theta*dt/(Float)A_steps, A_steps);
    B->evolve((1.0-2.0*lambda)*dt/(2.0*(Float)B_steps), B_steps);
    A->evolve((1.0-2.0*(theta+rho))*dt/(Float)A_steps, A_steps);
    B->evolve((1.0-2.0*lambda)*dt/(2.0*(Float)B_steps), B_steps);
    A->evolve(theta*dt/(Float)A_steps, A_steps);
    B->evolve(lambda*dt/(Float)B_steps, B_steps);

    if (i < steps-1) A->evolve(2.0*rho*dt/(Float)A_steps, A_steps);
    else A->evolve(rho*dt/(Float)A_steps, A_steps);
  }

  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    
}

CPS_END_NAMESPACE
