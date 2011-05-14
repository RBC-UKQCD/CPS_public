#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_int_omelyan_45.C
//
// AlgIntOmelyan45 is derived from AlgIntAB, it is an implementation
// of the 4th order Omelyan integrator using abstract operators.  This
// integrator uses 5 force evaluations per step, see AlgIntOmelyan44
// for the 4 force version.
// 
// To construct a QPQPQPQPQPQ integrator, the update to the coordinate
// must be the first argument, the momentum the second.
//
// To construct a PQPQPQPQPQP integrator, the update to the momentum
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


AlgIntOmelyan45::AlgIntOmelyan45(AlgInt &A, AlgInt &B, 
			     IntABArg &arg_ab) 
  : AlgIntAB(A,B,arg_ab)
{

  int_type = INT_OMELYAN_45;
  A_calls = 6;
  B_calls = 5;
  theta = 0.08398315262876693;
  rho = 0.2539785108410595;
  lambda = 0.6822365335719091;
  mu = -0.03230286765269967;
  
}

AlgIntOmelyan45::~AlgIntOmelyan45() {

}

void AlgIntOmelyan45::evolve(Float dt, int steps) 
{
  char * fname = "evolve(Float, int)";
  step_cnt = 0;
  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(step_cnt);

  A->evolve(theta*dt/(Float)A_steps, A_steps);

  for (int i=0; i<steps; i++) {
    if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    
    B->evolve(rho*dt/(Float)B_steps, B_steps);
    A->evolve(lambda*dt/(Float)A_steps, A_steps);
    B->evolve(mu*dt/(Float)B_steps, B_steps);
    A->evolve((1.0-2.0*(lambda+theta))*dt/(2.0*(Float)A_steps), A_steps);
    B->evolve((1.0-2.0*(mu+rho))*dt/(Float)B_steps, B_steps);
    A->evolve((1.0-2.0*(lambda+theta))*dt/(2.0*(Float)A_steps), A_steps);
    B->evolve(mu*dt/(Float)B_steps, B_steps);
    A->evolve(lambda*dt/(Float)A_steps, A_steps);
    B->evolve(rho*dt/(Float)B_steps, B_steps);

    if (i < steps-1) A->evolve(2.0*theta*dt/(Float)A_steps, A_steps);
    else A->evolve(theta*dt/(Float)A_steps, A_steps);
  }

  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    
}

CPS_END_NAMESPACE
