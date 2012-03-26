#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_int_force_grad.C
//
// AlgIntForceGrad is derived from AlgIntAB, it is an implementation of
// a fourth order force gradient integrator using abstract operators.
// 
// To construct a QPQPQ integrator, the update to the coordinate
// must be the first argument, the momentum the second.
//
// To construct a PQPQP integrator, the update to the momentum
// must be the second argument, the coordinate the first.
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
#include<util/lat_cont.h>
CPS_START_NAMESPACE


AlgIntForceGrad::AlgIntForceGrad(AlgInt &A, AlgInt &B, IntABArg &arg_ab) 
  : AlgIntAB(A,B,arg_ab)
{
  cname = "AlgIntForceGrad";

  int_type = arg_ab.type;
  A_calls = 3;
  B_calls = 2;

  // Can add more force gradient definitions if we want here
  if (int_type == INT_FORCE_GRAD_PQPQP) {
    lambda = 1.0/6.0;
    xi = 0.0;
    chi = 1.0/72.0;
    theta = 0.0;
  } else if (int_type == INT_FORCE_GRAD_QPQPQ) {
    lambda = 0.5*(1.0 - 1.0/sqrt(3.0));
    xi = 0.0;
    chi = 0.0;
    theta = (2.0-sqrt(3.0))/48.0;
  }

  // save the gauge size for future reference
  Lattice & lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  g_size = GJP.VolNodeSites() * lat.GsiteSize();
  LatticeFactory::Destroy();
}

AlgIntForceGrad::~AlgIntForceGrad() {

}

// the force gradient step, 'which_int' should be either A or B
void AlgIntForceGrad::evolve_fg(AlgInt * which_int, Float fg_dt, Float dt)
{
  const char fname[] = "evolve_fg()";
  // in principle steps here should always be 1, is there any reason
  // that people want to do multiple force gradient update sequentially?
  int steps = (which_int == A) ? A_steps : B_steps;

  Matrix * force = (Matrix*)smalloc(g_size*sizeof(Float), "force", fname, cname); 
  ((Vector*)force)->VecZero(g_size);

  which_int->prepare_fg(force, fg_dt);

  // evolve the gauge field temporarily to include the force gradient
  // contribution
  LatticeContainer lat_cont;
  {
    Lattice & lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    lat_cont.Get(lat);
    lat.EvolveGfield(force, 1.0);
    LatticeFactory::Destroy();
  }
  sfree(cname, fname, "force", force);
  
  // do the actual evolution
  which_int->evolve(dt/(Float)steps, steps);

  // restore the gauge field
  {
    Lattice & lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    lat_cont.Set(lat);
    LatticeFactory::Destroy();
  }
}

void AlgIntForceGrad::evolve(Float dt, int steps)
{
  const char fname[] = "evolve(Float, int)";

  //Float Xi = xi*dt*dt*dt;
  Float Theta = theta*dt*dt*dt;
  Float Chi = chi*dt*dt*dt;
  
  step_cnt = 0;
  if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(step_cnt);

  switch(int_type){
  case INT_FORCE_GRAD_PQPQP:
    A->evolve(lambda*dt/(Float)A_steps, A_steps);
    for (int i=0; i<steps; i++) {
      B->evolve(dt/(2.0*(Float)B_steps), B_steps);

      evolve_fg(A, 2*Chi/((1.0-2.0*lambda)*dt), (1.0-2.0*lambda)*dt);

      B->evolve(dt/(2.0*(Float)B_steps), B_steps);

      if (i < steps-1) A->evolve(2.0*lambda*dt/(Float)A_steps, A_steps);
      else A->evolve(lambda*dt/(Float)A_steps, A_steps);
      
      if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    }
    break;

  case INT_FORCE_GRAD_QPQPQ:
    A->evolve(lambda*dt/(Float)A_steps, A_steps);
    for (int i=0; i<steps; i++) {
      evolve_fg(B, 4.0*Theta/dt, dt/2.0);

      A->evolve((1.0-2.0*lambda)*dt/(Float)A_steps, A_steps);

      evolve_fg(B, 4.0*Theta/dt, dt/2.0);

      if (i < steps-1) A->evolve(2.0*lambda*dt/(Float)A_steps, A_steps);
      else A->evolve(lambda*dt/(Float)A_steps, A_steps);
      
      if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    }
    break;

  default:
    ERR.NotImplemented(cname, fname);
  }
}

CPS_END_NAMESPACE
