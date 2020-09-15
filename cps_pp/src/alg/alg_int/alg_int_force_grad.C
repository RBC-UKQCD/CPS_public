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
    if(!UniqueID()) printf("Initialized force gradient PQPQP\n");
  } else if (int_type == INT_FORCE_GRAD_QPQPQ) {
    lambda = 0.5*(1.0 - 1.0/sqrt(3.0));
    xi = 0.0;
    chi = 0.0;
    theta = (2.0-sqrt(3.0))/48.0;
    if(!UniqueID()) printf("Initialized force gradient QPQPQ\n");
  }else{
    ERR.General(cname,cname,"Unknown force gradient integrator type\n");
  }

  // save the gauge size for future reference
  Lattice & lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  g_size = GJP.VolNodeSites() * lat.GsiteSize();
  if(GJP.Gparity()) g_size *= 2;
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
#if 0
  long g_size_tlc = GJP.VolNodeSites()*18*4; if(GJP.Gparity()) g_size_tlc*=2;
  Float *test_lc = (Float*)pmalloc(g_size_tlc*sizeof(Float));
#endif

  LatticeContainer lat_cont;
  {
    Lattice & lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    lat_cont.Get(lat);
#if 0
    for(int i=0;i<g_size_tlc;i++) test_lc[i] = ((Float*)lat.GaugeField())[i];
#endif

    lat.EvolveGfield(force, 1.0, true); //For G-parity, ensure both the U and U* links are updated together
    LatticeFactory::Destroy();
  }
  sfree(cname, fname, "force", force);
  
  // do the actual evolution
  which_int->evolve(dt/(Float)steps, steps);

  // restore the gauge field
  {
    Lattice & lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    lat_cont.Set(lat);

#if 0
    for(int i=0;i<g_size_tlc;i++){
      if( test_lc[i] != ((Float*)lat.GaugeField())[i] ){
	printf("AlgIntForceGrad::evolve_fg lattice restore fail %d: %.9e %.9e\n",i,test_lc[i],((Float*)lat.GaugeField())[i]);
	exit(-1);
      }
    }
#endif

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
//    if( getVer(cname,fname)  <2) 
    A->evolve(lambda*dt/(Float)A_steps, A_steps);
    for (int i=0; i<steps; i++) {
      if (level == TOP_LEVEL_INTEGRATOR) checkpoint(cname,fname,i,steps);
      B->evolve(dt/(2.0*(Float)B_steps), B_steps);

      evolve_fg(A, 2*Chi/((1.0-2.0*lambda)*dt), (1.0-2.0*lambda)*dt);

      B->evolve(dt/(2.0*(Float)B_steps), B_steps);

      if (i < steps-1) A->evolve(2.0*lambda*dt/(Float)A_steps, A_steps);
      else A->evolve(lambda*dt/(Float)A_steps, A_steps);
      
      if (level == TOP_LEVEL_INTEGRATOR) CSM.SaveComment(++step_cnt);
    }
    break;

  case INT_FORCE_GRAD_QPQPQ:
//    if( getVer(cname,fname)  <2) 
    A->evolve(lambda*dt/(Float)A_steps, A_steps);
    for (int i=0; i<steps; i++) {
      if (level == TOP_LEVEL_INTEGRATOR) checkpoint(cname,fname,i,steps);
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

  if(GJP.Gparity()){
    /*C.Kelly 09/11:
     *For lowest level integrator, during evolution of momentum and gauge fields we can be more efficient by only updating the links and not their conjugate copies
     *(the links we pull across the boundary can be conjugated in place for very few, if any, extra flops by using alternate functions, eg. Trans rather than Dagger.)
     *However for higher level integrators involving fermion fields we need both the links and their conjugates stored, so do the copy-conjugation here.
     *Do this by calling a function which does nothing apart from on the alg_action_gauge instance, where it performs the copy-conjugation.
     */
    A->copyConjLattice();
    B->copyConjLattice();
  }

}

CPS_END_NAMESPACE
