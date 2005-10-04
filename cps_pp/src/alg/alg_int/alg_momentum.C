#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_momentum.C
//
// AlgMomentum is a class which defines the conjugate momentum
// contribution to the Hamiltonian
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<alg/alg_hmd.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_int.h>
CPS_START_NAMESPACE

AlgMomentum::AlgMomentum() : AlgHamiltonian()
{
  cname = "AlgMomentum";
  char *fname = "AlgMomentum()";

  md_time_str = "MD_time/step_size = ";

  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);  
  g_size = GJP.VolNodeSites() * lat.GsiteSize();
  LatticeFactory::Destroy();

  mom = (Matrix*)smalloc(g_size*sizeof(Float),cname,fname,"mom");
}

AlgMomentum::~AlgMomentum() {
  char *fname = "~AlgMomentum()";
  sfree(mom, cname, fname, "mom");
}

//!< Heat Bath for the conjugate momentum
void AlgMomentum::heatbath() {

  char *fname = "heatbath()";
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  lat.RandGaussAntiHermMatrix(mom, 1.0);

  //!< reset MD time in Lattice (a momentum refresh means a new trajectory)
  lat.MdTime(0.0);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
      
  LatticeFactory::Destroy();
  
}

//!< Calculate gauge contribution to the Hamiltonian
Float AlgMomentum::energy() {

  char *fname = "energy()";
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  Float h = lat.MomHamiltonNode(mom);
  LatticeFactory::Destroy();
  return h;

}

//!< evolve method evolves the gauge field due to the momentum
void AlgMomentum::evolve(Float dt, int steps) 
{

  char *fname = "evolve(Float,int)";
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  for (int i=0; i<steps; i++) lat.EvolveGfield(mom, dt);
  lat.MdTimeInc(dt*steps);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
  LatticeFactory::Destroy();

}

void AlgMomentum::cost(CgStats *cg_stats_global){

}

Matrix* AlgMomentum::getMom(){
  return mom;
}

void AlgMomentum::reverse(){
  ((Vector*)mom)->VecTimesEquFloat(-1.0, g_size);
}

void AlgMomentum::init(){

}

CPS_END_NAMESPACE
