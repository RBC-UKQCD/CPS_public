#include<config.h>
#include<string.h>
#include<string>
#include<sys/stat.h>
#include<errno.h>
#include<stdio.h>
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
#include<util/time_cps.h>
#include<util/timer.h>
#include<alg/alg_int.h>
CPS_START_NAMESPACE

AlgMomentum::AlgMomentum() : AlgHamiltonian()
{
  cname = "AlgMomentum";
  const char *fname = "AlgMomentum()";

  int_type = INT_MOM;
  md_time_str = "MD_time/step_size = ";

  mom = (Matrix*)smalloc(g_size*sizeof(Float),"mom",fname,cname);
#ifdef HAVE_VELOC
 int veloc_id;
 VELOC_Mem_protect (  (veloc_id= VeloCCounter() ) , mom, g_size, sizeof(Float) );
 VRB.Result(cname,fname,"mom VELOC id:%d\n",veloc_id);
#endif
}

AlgMomentum::~AlgMomentum() {
  const char *fname = "~AlgMomentum()";
  sfree(mom, "mom", fname, cname);
}

//!< Heat Bath for the conjugate momentum
void AlgMomentum::heatbath() {

  const char *fname = "heatbath()";
  Float dtime = -dclock();
  VRB.Func(cname, fname);

  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  lat.RandGaussAntiHermMatrix(mom, 1.0);

  //!< reset MD time in Lattice (a momentum refresh means a new trajectory)
  lat.MdTime(0.0);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
      
  LatticeFactory::Destroy();
  
  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
}

//mom = (Matrix*)smalloc(g_size*sizeof(Float),"mom",fname,cname);
void AlgMomentum::SaveState(std::string name){
  const char *fname = "SaveState()";
  VRB.Func(cname,fname);
    std::stringstream dirname;
    dirname << name;
    CPS_NAMESPACE::sync();
//    if(!UniqueID())
    if(mkdir((dirname.str()).c_str(),0777)!=0 && errno != EEXIST)
      ERR.General(cname,fname,"cannot create directory %s\n",(dirname.str()).c_str());
    CPS_NAMESPACE::sync();
    //ugly, but C++ file io is very slow on some systems
    std::stringstream filename;
    filename <<dirname.str()<<"/"<<dirname.str();
    VRB.Result(cname,fname,"opening %s\n",(filename.str()).c_str());
#if 1
    FILE *fp = Fopen(ADD_ID,(filename.str()).c_str(),"w");
    Fwrite(mom,sizeof(Float),g_size,fp);
    Fclose(fp);
#endif
}

void AlgMomentum::LoadState(std::string name){
  const char *fname = "LoadState()";
  VRB.Func(cname,fname);
    std::stringstream dirname;
    dirname << name;
    //ugly, but C++ file io is very slow on some systems
    std::stringstream filename;
    filename <<dirname.str()<<"/"<<dirname.str();
    VRB.Result(cname,fname,"opening %s\n",(filename.str()).c_str());
#if 1
    FILE *fp = Fopen(ADD_ID,(filename.str()).c_str(),"r");
    Fread(mom,sizeof(Float),g_size,fp);
    Fclose(fp);
#endif
}

//!< Calculate gauge contribution to the Hamiltonian
Float AlgMomentum::energy() {

  const char *fname = "energy()";
  VRB.Func(cname, fname);
  static Timer time(cname, fname);
  time.start(true);

  Float dtime = -dclock();
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  Float h = lat.MomHamiltonNode(mom);
  LatticeFactory::Destroy();

  Float total_h = h;
  glb_sum(&total_h);
  VRB.Result(cname, fname, "ham = %0.16e\n", total_h);

  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
  time.stop(true);

  return h;
}

//!< evolve method evolves the gauge field due to the momentum
void AlgMomentum::evolve(Float dt, int steps) 
{
  const char *fname = "evolve()";
  Float dtime = -dclock();

  VRB.Func(cname, fname);

  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  for (int i=0; i<steps; i++) lat.EvolveGfield(mom, dt);
  lat.MdTimeInc(dt*steps);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
  LatticeFactory::Destroy();

  dtime += dclock();
  print_flops(cname, fname, 1968. * 4. * GJP.VolNodeSites() * steps, dtime);
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
