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
  int veloc_id=-1;
#if 1
 veloc_id = protect(cname,fname,mom,g_size,sizeof(Float));
#else
#ifdef HAVE_VELOC
 Float dtime = -dclock();
 VELOC_Mem_protect (  (veloc_id= VeloCCounter() ) , mom, g_size, sizeof(Float) );
 VRB.Result(cname,fname,"mom VELOC id:%d\n",veloc_id);
 dtime +=dclock();
 print_flops(fname,"VeloC()",0,dtime);
#endif
#endif
  if (veloc_id>-1) md_veloc_all.push_back(veloc_id);


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

  int veloc_v=-1;
  std::string veloc_label;
  veloc_label="MD_traj"+std::to_string(traj_num);
#if 1
  veloc_v = getVer(veloc_label.c_str());
#else
#ifdef HAVE_VELOC
   Float dtime2 = -dclock();
  {
    std::stringstream veloc_label;
    veloc_label <<"MD_traj"<<traj_num<<std::endl;
    std::string veloc_str = veloc_label.str();
    std::cout <<"Veloc label: "<<(veloc_str).c_str()<<std::endl;
    VRB.Result(cname,fname,"Veloc label version: %s %d\n",
              (veloc_str).c_str(),veloc_v);
    int veloc_v = VELOC_Restart_test((veloc_str).c_str(),65536);
    VRB.Result(cname,fname,"Veloc label version: %s %d\n",
              (veloc_str).c_str(),veloc_v);
  }
   dtime2 += dclock();
   print_flops(fname,"VeloC_Restart_test()",0,dtime);
  exit(-4);
#endif
#endif

  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  lat.RandGaussAntiHermMatrix(mom, 1.0);

  if(UniqueID()==0) printf("Heatbath for conjugate momentum\n");

  if(GJP.Gparity() && GJP.Gparity1f2fComparisonCode()){
    //For comparison with 1f approach, run RNG over second field too
    //to keep RNG sync'd over evolution
    for(int n = 0; n < GJP.VolNodeSites(); n++) {
      LRG.AssignGenerator(n,1);
      for(int j = 0; j < 4; j++) {
	for(int i = 0; i < 8; ++i) {
	  LRG.Grand(FOUR_D);
	}
      }
    }
  }

  if(GJP.Gparity1fX()){
    //doubled lattice in X-direction
    //need | P | P* |
    //or quad lattice in XY-directions (if GJP.Gparity1fY() also - note GJP.Gparity1fY() cannot return true without Gparity1fX() also true)
    // | P* | P  |
    // | P  | P* |
    if(!UniqueID()){ printf("1f G-parity: copy-conjugating momentum field\n"); fflush(stdout); }
    Lattice::CopyConjMatrixField(mom,4);
  }

#if 0
  {
    unsigned int gcsum = lat.CheckSum(mom);

    //note: for 2f G-parity the above lat.CheckSum just checksums the flavour-0 part
    //      so for correct comparison between 1f and 2f we need to do both flavours
    //      this takes extra computation so make it optional

    if(GJP.Gparity() && GJP.Gparity1f2fComparisonCode()){
      if(!UniqueID()){ printf("2f G-parity: copy-conjugating momentum field for checksum\n"); fflush(stdout); }
      Lattice::CopyConjMatrixField(mom,4);
      gcsum += lat.CheckSum(mom + 4*GJP.VolNodeSites());
    }

    QioControl qc;
    gcsum = qc.globalSumUint(gcsum);

    if(UniqueID()==0) printf("Initial conjugate momentum checksum %u\n",gcsum);
  }
#endif

  //!< reset MD time in Lattice (a momentum refresh means a new trajectory)
  lat.MdTime(0.0);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
  h_init = lat.MomHamiltonNode(mom);
  Float h= h_init;
  glb_sum(&h);
  VRB.Result(cname, fname, "Initial ham = %0.16e\n", h);
      
  LatticeFactory::Destroy();
  
  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
}

//mom = (Matrix*)smalloc(g_size*sizeof(Float),"mom",fname,cname);
void AlgMomentum::SaveState(std::string name){
  const char *fname = "SaveState()";
  VRB.Func(cname,fname);
    std::string dirname=name;
    CPS_NAMESPACE::sync();
//    if(!UniqueID())
    if(mkdir((dirname).c_str(),0777)!=0 && errno != EEXIST)
      ERR.General(cname,fname,"cannot create directory %s\n",(dirname).c_str());
    CPS_NAMESPACE::sync();
    std::string filename;
    filename = dirname+"/"+dirname;
    VRB.Result(cname,fname,"opening %s\n",(filename).c_str());
#if 1
    FILE *fp = Fopen(ADD_ID,(filename).c_str(),"w");
    Fwrite(mom,sizeof(Float),g_size,fp);
    Fclose(fp);
#endif
}

void AlgMomentum::LoadState(std::string name){
  const char *fname = "LoadState()";
  VRB.Func(cname,fname);
    std::string dirname=name;
    std::string filename;
    filename = dirname+"/"+dirname;
    VRB.Result(cname,fname,"opening %s\n",(filename).c_str());
#if 1
    FILE *fp = Fopen(ADD_ID,(filename).c_str(),"r");
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

  if(GJP.Gparity1fX() && !GJP.Gparity1fY()) h/=2;
  else if(GJP.Gparity1fX() && GJP.Gparity1fY()) h/=4;

  {
    Float gsum_h(h);
    glb_sum(&gsum_h);
    if(UniqueID()==0) printf("AlgMomentum::energy() %e\n",gsum_h);
  }

  LatticeFactory::Destroy();

  Float total_h = h;
  glb_sum(&total_h);
  VRB.Result(cname, fname, "ham = %0.16e\n", total_h);
  total_h= h - h_init;
  glb_sum(&total_h);
  VRB.Result(cname, fname, "delta_ham = %0.16e\n", total_h);

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

  {
    Float pvals[4];
    for(int ii=0;ii<4;ii++){
      int off = 18 * ii + 2;
      pvals[ii] = ((Float*)mom)[off];
    }
    if(UniqueID()==0) printf("AlgMomentum evolve conj mom Px(0) = %.9e, Py(0) = %.9e, Pz(0) = %.9e, Pt(0) = %.9e\n",pvals[0],pvals[1],pvals[2],pvals[3]);
  }    

  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  for (int i=0; i<steps; i++) lat.EvolveGfield(mom, dt);
  

  lat.MdTimeInc(dt*steps);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));

  {
    Float linkvals[4];
    for(int ii=0;ii<4;ii++){
      int off = 18 * ii;
      linkvals[ii] = ((Float*)lat.GaugeField())[off];
    }

    if(UniqueID()==0) printf("Post AlgMomentum evolve gauge links Ux(0) = %.9e, Uy(0) = %.9e, Uz(0) = %.9e, Ut(0) = %.9e\n",linkvals[0],linkvals[1],linkvals[2],linkvals[3]);
  }

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
