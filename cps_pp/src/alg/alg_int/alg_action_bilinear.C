#include<config.h>
#include<math.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_action_bilinear.C
//
// AlgAction is a derived class which defines methods common to all
// bilinear actions.
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
#include<alg/alg_remez.h>
CPS_START_NAMESPACE

//!< Dummy constructor - does nothing
AlgActionBilinear::AlgActionBilinear() 
  : AlgAction() 
{
}

AlgActionBilinear::AlgActionBilinear(AlgMomentum &mom,
				     ActionBilinearArg &b_arg)
				     
  : AlgAction(mom, b_arg.action_arg) {

  cname = "AlgActionBiliniear";
  char *fname="AlgActionBilinear(FclassType, int, HmdArg*, Matrix*)";
  VRB.Func(cname,fname);

  //!< First copy required instance parameters
  bi_arg = &b_arg;

  //!< Initialize the number of masses
  //n_masses = bi_arg->n_masses;
  n_masses = bi_arg->bilinears.bilinears_len;
  fermion = bi_arg->fermion;

  if(n_masses > MAX_HMD_MASSES){
    ERR.General(cname,fname,
		"n_masses = %d is larger than MAX_HMD_MASSES = %d\n",
		n_masses, MAX_HMD_MASSES);
  }


  if (n_masses > 0) {
    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(bi_arg->fermion, G_CLASS_NONE);
    
    //!< Number of lattice sites
    f_sites = GJP.SnodeSites()*GJP.VolNodeSites() / (lat.FchkbEvl()+1);
    //!< Number of Vectors in a Vector array
    f_vec_count = f_sites * lat.SpinComponents();
    //!< Number of Floats in a Vector array
    f_size = f_vec_count * lat.Colors() * 2;

    if(lat.FchkbEvl() == 1) Ncb = 1;      //!< Half Checkerboard
    else if(lat.FchkbEvl() == 0) Ncb = 2; //!< Full Checkerboard
    
    LatticeFactory::Destroy();
    
    //!< Allocate memory for the phi field array
    phi = (Vector **) smalloc(n_masses * sizeof(Vector*),
			      "phi",fname,cname);
    
    for(int i=0; i<n_masses; i++) {
      phi[i] = (Vector *) smalloc(f_size*sizeof(Float),"phi[i]",fname,cname);
    }

    //! Copy over mass and max iteration parameters
    mass = (Float*)smalloc(n_masses * sizeof(Float), "mass", fname, cname);
    max_num_iter = (int*)smalloc(n_masses * sizeof(Float), 
				 "max_num_iter", fname, cname);
    for(int i=0; i<n_masses; i++) {
      mass[i] = bi_arg->bilinears.bilinears_val[i].mass;
      max_num_iter[i] = bi_arg->bilinears.bilinears_val[i].max_num_iter;
    }

  }

  init();
  VRB.FuncEnd(cname,fname);
}

void AlgActionBilinear::init() {

  cg_stats.init();
  md_steps = 0;

}

AlgActionBilinear::~AlgActionBilinear() {

  char *fname="~AlgActionBilinear()";

  if (n_masses > 0) {
    //!< Free memory
    sfree(mass, cname, fname, "mass");
    sfree(max_num_iter, cname, fname, "max_num_iter");
    for(int i=0; i<n_masses; i++) 
      sfree(phi[i], cname,fname, "phi[i]");
    sfree(phi, cname,fname, "phi");
  }

}

void AlgActionBilinear::cost(CgStats *cg_stats_global) {

  cg_stats_global -> cg_calls += cg_stats.cg_calls;

  cg_stats_global->cg_iter_total += cg_stats.cg_iter_total;
  if (cg_stats.cg_iter_min < cg_stats_global->cg_iter_min)
    cg_stats_global->cg_iter_min = cg_stats.cg_iter_min;
  if (cg_stats.cg_iter_max > cg_stats_global->cg_iter_max)
    cg_stats_global->cg_iter_max = cg_stats.cg_iter_max;
  if (cg_stats_global->cg_calls > 0)
    cg_stats_global->cg_iter_av = (cg_stats_global->cg_iter_total) / 
      (Float)(cg_stats_global->cg_calls);
  else
    cg_stats_global->cg_iter_av = 0.0;

  cg_stats_global->true_rsd_total += cg_stats.true_rsd_total;
  if (cg_stats.true_rsd_min < cg_stats_global->true_rsd_min)
    cg_stats_global->true_rsd_min = cg_stats.true_rsd_min;
  if (cg_stats.true_rsd_max > cg_stats_global->true_rsd_max)
    cg_stats_global->true_rsd_max = cg_stats.true_rsd_max;
  if (cg_stats_global->cg_calls > 0)
    cg_stats_global->true_rsd_av = (cg_stats_global->true_rsd_total) / 
      (Float)(cg_stats_global->cg_calls);
  else
    cg_stats_global->true_rsd_av = 0.0;

}

void AlgActionBilinear::updateCgStats(CgArg *cg_arg) {

  cg_stats.cg_calls++;      
  cg_stats.cg_iter_total += cg_iter;
  if(cg_iter < cg_stats.cg_iter_min) cg_stats.cg_iter_min = cg_iter;
  if(cg_iter > cg_stats.cg_iter_max) cg_stats.cg_iter_max = cg_iter;
  if (cg_stats.cg_calls > 0)
    cg_stats.cg_iter_av = (cg_stats.cg_iter_total) / 
      (Float)(cg_stats.cg_calls);
  else
    cg_stats.cg_iter_av = 0.0;

  cg_stats.true_rsd_total += cg_arg->true_rsd;
  if(cg_arg->true_rsd < cg_stats.true_rsd_min) 
    cg_stats.true_rsd_min = cg_arg->true_rsd;
  if(cg_arg->true_rsd > cg_stats.true_rsd_max) 
    cg_stats.true_rsd_max = cg_arg->true_rsd;
  if (cg_stats.cg_calls > 0)
    cg_stats.true_rsd_av = (cg_stats.true_rsd_total) / 
    (Float)(cg_stats.cg_calls);
  else
    cg_stats.true_rsd_av = 0.0;

}

int AlgActionBilinear::getNmass(){
  return n_masses;
}

Float AlgActionBilinear::getMass(int i){
  return mass[i];
}

FclassType AlgActionBilinear::getFermion(){
  return bi_arg->fermion;
}

CPS_END_NAMESPACE
