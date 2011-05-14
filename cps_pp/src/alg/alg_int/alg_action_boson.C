#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_action_boson.C
//
// AlgActionBoson is a class describing a bilinear bosonic action
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

AlgActionBoson::AlgActionBoson(AlgMomentum &mom, ActionBosonArg &b_arg)
			       
  : AlgActionBilinear(mom, b_arg.bi_arg)
{

  cname = "AlgActionBoson";
  char *fname = "AlgActionBoson()";

  int_type = INT_BOSON;
  bsn_arg = &b_arg;

  //!< First check n_masses bilinear = n_masses boson
  if (bsn_arg->bosons.bosons_len != bsn_arg->bi_arg.bilinears.bilinears_len)
    ERR.General(cname, fname,
		"Inconsistency between BosonArg and BilinearArg n_masses");
      // ~~bsn_cp_arg local to AlgActionBoson
  if(n_masses > 0){
    //!< Allocate memory for the boson CG arguments.
    bsn_cg_arg = (CgArg **) smalloc(n_masses * sizeof(CgArg*), 
				    "bsn_cg_arg", fname, cname);
    
    for(int i=0; i<n_masses; i++){
      bsn_cg_arg[i] = (CgArg *) 
	smalloc(sizeof(CgArg), "bsn_cg_arg[i]", fname, cname);
    } 

    //!< Initialize the boson CG arguments
    for(int i=0; i<n_masses; i++){
      bsn_cg_arg[i]->mass = mass[i];
      // ~~added for twisted mass Wilson fermions
      bsn_cg_arg[i]->epsilon = bsn_arg->bosons.bosons_val[i].epsilon;
      bsn_cg_arg[i]->max_num_iter = max_num_iter[i];
      bsn_cg_arg[i]->stop_rsd = bsn_arg->bosons.bosons_val[i].stop_rsd_hb;
    }

  }

}

AlgActionBoson::~AlgActionBoson() {

  char *fname = "~AlgActionBoson()";

  //!< Free memory for the boson CG arguments
  if(n_masses > 0){
    for(int i=0; i<n_masses; i++) {
      sfree(bsn_cg_arg[i], "bsn_cg_arg[i]", fname, cname);
    }
    sfree(bsn_cg_arg, "bsn_cg_arg", fname, cname);
  }

}

//!< Heat Bath for bosons
void AlgActionBoson::heatbath() {

  char *fname = "heatbath()";

  if (n_masses > 0) {
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
    Vector *tmp1 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1",fname,cname);
    Vector *tmp2 = (Vector*)smalloc(f_size*sizeof(Float),"tmp2",cname,fname);
    
    for(int i=0; i<n_masses; i++){
      lat.RandGaussVector(tmp1, 0.5, Ncb);
      lat.RandGaussVector(phi[i], 0.5, Ncb);
      // ~~changed for twisted mass Wilson fermions; also added DAG_YES parameter
      lat.SetPhi(tmp2, tmp1, phi[i], bsn_cg_arg[i]->mass, bsn_cg_arg[i]->epsilon, DAG_YES);
      phi[i] -> VecZero(f_size);
      cg_iter = lat.FmatEvlInv(phi[i], tmp2, bsn_cg_arg[i], CNV_FRM_NO);
      
      updateCgStats(bsn_cg_arg[i]);
    }
    
    sfree(tmp2, "tmp2", cname, fname);
    sfree(tmp1, "tmp1", cname, fname);
    
    LatticeFactory::Destroy();

    traj++;
  }

}

//!< Calculate boson contribution to the Hamiltonian
Float AlgActionBoson::energy() {

  char *fname = "energy()";
  Float h = 0.0;
  
  if (n_masses > 0) {
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

    for(int i=0; i<n_masses; i++)
      // ~~changed for twisted mass Wilson fermions
      h += lat.BhamiltonNode(phi[i], bsn_cg_arg[i]->mass, bsn_cg_arg[i]->epsilon);

    LatticeFactory::Destroy();
  }

  return h;

}

void AlgActionBoson::prepare_fg(Matrix * force, Float dt_ratio)
{
  const char fname[] = "prepare_fg(M*,F)";
  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

  for (int i = 0; i<n_masses; i++) {
    //!< Need to include this hack for stag force to be correct;
    //!< implemented in F_none.C
    lat.BforceVector(phi[i], bsn_cg_arg[i]);
    // ~~changed for twisted mass Wilson fermions
    Fdt = lat.EvolveMomFforce(force, phi[i], bsn_cg_arg[i]->mass, bsn_cg_arg[i]->epsilon, -dt_ratio);
        
    if (force_measure == FORCE_MEASURE_YES) {
      char label[200];
      sprintf(label, "%s, mass = %e:", force_label, mass[i]);
      Fdt.print(dt_ratio, label);
    }
  }
  LatticeFactory::Destroy();
}

//!< run method evolves the momentum due to the boson force
void AlgActionBoson::evolve(Float dt, int nsteps) 
{
  char *fname = "evolve(Float, int)";
  if (n_masses <= 0) return;
  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
  for (int steps=0; steps<nsteps; steps++){
    for (int i = 0; i<n_masses; i++) {
      //!< Need to include this hack for stag force to be correct;
      //!< implemented in F_none.C
      lat.BforceVector(phi[i], bsn_cg_arg[i]);
      // ~~changed for twisted mass Wilson fermions
      Fdt = lat.EvolveMomFforce(mom, phi[i], bsn_cg_arg[i]->mass, bsn_cg_arg[i]->epsilon, -dt);

      if (force_measure == FORCE_MEASURE_YES) {
        char label[200];
        sprintf(label, "%s, mass = %e:", force_label, mass[i]);
        Fdt.print(dt, label);
      }
    }
  }
  
  LatticeFactory::Destroy();
}

CPS_END_NAMESPACE
