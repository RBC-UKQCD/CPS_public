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
  bsn_arg = &b_arg;

  //!< First check n_masses bilinear = n_masses boson
  if (bsn_arg->bosons.bosons_len != bsn_arg->bi_arg.bilinears.bilinears_len)
    ERR.General(cname, fname,
		"Inconsistency between BosonArg and BilinearArg n_masses");

  if(n_masses > 0){
    //!< Allocate memory for the boson CG arguments.
    bsn_cg_arg = (CgArg **) smalloc(n_masses * sizeof(CgArg*), 
				    cname, fname, "bsn_cg_arg");
    
    for(int i=0; i<n_masses; i++){
      bsn_cg_arg[i] = (CgArg *) 
	smalloc(sizeof(CgArg), cname, fname, "bsn_cg_arg[i]");
    } 

    //!< Initialize the boson CG arguments
    for(int i=0; i<n_masses; i++){
      bsn_cg_arg[i]->mass = mass[i];
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
      sfree(bsn_cg_arg[i], cname,fname, "bsn_cg_arg[i]");
    }
    sfree(bsn_cg_arg, cname,fname, "bsn_cg_arg");
  }

}

//!< Heat Bath for bosons
void AlgActionBoson::heatbath() {

  char *fname = "heatbath()";

  if (n_masses > 0) {
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
    Vector *tmp1 = (Vector*)smalloc(f_size*sizeof(Float),cname,fname,"tmp1");
    Vector *tmp2 = (Vector*)smalloc(f_size*sizeof(Float),cname,fname,"tmp2");
    
    for(int i=0; i<n_masses; i++){
      lat.RandGaussVector(tmp1, 0.5, Ncb);
      lat.RandGaussVector(phi[i], 0.5, Ncb);
      lat.SetPhi(tmp2, tmp1, phi[i], mass[i]);
      phi[i] -> VecZero(f_size);
      cg_iter = lat.FmatEvlInv(phi[i], tmp2, bsn_cg_arg[i], CNV_FRM_NO);
      
      updateCgStats(bsn_cg_arg[i]);
    }
    
    sfree(tmp2, cname, fname, "tmp2");
    sfree(tmp1, cname, fname, "tmp1");
    
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
      h += lat.BhamiltonNode(phi[i], mass[i]);

    LatticeFactory::Destroy();
  }

  return h;

}

//!< run method evolves the momentum due to the boson force
void AlgActionBoson::evolve(Float dt, int nsteps) 
{

  char *fname = "run(Float,int)";

  if (n_masses > 0){
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
    for (int steps=0; steps<nsteps; steps++) 
      for (int i = 0; i<n_masses; i++) {
	//!< Need to include this hack for stag force to be correct
	lat.BforceVector(phi[i], bsn_cg_arg[i]);
	Fdt = lat.EvolveMomFforce(mom, phi[i], mass[i], -dt);

	if (force_measure == FORCE_MEASURE_YES) {
	  sprintf(force_label, "Boson, mass = %e:", mass[i]);
	  printForce(Fdt, dt, force_label);
	}
      }
    
    LatticeFactory::Destroy();
  }

}

CPS_END_NAMESPACE
