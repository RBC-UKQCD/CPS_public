#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_action_fermion.C
//
// AlgActionFermion is a class describing a bilinear fermionic action
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
#include<util/dirac_op.h>
CPS_START_NAMESPACE

AlgActionFermion::AlgActionFermion(AlgMomentum &mom,
				   ActionFermionArg &f_arg)
				   
  : AlgActionBilinear(mom, f_arg.bi_arg)
{

  cname = "AlgActionFermion";
  char *fname = "AlgActionFermion(M*, L&, HmdArg*)";

  frm_arg = &f_arg;

  //!< First check n_masses bilinear = n_masses boson
  if (frm_arg->fermions.fermions_len != 
      frm_arg->bi_arg.bilinears.bilinears_len)
    ERR.General(cname, fname,
		"Inconsistency between BosonArg and BilinearArg n_masses\n");

  if(n_masses > 0){
    //!< Allocate memory for the fermion CG arguments.
    frm_cg_arg_md = (CgArg **) smalloc(n_masses * sizeof(CgArg*), 
				       cname, fname, "frm_cg_arg_md");
    frm_cg_arg_mc = (CgArg **) smalloc(n_masses * sizeof(CgArg*), 
				       cname, fname, "frm_cg_arg_mc");
    
    for(int i=0; i<n_masses; i++){
      frm_cg_arg_md[i] = (CgArg *) 
	smalloc(sizeof(CgArg), cname, fname, "frm_cg_arg_md[i]");
      frm_cg_arg_mc[i] = (CgArg *) 
	smalloc(sizeof(CgArg), cname, fname, "frm_cg_arg_mc[i]");
    } 

    //!< Initialize the fermion CG arguments
    for(int i=0; i<n_masses; i++){
      frm_cg_arg_md[i]->mass = mass[i];
      frm_cg_arg_md[i]->max_num_iter = max_num_iter[i];
      frm_cg_arg_md[i]->stop_rsd = frm_arg->fermions.fermions_val[i].stop_rsd_md;
      frm_cg_arg_mc[i]->mass = mass[i];
      frm_cg_arg_mc[i]->max_num_iter = max_num_iter[i];
      frm_cg_arg_mc[i]->stop_rsd = frm_arg->fermions.fermions_val[i].stop_rsd_mc;
    }

    evolved = 1;

    //!< Copy over chronological parameters
    chrono = (int*)smalloc(n_masses*sizeof(int), cname, fname, "chrono");
    for (int i=0; i<n_masses; i++) 
      chrono[i] = frm_arg->fermions.fermions_val[i].chrono;

    //!< Vectors used to store solution history
    v = (Vector***) smalloc(n_masses*sizeof(Vector**),
				  cname, fname, "v");
    cg_sol_old = (Vector***) smalloc(n_masses*sizeof(Vector**),
				     cname, fname, "cg_sol_old");
    vm = (Vector***) smalloc(n_masses*sizeof(Vector**),
			     cname, fname, "vm");

    for (int i=0; i<n_masses; i++) {
      int deg;
      if (chrono[i] > 0) deg = chrono[i];
      else if (chrono[i] == 0) deg = 1;
      else ERR.General(cname,fname,"Cannot have negative chronology\n");

      v[i] = (Vector**) smalloc(deg*sizeof(Vector*),
				      cname, fname, "v[i]");
      cg_sol_old[i] = (Vector**) smalloc(deg*sizeof(Vector*),
					 cname, fname, "cg_sol_old[i]");
      vm[i] = (Vector**) smalloc(deg*sizeof(Vector*), cname, fname, "vm[i]");
      for (int j=0; j<deg; j++) {
	v[i][j] = (Vector*) smalloc(f_size*sizeof(Float),
					   cname, fname, "v[i][j]");
	vm[i][j] = (Vector*) smalloc(f_size*sizeof(Float),
				     cname, fname, "vm[i][j]");
      }
    }

  }

  init();

}

void AlgActionFermion::init() {

  AlgActionBilinear::init();
  evolved = 1;
  for (int i=0; i<n_masses; i++) v[i][0]->VecZero(f_size);

}

AlgActionFermion::~AlgActionFermion() {

  char *fname = "~AlgActionFermion()";

  if(n_masses > 0){
    //!< Free chronology
    for (int i=0; i<n_masses; i++) {
      int deg;
      if (chrono[i] > 0) deg = chrono[i];
      else if (chrono[i] == 0) deg = 1;

      for (int j=0; j<deg; j++) {
	sfree(vm[i][j],cname, fname, "vm[i][j]");
	sfree(v[i][j],cname, fname, "v[i][j]");
      }
      sfree(cg_sol_old[i],cname, fname, "cg_sol_old[i]");
      sfree(vm[i],cname, fname, "vm[i]");
      sfree(v[i],cname, fname, "v[i]");
    }
    sfree(cg_sol_old,cname, fname, "cg_sol_old");
    sfree(vm,cname, fname, "vm");
    sfree(v,cname, fname, "v");

    sfree(chrono,cname,fname,"chrono");

    //!< Free memory for the fermion CG arguments
    for(int i=0; i<n_masses; i++) {
      sfree(frm_cg_arg_mc[i], cname,fname, "frm_cg_arg_mc[i]");
      sfree(frm_cg_arg_md[i], cname,fname, "frm_cg_arg_md[i]");
    }
    sfree(frm_cg_arg_mc, cname,fname, "frm_cg_arg_mc");
    sfree(frm_cg_arg_md, cname,fname, "frm_cg_arg_md");
  }

}

//!< Heat Bath for fermions
void AlgActionFermion::heatbath() {

  char *fname = "heatbath()";

  if (n_masses > 0) {
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
    Vector *tmp1 = (Vector*)smalloc(f_size*sizeof(Float),cname,fname,"tmp1");
    Vector *tmp2 = (Vector*)smalloc(f_size*sizeof(Float),cname,fname,"tmp2");
    
    h_init = 0.0;

    for(int i=0; i<n_masses; i++){
      lat.RandGaussVector(tmp1, 0.5, Ncb);
      lat.RandGaussVector(tmp2, 0.5, Ncb);
      h_init += lat.SetPhi(phi[i], tmp1, tmp2, mass[i]);
    }
    
    sfree(tmp2, cname, fname, "tmp2");
    sfree(tmp1, cname, fname, "tmp1");
    
    LatticeFactory::Destroy();

    evolved = 0;
    traj++;
  }

}

//!< Calculate fermion contribution to the Hamiltonian
Float AlgActionFermion::energy() {

  char *fname = "energy()";
  Float h = 0.0;
  
  if (n_masses > 0) {
    if (!evolved && h_init != 0.0) {
      return h_init;
    } else {
      Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

      Vector *cg_sol = 
	(Vector*)smalloc(f_size*sizeof(Float),cname,fname,"cg_sol");
      
      for(int i=0; i<n_masses; i++) {
	cg_sol -> VecZero(f_size);
	cg_iter = 
	  lat.FmatEvlInv(cg_sol, phi[i], frm_cg_arg_mc[i], CNV_FRM_NO);
	
	updateCgStats(frm_cg_arg_mc[i]);
	  
	h += lat.FhamiltonNode(phi[i], cg_sol);
      }
      
      sfree(cg_sol, cname, fname, "cg_sol");

      LatticeFactory::Destroy();
    }
  }

  return h;

}

//!< run method evolves the momentum due to the fermion force
void AlgActionFermion::evolve(Float dt, int nsteps) 
{

  char *fname = "run(Float,int)";
  int chronoDeg;

  if (n_masses > 0){
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
    for (int steps=0; steps<nsteps; steps++) {
      for(int i=0; i<n_masses; i++) {

	if (md_steps > chrono[i]) chronoDeg = chrono[i];
	else chronoDeg = md_steps;
	
	//!< Perform pointer arithmetic to avoid unnecessary copying
	int isz=0;
	if (chrono[i] > 0) isz = md_steps%chrono[i];
	cg_sol = v[i][isz];
	
	for (int j=0; j<chrono[i]; j++) {
	  int shift = isz - (j+1);
	  if (shift<0) shift += chrono[i];
	  cg_sol_old[i][j] = v[i][shift];
	}

	//!< Construct the initial guess
	lat.FminResExt(cg_sol, phi[i], cg_sol_old[i], vm[i], 
		       chronoDeg, frm_cg_arg_md[i], CNV_FRM_NO);

	cg_iter = 
	  lat.FmatEvlInv(cg_sol, phi[i], frm_cg_arg_md[i], CNV_FRM_NO);

	updateCgStats(frm_cg_arg_md[i]);

	Fdt = lat.EvolveMomFforce(mom, cg_sol, mass[i], dt);

	if (force_measure == FORCE_MEASURE_YES) {
	  sprintf(force_label, "Fermion, mass = %e:", mass[i]);
	  printForce(Fdt, dt, force_label);
	}

      }

      md_steps++;
    }

    LatticeFactory::Destroy();

    evolved = 1;
  }

}

CPS_END_NAMESPACE
