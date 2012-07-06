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
#include<alg/alg_int.h>
#include<alg/alg_hmd.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/dirac_op.h>
CPS_START_NAMESPACE

AlgActionFermion::AlgActionFermion(AlgMomentum &mom,
				   ActionFermionArg &f_arg)
				   
  : AlgActionBilinear(mom, f_arg.bi_arg)
{

  cname = "AlgActionFermion";
  char *fname = "AlgActionFermion(M*, L&, HmdArg*)";

  int_type = INT_FERMION;
  frm_arg = &f_arg;

  //!< First check n_masses bilinear = n_masses fermion
  if (frm_arg->fermions.fermions_len != 
      frm_arg->bi_arg.bilinears.bilinears_len)
    ERR.General(cname, fname,
		"Inconsistency between FermionArg and BilinearArg n_masses\n");

  if(n_masses > 0){
    //!< Allocate memory for the fermion CG arguments.
    frm_cg_arg_md = (CgArg **) smalloc(n_masses * sizeof(CgArg*), "frm_cg_arg_md", fname, cname);
    frm_cg_arg_fg = (CgArg **) smalloc(n_masses * sizeof(CgArg*), "frm_cg_arg_fg", fname, cname);
    frm_cg_arg_mc = (CgArg **) smalloc(n_masses * sizeof(CgArg*), "frm_cg_arg_mc", fname, cname);
    
    for(int i=0; i<n_masses; i++){
      frm_cg_arg_md[i] = (CgArg *) smalloc(sizeof(CgArg), "frm_cg_arg_md[i]", fname, cname);
      frm_cg_arg_fg[i] = (CgArg *) smalloc(sizeof(CgArg), "frm_cg_arg_fg[i]", fname, cname);
      frm_cg_arg_mc[i] = (CgArg *) smalloc(sizeof(CgArg), "frm_cg_arg_mc[i]", fname, cname);
    }

    //!< Initialize the fermion CG arguments
    for(int i=0; i<n_masses; i++){
      frm_cg_arg_md[i]->mass = mass[i];
      //~~ added for twisted mass Wilson fermions
      frm_cg_arg_md[i]->epsilon = frm_arg->fermions.fermions_val[i].epsilon;
      frm_cg_arg_md[i]->max_num_iter = max_num_iter[i];
      frm_cg_arg_md[i]->stop_rsd = frm_arg->fermions.fermions_val[i].stop_rsd_md;

      frm_cg_arg_fg[i]->mass = mass[i];
      //~~ added for twisted mass Wilson fermions
      frm_cg_arg_fg[i]->epsilon = frm_arg->fermions.fermions_val[i].epsilon;
      frm_cg_arg_fg[i]->max_num_iter = max_num_iter[i];
      frm_cg_arg_fg[i]->stop_rsd = frm_arg->fermions.fermions_val[i].stop_rsd_md
        *(frm_arg->fermions.fermions_val[i].stop_rsd_fg_mult);

      frm_cg_arg_mc[i]->mass = mass[i];
      //~~ added for twisted mass Wilson fermions
      frm_cg_arg_mc[i]->epsilon = frm_arg->fermions.fermions_val[i].epsilon;
      frm_cg_arg_mc[i]->max_num_iter = max_num_iter[i];
      frm_cg_arg_mc[i]->stop_rsd = frm_arg->fermions.fermions_val[i].stop_rsd_mc;
    }

    evolved = 1;

    //!< Copy over chronological parameters
    chrono = (int*)smalloc(n_masses*sizeof(int), cname, fname, "chrono");
    for (int i=0; i<n_masses; i++) 
      chrono[i] = frm_arg->fermions.fermions_val[i].chrono;

    //!< Vectors used to store solution history
    v = (Vector***) smalloc(n_masses*sizeof(Vector**), "v", fname, cname);
    cg_sol_old = (Vector***) smalloc(n_masses*sizeof(Vector**), "cg_sol_old", fname, cname);
    vm = (Vector***) smalloc(n_masses*sizeof(Vector**), "vm", fname, cname);

    for (int i=0; i<n_masses; i++) {
      int deg=0;
      if (chrono[i] > 0) deg = chrono[i];
      else if (chrono[i] == 0) deg = 1;
      else ERR.General(cname,fname,"Cannot have negative chronology\n");

      v[i] = (Vector**) smalloc(deg*sizeof(Vector*), "v[i]", fname, cname);
      cg_sol_old[i] = (Vector**) smalloc(deg*sizeof(Vector*), "cg_sol_old[i]", fname, cname);
      vm[i] = (Vector**) smalloc(deg*sizeof(Vector*), "vm[i]", fname, cname);
      for (int j=0; j<deg; j++) {
	v[i][j] = (Vector*) smalloc(f_size*sizeof(Float), "v[i][j]", fname, cname);
	vm[i][j] = (Vector*) smalloc(f_size*sizeof(Float), "vm[i][j]", fname, cname);
      }
    }

  }

  init();

  fg_forecast = false;
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
      int deg=0;
      if (chrono[i] > 0) deg = chrono[i];
      else if (chrono[i] == 0) deg = 1;

      for (int j=0; j<deg; j++) {
	sfree(vm[i][j],"vm[i][j]",fname, cname);
	sfree(v[i][j],"v[i][j]",fname, cname);
      }
      sfree(cg_sol_old[i],"cg_sol_old[i]", fname, cname);
      sfree(vm[i],"vm[i]", fname, cname);
      sfree(v[i],"v[i]", fname, cname);
    }
    sfree(cg_sol_old, "cg_sol_old", fname, cname);
    sfree(vm,"vm", fname, cname);
    sfree(v,"v", fname, cname);

    sfree(chrono,"chrono",fname,cname);

    //!< Free memory for the fermion CG arguments
    for(int i=0; i<n_masses; i++) {
      sfree(frm_cg_arg_mc[i], "frm_cg_arg_mc[i]", fname,cname);
      sfree(frm_cg_arg_md[i], "frm_cg_arg_md[i]", fname,cname);
      sfree(frm_cg_arg_fg[i], "frm_cg_arg_fg[i]", fname,cname);
    }
    sfree(frm_cg_arg_mc, "frm_cg_arg_mc", fname,cname);
    sfree(frm_cg_arg_md, "frm_cg_arg_md", fname,cname);
    sfree(frm_cg_arg_fg, "frm_cg_arg_fg", fname,cname);
  }

}

//!< Heat Bath for fermions
void AlgActionFermion::heatbath() {

  char *fname = "heatbath()";

  if (n_masses > 0) {
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
    Vector *tmp1 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1",fname,cname);
    Vector *tmp2 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1",fname,cname);
    
    h_init = 0.0;

    for(int i=0; i<n_masses; i++){
      lat.RandGaussVector(tmp1, 0.5, Ncb);
      lat.RandGaussVector(tmp2, 0.5, Ncb);
      //~~ changed for twisted mass Wilson fermions; also added DAG_YES parameter
      h_init += (lat.Fclass() == F_CLASS_WILSON_TM) ?
			lat.SetPhi(phi[i], tmp1, tmp2, frm_cg_arg_mc[i]->mass, frm_cg_arg_mc[i]->epsilon, DAG_YES) :
   			lat.SetPhi(phi[i], tmp1, tmp2, frm_cg_arg_mc[i]->mass, DAG_YES);
    }
    
    sfree(tmp2, "tmp2", fname, cname);
    sfree(tmp1, "tmp1", fname, cname);
    
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
                (Vector*)smalloc(f_size*sizeof(Float),"cg_sol",fname,cname);
      
            for(int i=0; i<n_masses; i++) {
                cg_sol -> VecZero(f_size);
                cg_iter = 
                    lat.FmatEvlInv(cg_sol, phi[i], frm_cg_arg_mc[i], CNV_FRM_NO);
	
                updateCgStats(frm_cg_arg_mc[i]);
	  
                h += lat.FhamiltonNode(phi[i], cg_sol);
            }
      
            sfree(cg_sol, "cg_sol", cname, fname);
            LatticeFactory::Destroy();
        }
    }

    return h;

}


void AlgActionFermion::prepare_fg(Matrix * force, Float dt_ratio)
{
  const char fname[] = "prepare_fg(M*,F)";
  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

  int chronoDeg;
  for(int i=0; i<n_masses; i++) {
    chronoDeg = ( md_steps > chrono[i] ) ? chrono[i] : md_steps ;

    //!< Perform pointer arithmetic to avoid unnecessary copying
    int isz = ( chrono[i] > 0 ) ? ( md_steps % chrono[i] ) : 0 ;
    cg_sol = v[i][isz];
	
    for (int j=0; j<chrono[i]; j++) {
      int shift = isz - (j+1);
      if (shift<0) shift += chrono[i];
      cg_sol_old[i][j] = v[i][shift];
    }

    //!< Construct the initial guess
    lat.FminResExt(cg_sol, phi[i], cg_sol_old[i], vm[i], chronoDeg, frm_cg_arg_fg[i], CNV_FRM_NO);
    cg_iter = lat.FmatEvlInv(cg_sol, phi[i], frm_cg_arg_fg[i], CNV_FRM_NO);

    updateCgStats(frm_cg_arg_fg[i]);

    //~~ changed for twisted mass Wilson fermions; also added DAG_YES parameter
    Fdt = (lat.Fclass() == F_CLASS_WILSON_TM) ?
      lat.EvolveMomFforce(force, cg_sol, mass[i], frm_cg_arg_mc[i]->epsilon, dt_ratio) :
      lat.EvolveMomFforce(force, cg_sol, mass[i], dt_ratio);

    if (force_measure == FORCE_MEASURE_YES) {
      char label[200];
      sprintf(label, "%s, mass = %e:", force_label, mass[i]);
      Fdt.print(dt_ratio, label);
    }
  }
  // We now have a solution to forecast the next normal solve.
  fg_forecast = true;
  
  md_steps++;
  LatticeFactory::Destroy();
}

//!< run method evolves the momentum due to the fermion force
void AlgActionFermion::evolve(Float dt, int nsteps)
{
  char *fname = "evolve(Float, int)";
  
  int chronoDeg;
  if (n_masses <= 0) return;
  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

  for (int steps=0; steps<nsteps; steps++) {
    for(int i=0; i<n_masses; i++) {

      chronoDeg = ( md_steps > chrono[i] ) ? chrono[i] : md_steps ;

      //!< Perform pointer arithmetic to avoid unnecessary copying
      int isz = ( chrono[i] > 0 ) ? ( md_steps % chrono[i] ) : 0 ;
      cg_sol = v[i][isz];
	
      for (int j=0; j<chrono[i]; j++) {
        int shift = isz - (j+1);
        if (shift<0) shift += chrono[i];
        cg_sol_old[i][j] = v[i][shift];
      }

      //!< Construct the initial guess
      //
      // Branch condition added by Hantao: only if we are NOT using
      // chronological inverter as well as we have a previous
      // fg-solution in hand, we skip the function FminResExt().
      //
      // If chronoDeg == 0 and we don't have a fg forecast, we still
      // want this function to zero the initial guess for us.
      if (fg_forecast == false || chronoDeg != 0) {
        lat.FminResExt(cg_sol, phi[i], cg_sol_old[i], vm[i], chronoDeg, frm_cg_arg_md[i], CNV_FRM_NO);
      }else{
        VRB.Result(cname, fname, "Using force gradient forecasting.\n");
      }
      cg_iter = lat.FmatEvlInv(cg_sol, phi[i], frm_cg_arg_md[i], CNV_FRM_NO);

      updateCgStats(frm_cg_arg_md[i]);

      //~~ changed for twisted mass Wilson fermions; also added DAG_YES parameter
      Fdt = (lat.Fclass() == F_CLASS_WILSON_TM) ?
        lat.EvolveMomFforce(mom, cg_sol, mass[i], frm_cg_arg_mc[i]->epsilon, dt) :
        lat.EvolveMomFforce(mom, cg_sol, mass[i], dt);

      if (force_measure == FORCE_MEASURE_YES) {
        char label[200];
        sprintf(label, "%s, mass = %e:", force_label, mass[i]);
        Fdt.print(dt, label);
      }

    }
    // Note that as long as the last solve in a trajectory is NOT a
    // force gradient solve (which should always be the case), we
    // won't bring our solution across trajectories by using the
    // statement here.
    fg_forecast = false;

    md_steps++;
  }

  LatticeFactory::Destroy();
  evolved = 1;
}

CPS_END_NAMESPACE
