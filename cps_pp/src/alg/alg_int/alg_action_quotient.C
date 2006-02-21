#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_action_quotient.C
//
// AlgActionQuotient is a class describing a bilinear action, with a
// kernel given as a matrix quotient, i.e.,
//
//   S = phi^dagger M_1 (M_2^dagger M_2)^{-1} M_1^dagger phi.
//
// This action covers 2 flavour domain wall theories with the
// Pauli-Villars cancellation and also Hasenbusch type actions.  The
// mass parameters for M_1 and M_2 are specified separately.
//
// This action also uses a chronological predictor to reduce the
// number of cg iterations.
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

AlgActionQuotient::AlgActionQuotient(AlgMomentum &mom,
				   ActionQuotientArg &q_arg)
				   
  : AlgActionBilinear(mom, q_arg.bi_arg)
{

  cname = "AlgActionQuotient";
  char *fname = "AlgActionQuotient(AlgMomentum &, ActionQuotientArg &)";

  int_type = INT_QUOTIENT;
  quo_arg = &q_arg;

  //!< First check n_masses bilinear = n_masses quotient
  if (quo_arg->quotients.quotients_len != 
      quo_arg->bi_arg.bilinears.bilinears_len)
    ERR.General(cname, fname,
		"Inconsistency between QuotientArg and BilinearArg n_masses\n");

  if(n_masses > 0){
    bsn_mass = (Float*) smalloc(n_masses * sizeof(Float), 
				"bsn_mass", fname, cname);
    frm_mass = (Float*) smalloc(n_masses * sizeof(Float),
				"frm_mass", fname, cname);

    //!< Allocate memory for the CG arguments.
    bsn_cg_arg = (CgArg **) smalloc(n_masses * sizeof(CgArg*), 
				    "bsn_cg_arg", fname, cname);
    frm_cg_arg_md = (CgArg **) smalloc(n_masses * sizeof(CgArg*), 
				       "frm_cg_arg_md", fname, cname);
    frm_cg_arg_mc = (CgArg **) smalloc(n_masses * sizeof(CgArg*), 
				       "frm_cg_arg_mc", fname, cname);
    
    for(int i=0; i<n_masses; i++){
      bsn_cg_arg[i] = (CgArg *) 
	smalloc(sizeof(CgArg), "bsn_cg_arg[i]", fname, cname);
      frm_cg_arg_md[i] = (CgArg *) 
	smalloc(sizeof(CgArg), "frm_cg_arg_md[i]", fname, cname);
      frm_cg_arg_mc[i] = (CgArg *) 
	smalloc(sizeof(CgArg), "frm_cg_arg_mc[i]", fname, cname);
    } 

    //!< Initialize the CG arguments
    for(int i=0; i<n_masses; i++){
      bsn_mass[i] = quo_arg->quotients.quotients_val[i].bsn_mass;
      frm_mass[i] = quo_arg->quotients.quotients_val[i].frm_mass;

      bsn_cg_arg[i]->mass = quo_arg->quotients.quotients_val[i].bsn_mass;
      bsn_cg_arg[i]->max_num_iter = max_num_iter[i];
      bsn_cg_arg[i]->stop_rsd = quo_arg->quotients.quotients_val[i].stop_rsd_hb;

      frm_cg_arg_md[i]->mass = quo_arg->quotients.quotients_val[i].frm_mass;
      frm_cg_arg_md[i]->max_num_iter = max_num_iter[i];
      frm_cg_arg_md[i]->stop_rsd = quo_arg->quotients.quotients_val[i].stop_rsd_md;

      frm_cg_arg_mc[i]->mass = quo_arg->quotients.quotients_val[i].frm_mass;
      frm_cg_arg_mc[i]->max_num_iter = max_num_iter[i];
      frm_cg_arg_mc[i]->stop_rsd = quo_arg->quotients.quotients_val[i].stop_rsd_mc;
    }

    evolved = 1;

    //!< Copy over chronological parameters
    chrono = (int*)smalloc(n_masses*sizeof(int), "chrono", fname, cname);
    for (int i=0; i<n_masses; i++) 
      chrono[i] = quo_arg->quotients.quotients_val[i].chrono;

    //!< Vectors used to store solution history
    v = (Vector***) smalloc(n_masses*sizeof(Vector**),
			    "v", fname, cname);
    cg_sol_old = (Vector***) smalloc(n_masses*sizeof(Vector**),
				     "cg_sol_old", fname, cname);
    vm = (Vector***) smalloc(n_masses*sizeof(Vector**),
			     "vm", fname, cname);

    tmp1 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1",fname,cname);
    tmp2 = (Vector*)smalloc(f_size*sizeof(Float),"tmp2",fname,cname);
    
    for (int i=0; i<n_masses; i++) {
      int deg;
      if (chrono[i] > 0) deg = chrono[i];
      else if (chrono[i] == 0) deg = 1;
      else ERR.General(cname,fname,"Cannot have negative chronology\n");

      v[i] = (Vector**) smalloc(deg*sizeof(Vector*),"v[i]", fname, cname);
      cg_sol_old[i] = (Vector**) smalloc(deg*sizeof(Vector*),
					  "cg_sol_old[i]", fname, cname);
      vm[i] = (Vector**) smalloc(deg*sizeof(Vector*), "vm[i]", fname, cname);
      for (int j=0; j<deg; j++) {
	v[i][j] = (Vector*) smalloc(f_size*sizeof(Float),
				    "v[i][j]", fname, cname);
	vm[i][j] = (Vector*) smalloc(f_size*sizeof(Float),
				     "vm[i][j]", fname, cname);
      }
    }

  }

  init();

}

void AlgActionQuotient::init() {

  AlgActionBilinear::init();
  evolved = 1;
  for (int i=0; i<n_masses; i++) v[i][0]->VecZero(f_size);

}

AlgActionQuotient::~AlgActionQuotient() {

  char *fname = "~AlgActionQuotient()";

  if(n_masses > 0){
    //!< Free chronology
    for (int i=0; i<n_masses; i++) {
      int deg;
      if (chrono[i] > 0) deg = chrono[i];
      else if (chrono[i] == 0) deg = 1;

      for (int j=0; j<deg; j++) {
	sfree(vm[i][j],"vm[i][j]", fname, cname);
	sfree(v[i][j],"v[i][j]", fname, cname);
      }
      sfree(cg_sol_old[i], "cg_sol_old[i]", fname, cname);
      sfree(vm[i],"vm[i]", fname, cname);
      sfree(v[i],"v[i]", fname, cname);
    }
    sfree(tmp2, "tmp2", fname, cname);
    sfree(tmp1, "tmp1", fname, cname);
    sfree(cg_sol_old, "cg_sol_old", fname, cname);
    sfree(vm, "vm", fname, cname);
    sfree(v, "v", fname, cname);

    sfree(chrono,"chrono",fname,cname);

    //!< Free memory for the fermion CG arguments
    for(int i=0; i<n_masses; i++) {
      sfree(frm_cg_arg_mc[i], "frm_cg_arg_mc[i]", fname,cname);
      sfree(frm_cg_arg_md[i], "frm_cg_arg_md[i]", fname,cname);
      sfree(bsn_cg_arg[i], "bsn_cg_arg[i]", fname,cname);
    }
    sfree(frm_cg_arg_mc, "frm_cg_arg_mc", fname,cname);
    sfree(frm_cg_arg_md, "frm_cg_rg_md", fname,cname);
    sfree(bsn_cg_arg, "bsn_cg_arg", fname,cname);

    sfree(frm_mass, cname, fname, "frm_mass");
    sfree(bsn_mass, cname, fname, "bsn_mass");
  }

}

//!< Heat Bath for quotients
void AlgActionQuotient::heatbath() {

  char *fname = "heatbath()";

  if (n_masses > 0) {
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
    h_init = 0.0;

    for(int i=0; i<n_masses; i++){
      lat.RandGaussVector(tmp1, 0.5, Ncb);
      lat.RandGaussVector(tmp2, 0.5, Ncb);

      h_init += lat.SetPhi(phi[i], tmp1, tmp2, frm_mass[i], DAG_YES);

      tmp2 -> VecZero(f_size);
      cg_iter = lat.FmatEvlInv(tmp2, phi[i], bsn_cg_arg[i], CNV_FRM_NO);
      
      lat.SetPhi(phi[i], tmp2, tmp1, bsn_mass[i], DAG_NO);

      updateCgStats(bsn_cg_arg[i]);
    }
    
    LatticeFactory::Destroy();

    evolved = 0;

  }

}

//!< Calculate fermion contribution to the Hamiltonian
Float AlgActionQuotient::energy() {

  char *fname = "energy()";
  Float h = 0.0;
  
  if (n_masses > 0) {
    if (!evolved && h_init != 0.0) {
      return h_init;
    } else {
      Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);

      for(int i=0; i<n_masses; i++) {
	lat.SetPhi(tmp1, phi[i], tmp2, bsn_mass[i], DAG_YES);

	tmp2 -> VecZero(f_size);
	cg_iter = 
	  lat.FmatEvlInv(tmp2, tmp1, frm_cg_arg_mc[i], CNV_FRM_NO);
	
	updateCgStats(frm_cg_arg_mc[i]);
	  
	h += lat.FhamiltonNode(tmp1, tmp2);
      }
      
      LatticeFactory::Destroy();
    }
  }

  return h;

}

//!< run method evolves the momentum due to the fermion force
void AlgActionQuotient::evolve(Float dt, int nsteps) 
{

  char *fname = "run(Float,int)";
  int chronoDeg;
  Float Fdt;

  if (n_masses > 0){
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
    
    for (int steps=0; steps<nsteps; steps++) {
      for(int i=0; i<n_masses; i++) {

	//!< tmp1 = M^{dagger}(m_bsn) phi
	lat.SetPhi(tmp1, phi[i], tmp1, bsn_mass[i], DAG_YES);

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
	lat.FminResExt(cg_sol, tmp1, cg_sol_old[i], vm[i], 
		       chronoDeg, frm_cg_arg_md[i], CNV_FRM_NO);

	//!< cg_sol = (M_f^{dagger} M_f)^{-1} M_b^{dagger} phi
	cg_iter = 
	  lat.FmatEvlInv(cg_sol, tmp1, frm_cg_arg_md[i], CNV_FRM_NO);

	updateCgStats(frm_cg_arg_md[i]);

	int g_size = GJP.VolNodeSites() * lat.GsiteSize();
	Matrix *mom_tmp;

	if (force_measure == FORCE_MEASURE_YES) {
	  mom_tmp = (Matrix*)smalloc(g_size*sizeof(Float),"mom_tmp", fname, cname);
	  ((Vector*)mom_tmp) -> VecZero(g_size);
	} else {
	  mom_tmp = mom;
	}

	//!< Evolve mom using fermion force
	Fdt = lat.EvolveMomFforce(mom_tmp, cg_sol, frm_mass[i], dt);
	if (force_measure == FORCE_MEASURE_YES) {
	  char label[200];
	  sprintf(label, "%s (fermion), mass = %e:", force_label, frm_mass[i]);
	  printForce(Fdt, dt, label);
	}

	//!< Evolve mom using boson force
	Fdt = lat.EvolveMomFforce(mom_tmp, phi[i], cg_sol, bsn_mass[i], dt);
	if (force_measure == FORCE_MEASURE_YES) {
	  char label[200];
	  sprintf(label, "%s (boson), mass = %e:", force_label, bsn_mass[i]);
	  printForce(Fdt, dt, label);
	}

	// If measuring the force, need to measure and then sum to mom
	if (force_measure == FORCE_MEASURE_YES) {
	  Fdt = dotProduct((IFloat*)mom_tmp, (IFloat*)mom_tmp, g_size);
	  glb_sum(&Fdt);
	  fTimesV1PlusV2((IFloat*)mom, 1.0, (IFloat*)mom_tmp, 
			 (IFloat*)mom, g_size);

	  char label[200];
	  sprintf(label, "%s (total), mass = (%e,%e):", 
		  force_label, frm_mass[i], bsn_mass[i]);
	  printForce(sqrt(Fdt), dt, label);

	  sfree(mom_tmp, "mom_tmp", fname, cname);
	}

      }

      md_steps++;
    }

    LatticeFactory::Destroy();

    evolved = 1;
  }

}

CPS_END_NAMESPACE
