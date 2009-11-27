#include<config.h>
#include<math.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_action_rational.C
//
// AlgActionRationalQuotient is a bilinear action where the matrix is represented
// by a rational approximation.
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

//!< Dummy contructor - does nothing
AlgActionRationalQuotient::AlgActionRationalQuotient()
  : AlgActionRational()
{

}

//!< Need to add restart method or similar so do not need to
//reconstruct actions - add calls to restart methods to hmc class

AlgActionRationalQuotient::AlgActionRationalQuotient(AlgMomentum &mom,
				     ActionRationalQuotientArg &r_arg, int traj_num)
				     
  : AlgActionRational(mom, r_arg.bi_arg)
{

  cname = "AlgActionRationalQuotient";
  char *fname = "AlgActionRationalQuotient()";

  int_type = INT_RATIONAL_QUOTIENT;
  rat_quo_arg = &r_arg;

  //!< First check n_masses bosons = n_masses fermions
  if (rat_quo_arg->bosons.bosons_len != rat_quo_arg->fermions.fermions_len)
    ERR.General(cname, fname,
	      "Inconsistency between number of fermions and bosons\n");
  
  //!< Also check n_masses bilinear = n_masses fermions
  if (rat_quo_arg->bi_arg.bilinears.bilinears_len != 
      rat_quo_arg->fermions.fermions_len)
    ERR.General(cname, fname,
		"Inconsistency between number of fermions and bilinears\n");
  
  //!< Also check we have bosons and fermions correctly specified
  for (int i=0; i<n_masses; i++) {
    if (rat_quo_arg->bosons.bosons_val[i].field_type != BOSON)
      ERR.General(cname,fname,"Boson %d not set as a Boson\n", i);
    if (rat_quo_arg->fermions.fermions_val[i].field_type != FERMION)
      ERR.General(cname,fname,"Fermion %d not set as a Fermion\n", i);
  }

  // RationalQuotient force term not implemented for asqtad or p4 yet
  if( fermion == F_CLASS_ASQTAD || fermion == F_CLASS_P4)
    ERR.General(cname,fname,"Force not implemented for Fasqtad or Fp4\n");

  //!< Allocate memory for the fermion CG arguments.
  if(n_masses > 0){
    bsn_mass = (Float*) smalloc(n_masses * sizeof(Float), 
				"bsn_mass", fname, cname);
    frm_mass = (Float*) smalloc(n_masses * sizeof(Float),
				"frm_mass", fname, cname);

    for(int i=0; i<n_masses; i++) {
      bsn_mass[i] = rat_quo_arg->bsn_mass.bsn_mass_val[i];
      frm_mass[i] = rat_quo_arg->frm_mass.frm_mass_val[i];
    }
    
    //!< construct approximation if necessary
    generateApprox(frm_mass,&frm_remez_arg_md,&frm_remez_arg_mc,
		   rat_quo_arg->fermions.fermions_val);
    generateApprox(bsn_mass,&bsn_remez_arg_md,&bsn_remez_arg_mc,
		   rat_quo_arg->bosons.bosons_val);
    generateCgArg(frm_mass,&frm_cg_arg_md,&frm_cg_arg_mc,"frm_cg_arg",
		  rat_quo_arg->fermions.fermions_val);
    generateCgArg(bsn_mass,&bsn_cg_arg_md,&bsn_cg_arg_mc,"bsn_cg_arg",
		  rat_quo_arg->bosons.bosons_val);

    max_size = 0;
    for (int i=0; i<n_masses; i++) {
      int tmp = frm_remez_arg_md[i].degree + 3*bsn_remez_arg_md[i].degree;
      if (max_size < tmp) max_size = tmp;
    }
    
    //!< Allocate memory for the frmn solution fermion fields.
      
    //!< For dwf we need the solution vector contiguous in memory
    frmn = (Vector**) smalloc(max_size*sizeof(Vector*), "frmn", fname, cname);
    
    frmn[0] = (Vector*) smalloc(f_size*max_size*sizeof(Float), 
				"frmn[0]", fname, cname);
    
    for (int i=1; i<max_size; i++) frmn[i] = frmn[0] + i*f_vec_count;
    
    frmn_d = 0;

    eta = (Vector**)smalloc(2*sizeof(Vector*), "eta", fname, cname);
    eta[0] = (Vector*)smalloc(f_size*sizeof(Float), "eta[0]", fname, cname);
    eta[1] = (Vector*)smalloc(f_size*sizeof(Float), "eta[1]", fname, cname);

    all_res = (Float*)smalloc(max_size*sizeof(Float), "all_res", fname, cname);
    frmn_tmp = (Vector**)smalloc(max_size*sizeof(Vector*), "frmn_tmp", fname, cname);

  }

  init(traj_num);

  if (rat_quo_arg->eigen.eigen_measure == EIGEN_MEASURE_YES) 
    generateEigArg(rat_quo_arg->eigen);
}

AlgActionRationalQuotient::~AlgActionRationalQuotient() {

  char *fname = "~AlgActionRationalQuotient()" ;
  VRB.Func(cname,fname);

  //!< Free memory for timescale split partial fraction
  if (n_masses > 0) {
    sfree(frmn_tmp,"frmn_tmp",fname,cname);
    sfree(all_res,"all_res",fname,cname);

    sfree(eta[1],"eta[1]",fname,cname);
    sfree(eta[0],"eta[0]",fname,cname);
    sfree(eta,"eta",fname,cname);

    //!< Free memory for the frmn (pseudo fermion) solution fields.
    sfree(frmn[0], "frmn[0]", fname, cname);
    sfree(frmn, "frmn", fname, cname);

    //!< Free memory for the fermion CG arguments
    destroyCgArg(bsn_cg_arg_md, bsn_cg_arg_mc, "bsn_cg_arg", 
		 bsn_remez_arg_md, bsn_remez_arg_mc);
    destroyCgArg(frm_cg_arg_md, frm_cg_arg_mc, "frm_cg_arg", 
		 frm_remez_arg_md, frm_remez_arg_mc);

    //!< Must not free these until CG args are freed.
    destroyApprox(bsn_remez_arg_md,bsn_remez_arg_mc);
    destroyApprox(frm_remez_arg_md,frm_remez_arg_mc);

    sfree(bsn_mass,"bsn_mass",fname,cname);
    sfree(frm_mass,"frm_mass",fname,cname);

    if (rat_quo_arg->eigen.eigen_measure == EIGEN_MEASURE_YES) destroyEigArg();
  }

}

//!< Heat Bath for the pseudo-fermions (phi)
void AlgActionRationalQuotient::reweight(Float *rw_fac,Float *norm) {

  char *fname = "reweight()";


//    if (n_masses>1)
//      ERR.General(cname,fname,"not implemented for n_mass>1(%d)\n",n_masses);

    //!< Before energy is measured, do we want to check bounds?
    if (rat_quo_arg->eigen.eigen_measure == EIGEN_MEASURE_YES) {
      checkApprox(bsn_mass, bsn_remez_arg_mc, rat_quo_arg->eigen);
      checkApprox(frm_mass, frm_remez_arg_mc, rat_quo_arg->eigen);
    }

    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  
    
    for(int i=0; i<n_masses; i++){

      lat.RandGaussVector(phi[i], 0.5, Ncb);
      norm[i] = lat.FhamiltonNode(phi[i],phi[i]);

      //!< First apply boson rational
      frmn[0] -> VecEqualsVecTimesEquFloat(phi[i],bsn_remez_arg_mc[i].norm,
					   f_size);
      
      cg_iter = lat.FmatEvlMInv(frmn, phi[i], bsn_remez_arg_mc[i].pole, 
				bsn_remez_arg_mc[i].degree, 0, 
				bsn_cg_arg_mc[i], CNV_FRM_NO, SINGLE, 
				bsn_remez_arg_mc[i].residue);
      
      //!< Now apply fermion rational
      frmn[1] -> VecEqualsVecTimesEquFloat(frmn[0],frm_remez_arg_mc[i].norm,
					   f_size);
      
      cg_iter = lat.FmatEvlMInv(frmn+1, frmn[0], frm_remez_arg_mc[i].pole, 
				frm_remez_arg_mc[i].degree, 0, 
				frm_cg_arg_mc[i], CNV_FRM_NO, SINGLE, 
				frm_remez_arg_mc[i].residue);
      
      updateCgStats(bsn_cg_arg_mc[i][0]);
      updateCgStats(frm_cg_arg_mc[i][0]);

      // shift this evaluation into minvcg?
      rw_fac[i] = lat.FhamiltonNode(frmn[1], frmn[1]);
      rw_fac[i] -= norm[i];
      glb_sum(rw_fac+i);
      glb_sum(norm+i);

    }

    LatticeFactory::Destroy();

//    evolved = 0;
//    heatbathEval = 1;
//    energyEval = 0;
//    traj++;

}

//!< Heat Bath for the pseudo-fermions (phi)
void AlgActionRationalQuotient::heatbath() {

  char *fname = "heatbath()";

  //!< Only evaluate heatbath if necessary
  if (!heatbathEval) {

    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  
    h_init = 0.0;
    
    for(int i=0; i<n_masses; i++){

      lat.RandGaussVector(phi[i], 0.5, Ncb);
      h_init += lat.FhamiltonNode(phi[i],phi[i]);

      //!< First apply the fermion rational
      frmn[0] -> 
	VecEqualsVecTimesEquFloat(phi[i],frm_remez_arg_mc[i].norm_inv,f_size);
      
      cg_iter = lat.FmatEvlMInv(frmn, phi[i], frm_remez_arg_mc[i].pole_inv, 
				frm_remez_arg_mc[i].degree, 0, 
				frm_cg_arg_mc[i], CNV_FRM_NO, SINGLE,
				frm_remez_arg_mc[i].residue_inv);
      
      //!< Now apply the boson rational
      phi[i] -> 
	VecEqualsVecTimesEquFloat(frmn[0],bsn_remez_arg_mc[i].norm_inv,f_size);

      cg_iter = lat.FmatEvlMInv(phi+i, frmn[0], bsn_remez_arg_mc[i].pole_inv, 
				bsn_remez_arg_mc[i].degree, 0, 
				bsn_cg_arg_mc[i], CNV_FRM_NO, SINGLE,
				bsn_remez_arg_mc[i].residue_inv);
      
      updateCgStats(frm_cg_arg_mc[i][0]);
      updateCgStats(bsn_cg_arg_mc[i][0]);
    }

    LatticeFactory::Destroy();

    evolved = 0;
    heatbathEval = 1;
    energyEval = 0;
    traj++;
  }

}

// Calculate rhmc fermion contribution to the Hamiltonian
Float AlgActionRationalQuotient::energy() {

  char *fname="energy()";

  if (energyEval) {
    return 0.0;
  } else if (!evolved) {
    energyEval = 1;
    return h_init;
  } else {
    int shift = 0;
    Float h = 0.0;

    //!< Before energy is measured, do we want to check bounds?
    if (rat_quo_arg->eigen.eigen_measure == EIGEN_MEASURE_YES) {
      checkApprox(bsn_mass, bsn_remez_arg_mc, rat_quo_arg->eigen);
      checkApprox(frm_mass, frm_remez_arg_mc, rat_quo_arg->eigen);
    }

    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  

    for (int i=0; i<n_masses; i++) {

      //!< First apply boson rational
      frmn[0] -> VecEqualsVecTimesEquFloat(phi[i],bsn_remez_arg_mc[i].norm,
					   f_size);
      
      cg_iter = lat.FmatEvlMInv(frmn, phi[i], bsn_remez_arg_mc[i].pole, 
				bsn_remez_arg_mc[i].degree, 0, 
				bsn_cg_arg_mc[i], CNV_FRM_NO, SINGLE, 
				bsn_remez_arg_mc[i].residue);
      
      //!< Now apply fermion rational
      frmn[1] -> VecEqualsVecTimesEquFloat(frmn[0],frm_remez_arg_mc[i].norm,
					   f_size);
      
      cg_iter = lat.FmatEvlMInv(frmn+1, frmn[0], frm_remez_arg_mc[i].pole, 
				frm_remez_arg_mc[i].degree, 0, 
				frm_cg_arg_mc[i], CNV_FRM_NO, SINGLE, 
				frm_remez_arg_mc[i].residue);
      
      updateCgStats(bsn_cg_arg_mc[i][0]);
      updateCgStats(frm_cg_arg_mc[i][0]);

      // shift this evaluation into minvcg?
      h += lat.FhamiltonNode(frmn[1], frmn[1]);
    }

    LatticeFactory::Destroy();

    energyEval = 1;

    return h;
  }

  
}

//!< run method evolves the integrator
void AlgActionRationalQuotient::evolve(Float dt, int nsteps)
{

  char *fname = "evolve(Float,int)";
  int shift = 0;
  int isz = 0;

  if (n_masses > 0) {
    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  
    
    for(int steps = 0; steps<nsteps; steps++) {
      for(int i=0; i<n_masses; i++){	

	int bsn_deg = bsn_remez_arg_md[i].degree;
	int frm_deg = frm_remez_arg_md[i].degree;

	//! First apply boson rational
	shift = 0;
	cg_iter = lat.FmatEvlMInv(frmn, phi[i], 
				  bsn_remez_arg_md[i].pole+isz, 
				  bsn_deg, isz, 
				  bsn_cg_arg_md[i]+isz, CNV_FRM_NO, 
				  frmn_d+shift);	

	updateCgStats(bsn_cg_arg_md[i][isz]);

	//!< Construct rhs
	eta[0] -> VecEqualsVecTimesEquFloat(phi[i],bsn_remez_arg_md[i].norm,
					    f_size);
	for (int j=0; j<bsn_deg; j++)
	  eta[0] -> FTimesV1PlusV2(bsn_remez_arg_md[i].residue[j],
				   frmn[j],eta[0],f_size);

	//!< Now apply fermion rational
	shift += bsn_deg;
	cg_iter = lat.FmatEvlMInv(frmn+shift, eta[0], 
				  frm_remez_arg_md[i].pole+isz, 
				  frm_deg, isz, 
				  frm_cg_arg_md[i]+isz, CNV_FRM_NO, 
				  frmn_d+shift+isz);	

	updateCgStats(frm_cg_arg_md[i][isz]);
	
	//!< Construct rhs
	eta[1] -> VecEqualsVecTimesEquFloat(eta[0],frm_remez_arg_md[i].norm,
					    f_size);
	for (int j=0; j<frm_deg; j++)
	  eta[1] -> FTimesV1PlusV2(frm_remez_arg_md[i].residue[j],
				   frmn[j+shift],eta[1],f_size);

	//!< Apply final boson rational
	shift += frm_deg;
	cg_iter = lat.FmatEvlMInv(frmn+shift, eta[1], 
				  bsn_remez_arg_md[i].pole+isz, 
				  bsn_deg, isz, 
				  bsn_cg_arg_md[i]+isz, CNV_FRM_NO, 
				  frmn_d+shift);	
	updateCgStats(bsn_cg_arg_md[i][isz]);

	//!< Now construct additional vectors needed
	for (int j=0; j<bsn_deg; j++) {
      Float one=1.;
	  frmn[j+shift+bsn_deg] -> 
	    FTimesV1PlusV2(one, frmn[j], frmn[j+shift], f_size);
	}

	//!< Copy over required residues and setup pointers for bosonic force
	for (int j=0; j<bsn_deg; j++) {
	  all_res[j] = -bsn_remez_arg_md[i].residue[j];
	  all_res[j+bsn_deg] = -bsn_remez_arg_md[i].residue[j];
	  all_res[j+2*bsn_deg] = bsn_remez_arg_md[i].residue[j];

	  frmn_tmp[j] = frmn[j];
	  frmn_tmp[j+bsn_deg] = frmn[j+bsn_deg+frm_deg];
	  frmn_tmp[j+2*bsn_deg] = frmn[j+2*bsn_deg+frm_deg];
	}

	Matrix *mom_tmp;
	if (force_measure == FORCE_MEASURE_YES) {
	  mom_tmp = (Matrix*)smalloc(g_size*sizeof(Float),cname, fname, "mom_tmp");
	  ((Vector*)mom_tmp)->VecZero(g_size);
	} else {
	  mom_tmp = mom;
	}

	//!< Do bosonic force contribution
	Fdt = lat.RHMC_EvolveMomFforce(mom_tmp, frmn_tmp, 3*bsn_deg, 0,
				       all_res, bsn_mass[i], dt, frmn_d, 
				       force_measure);
	if (force_measure == FORCE_MEASURE_YES) {	  
	  char label[200];
	  sprintf(label, "%s (boson), mass = %e:", 
		  force_label, bsn_mass[i]);
	  Fdt.print(dt, label);
	}

	//!< Do fermionic force contribution
	Fdt = lat.RHMC_EvolveMomFforce(mom_tmp, frmn+bsn_deg, frm_deg, 0,
				       frm_remez_arg_md[i].residue, frm_mass[i], 
				       dt, frmn_d, force_measure);
	if (force_measure == FORCE_MEASURE_YES) {	  
	  char label[200];
	  sprintf(label, "%s (fermion), mass = %e:", 
		  force_label, frm_mass[i]);
	  Fdt.print(dt, label);
	}

	//!< Monitor total force contribution
	if (force_measure == FORCE_MEASURE_YES) {
      Float L1 = 0.0;
      Float L2 = 0.0;
      Float Linf = 0.0;
	  for (int k=0; k<g_size/18; k++) {
	    Float norm = (mom_tmp+k)->norm();
	    Float tmp = sqrt(norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp>Linf ? tmp : Linf);
	  }
	  glb_sum(&L1);
	  glb_sum(&L2);
	  glb_max(&Linf);

	  L1 /= 4.0*GJP.VolSites();
	  L2 /= 4.0*GJP.VolSites();	 

	  fTimesV1PlusV2((IFloat*)mom,1.0,(IFloat*)mom_tmp,(IFloat*)mom,g_size);
	  sfree(mom_tmp);

	  char label[200];
	  sprintf(label, "%s (total), mass = (%e,%e):", 
		  force_label, frm_mass[i], bsn_mass[i]);

	  Fdt = ForceArg(L1, sqrt(L2), Linf);
	  Fdt.print(dt, label);

	}

      }
      
      evolved = 1;
      heatbathEval = 0;
      energyEval = 0;
      md_steps++;
    }
    
    LatticeFactory::Destroy();

  }

}

CPS_END_NAMESPACE
