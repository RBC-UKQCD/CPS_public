#include<config.h>
#include<string.h>
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
#include<util/time_cps.h>
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
    bsn_mass = (Float*) smalloc(n_masses * sizeof(Float), "bsn_mass", fname, cname);
    frm_mass = (Float*) smalloc(n_masses * sizeof(Float), "frm_mass", fname, cname);

    for(int i=0; i<n_masses; i++) {
      bsn_mass[i] = rat_quo_arg->bsn_mass.bsn_mass_val[i];
      frm_mass[i] = rat_quo_arg->frm_mass.frm_mass_val[i];
    }
    
    //!< construct approximation if necessary
    if(!loadPoles()) {
      generateApprox(frm_mass,&frm_remez_arg_md,&frm_remez_arg_mc,
                     rat_quo_arg->fermions.fermions_val);
      generateApprox(bsn_mass,&bsn_remez_arg_md,&bsn_remez_arg_mc,
                     rat_quo_arg->bosons.bosons_val);
      savePoles();
    }
    generateCgArg(frm_mass,&frm_cg_arg_fg, &frm_cg_arg_md,&frm_cg_arg_mc,"frm_cg_arg",
		  rat_quo_arg->fermions.fermions_val);
    generateCgArg(bsn_mass,&bsn_cg_arg_fg, &bsn_cg_arg_md,&bsn_cg_arg_mc,"bsn_cg_arg",
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
    destroyCgArg(bsn_cg_arg_fg, bsn_cg_arg_md, bsn_cg_arg_mc, "bsn_cg_arg", 
		 bsn_remez_arg_md, bsn_remez_arg_mc);
    destroyCgArg(frm_cg_arg_fg, frm_cg_arg_md, frm_cg_arg_mc, "frm_cg_arg", 
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

  Float dtime = -dclock();

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

  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
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
      Float dtime = -dclock();

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

    dtime += dclock();
    print_flops(cname, fname, 0, dtime);
    return h;
  }
}

void AlgActionRationalQuotient::prepare_fg(Matrix * force, Float dt_ratio)
{
  char * fname = "prepare_fg(M*,F)";
  Float dtime = -dclock();
  Float dtime_cg = 0.;
  Float dtime_force = 0.;

  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  

  for(int i=0; i<n_masses; i++){	
    int bsn_deg = bsn_remez_arg_md[i].degree;
    int frm_deg = frm_remez_arg_md[i].degree;

    //! First apply boson rational
    dtime_cg -= dclock();
    int shift = 0;
    cg_iter = lat.FmatEvlMInv(frmn, phi[i],
                              bsn_remez_arg_md[i].pole,
                              bsn_deg, 0,
                              bsn_cg_arg_fg[i], CNV_FRM_NO,
                              frmn_d+shift);
    dtime_cg += dclock();

    updateCgStats(bsn_cg_arg_fg[i][0]);

    //!< Construct rhs
    eta[0] -> VecEqualsVecTimesEquFloat(phi[i],bsn_remez_arg_md[i].norm,
                                        f_size);
    for (int j=0; j<bsn_deg; j++)
      eta[0] -> FTimesV1PlusV2(bsn_remez_arg_md[i].residue[j],
                               frmn[j],eta[0],f_size);

    //!< Now apply fermion rational
    dtime_cg -= dclock();
    shift += bsn_deg;
    cg_iter = lat.FmatEvlMInv(frmn+shift, eta[0], 
                              frm_remez_arg_md[i].pole, 
                              frm_deg, 0,
                              frm_cg_arg_fg[i], CNV_FRM_NO, 
                              frmn_d+shift);
    dtime_cg += dclock();

    updateCgStats(frm_cg_arg_fg[i][0]);
	
    //!< Construct rhs
    eta[1] -> VecEqualsVecTimesEquFloat(eta[0],frm_remez_arg_md[i].norm,
                                        f_size);
    for (int j=0; j<frm_deg; j++)
      eta[1] -> FTimesV1PlusV2(frm_remez_arg_md[i].residue[j],
                               frmn[j+shift],eta[1],f_size);

    //!< Apply final boson rational
    dtime_cg -= dclock();
    shift += frm_deg;
    cg_iter = lat.FmatEvlMInv(frmn+shift, eta[1], 
                              bsn_remez_arg_md[i].pole, 
                              bsn_deg, 0, 
                              bsn_cg_arg_fg[i], CNV_FRM_NO, 
                              frmn_d+shift);	
    dtime_cg += dclock();

    updateCgStats(bsn_cg_arg_fg[i][0]);

    //!< Now construct additional vectors needed
    for (int j=0; j<bsn_deg; j++) {
      Float one=1.;
      frmn[j+shift+bsn_deg] -> FTimesV1PlusV2(one, frmn[j], frmn[j+shift], f_size);
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

    Matrix * mom_tmp = force;
    if (force_measure == FORCE_MEASURE_YES) {
      mom_tmp = (Matrix*)smalloc(g_size*sizeof(Float),cname, fname, "mom_tmp");
      ((Vector*)mom_tmp)->VecZero(g_size);
    }
    
    dtime_force -= dclock();
    //!< Do bosonic force contribution
    Fdt = lat.RHMC_EvolveMomFforce(mom_tmp, frmn_tmp, 3*bsn_deg, 0,
                                   all_res, bsn_mass[i], dt_ratio, frmn_d, 
                                   force_measure);
    if (force_measure == FORCE_MEASURE_YES) {	  
      char label[200];
      sprintf(label, "%s (boson), mass = %e:", 
              force_label, bsn_mass[i]);
      Fdt.print(dt_ratio, label);
    }

    //!< Do fermionic force contribution
    Fdt = lat.RHMC_EvolveMomFforce(mom_tmp, frmn+bsn_deg, frm_deg, 0,
                                   frm_remez_arg_md[i].residue, frm_mass[i], 
                                   dt_ratio, frmn_d, force_measure);
    dtime_force += dclock();

    if (force_measure == FORCE_MEASURE_YES) {
      char label[200];
      sprintf(label, "%s (fermion), mass = %e:", 
              force_label, frm_mass[i]);
      Fdt.print(dt_ratio, label);
    }

    //!< Monitor total force contribution
    if (force_measure == FORCE_MEASURE_YES) {
        Fdt.measure(mom_tmp);
        Fdt.glb_reduce();

        ((Vector *)force)->VecAddEquVec((Vector *)mom_tmp, g_size);

        char label[200];
        sprintf(label, "%s (total), mass = (%e,%e):", 
                force_label, frm_mass[i], bsn_mass[i]);

        Fdt.print(dt_ratio, label);

        sfree(mom_tmp, "mom_tmp", fname, cname);
    }
  }
  LatticeFactory::Destroy();

  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
  print_flops(cname, "prepare_fg::cg()", 0, dtime_cg);
  print_flops(cname, "prepare_fg::force()", 0, dtime_force);
}

//!< run method evolves the integrator
void AlgActionRationalQuotient::evolve(Float dt, int nsteps)
{
  char * fname = "evolve(Float, int)";

  Float dtime = -dclock();
  Float dtime_cg = 0.;
  Float dtime_force = 0.;

  //!< Create an appropriate lattice
  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  
    
  for(int steps = 0; steps<nsteps; steps++) {
    for(int i=0; i<n_masses; i++){
      int bsn_deg = bsn_remez_arg_md[i].degree;
      int frm_deg = frm_remez_arg_md[i].degree;

      //! First apply boson rational
      int shift = 0;
      dtime_cg -= dclock();
      cg_iter = lat.FmatEvlMInv(frmn, phi[i], 
                                bsn_remez_arg_md[i].pole, 
                                bsn_deg, 0, 
                                bsn_cg_arg_md[i], CNV_FRM_NO, 
                                frmn_d+shift);	
      dtime_cg += dclock();

      updateCgStats(bsn_cg_arg_md[i][0]);

      //!< Construct rhs
      eta[0] -> VecEqualsVecTimesEquFloat(phi[i],bsn_remez_arg_md[i].norm,
                                          f_size);
      for (int j=0; j<bsn_deg; j++)
        eta[0] -> FTimesV1PlusV2(bsn_remez_arg_md[i].residue[j],
                                 frmn[j],eta[0],f_size);

      dtime_cg -= dclock();
      //!< Now apply fermion rational
      shift += bsn_deg;
      cg_iter = lat.FmatEvlMInv(frmn+shift, eta[0], 
                                frm_remez_arg_md[i].pole, 
                                frm_deg, 0, 
                                frm_cg_arg_md[i], CNV_FRM_NO, 
                                frmn_d+shift);
      dtime_cg += dclock();

      updateCgStats(frm_cg_arg_md[i][0]);
	
      //!< Construct rhs
      eta[1] -> VecEqualsVecTimesEquFloat(eta[0],frm_remez_arg_md[i].norm,
                                          f_size);
      for (int j=0; j<frm_deg; j++)
        eta[1] -> FTimesV1PlusV2(frm_remez_arg_md[i].residue[j],
                                 frmn[j+shift],eta[1],f_size);

      dtime_cg -= dclock();
      //!< Apply final boson rational
      shift += frm_deg;
      cg_iter = lat.FmatEvlMInv(frmn+shift, eta[1], 
                                bsn_remez_arg_md[i].pole, 
                                bsn_deg, 0, 
                                bsn_cg_arg_md[i], CNV_FRM_NO, 
                                frmn_d+shift);	
      dtime_cg += dclock();

      updateCgStats(bsn_cg_arg_md[i][0]);

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

      dtime_force -= dclock();
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
      dtime_force += dclock();

      if (force_measure == FORCE_MEASURE_YES) {	  
        char label[200];
        sprintf(label, "%s (fermion), mass = %e:", force_label, frm_mass[i]);
        Fdt.print(dt, label);
      }

      //!< Monitor total force contribution
      if (force_measure == FORCE_MEASURE_YES) {
          Fdt.measure(mom_tmp);
          Fdt.glb_reduce();

          fTimesV1PlusV2((IFloat*)mom,1.0,(IFloat*)mom_tmp,(IFloat*)mom,g_size);
          sfree(mom_tmp);

          char label[200];
          sprintf(label, "%s (total), mass = (%e,%e):", 
                  force_label, frm_mass[i], bsn_mass[i]);
          
          Fdt.print(dt, label);
      }
    }

    evolved = 1;
    heatbathEval = 0;
    energyEval = 0;
    md_steps++;
  }
    
  LatticeFactory::Destroy();

  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
  print_flops(cname, "evolve::cg()", 0, dtime_cg);
  print_flops(cname, "evolve::force()", 0, dtime_force);
}

bool AlgActionRationalQuotient::checkPolesFile(const RemezArg &md, const RemezArg &mc, const RationalDescr &r)
{
  if(md.field_type != r.field_type) return false;
  if(mc.field_type != r.field_type) return false;

  if(mc.power_num != r.power_num) return false;
  if(mc.power_den != r.power_den * 2) return false;

  if(md.power_num != r.power_num) return false;
  if(r.field_type == BOSON) {
    if(md.power_den != r.power_den * 2) return false;
  } else {
    if(md.power_den != r.power_den) return false;
  }

  if(md.degree != r.md_approx.stop_rsd.stop_rsd_len) return false;
  if(mc.degree != r.mc_approx.stop_rsd.stop_rsd_len) return false;

  if(fabs(md.lambda_low  - r.md_approx.lambda_low ) > 1e-3 * fabs(r.md_approx.lambda_low ) ) return false;
  if(fabs(md.lambda_high - r.md_approx.lambda_high) > 1e-3 * fabs(r.md_approx.lambda_high) ) return false;
  if(fabs(mc.lambda_low  - r.mc_approx.lambda_low ) > 1e-3 * fabs(r.mc_approx.lambda_low ) ) return false;
  if(fabs(mc.lambda_high - r.mc_approx.lambda_high) > 1e-3 * fabs(r.mc_approx.lambda_high) ) return false;

  return true;
}

// return true if we successfully loaded from a file.
bool AlgActionRationalQuotient::loadPoles(void)
{
  const char *fname = "loadPoles()";
  if(rat_quo_arg->remez_generate) return false;
  if(strlen(rat_quo_arg->rat_poles_file) == 0) return false;

  FILE *fp = fopen(rat_quo_arg->rat_poles_file, "r");
  if(fp == NULL) return false;
  fclose(fp);

  RationalQuotientRemezArg rq;
  if(!rq.Decode(rat_quo_arg->rat_poles_file, "rq")) return false;

  // a bunch of check
  if(rq.bsn_md.bsn_md_len != n_masses) return false;
  if(rq.bsn_mc.bsn_mc_len != n_masses) return false;
  if(rq.frm_md.frm_md_len != n_masses) return false;
  if(rq.frm_mc.frm_mc_len != n_masses) return false;

  frm_remez_arg_md = new RemezArg[n_masses];
  frm_remez_arg_mc = new RemezArg[n_masses];
  bsn_remez_arg_md = new RemezArg[n_masses];
  bsn_remez_arg_mc = new RemezArg[n_masses];
  
  // we don't try to make a deep copy since VML does not free space
  // used anyway.
  for(int i = 0; i < n_masses; ++i) {
    frm_remez_arg_md[i] = rq.frm_md.frm_md_val[i];
    frm_remez_arg_mc[i] = rq.frm_mc.frm_mc_val[i];
    bsn_remez_arg_md[i] = rq.bsn_md.bsn_md_val[i];
    bsn_remez_arg_mc[i] = rq.bsn_mc.bsn_mc_val[i];

    if(! checkPolesFile(frm_remez_arg_md[i], frm_remez_arg_mc[i],
                        rat_quo_arg->fermions.fermions_val[i])) return false;
    if(! checkPolesFile(bsn_remez_arg_md[i], bsn_remez_arg_mc[i],
                        rat_quo_arg->bosons.bosons_val[i])) return false;
  }

  VRB.Result(cname, fname, "Successfully loaded poles file %s.\n", rat_quo_arg->rat_poles_file);
  return true;
}

bool AlgActionRationalQuotient::savePoles(void)
{
  if(strlen(rat_quo_arg->rat_poles_file) == 0) return false;

  RationalQuotientRemezArg rq;

  rq.bsn_md.bsn_md_len = n_masses;
  rq.bsn_mc.bsn_mc_len = n_masses;
  rq.frm_md.frm_md_len = n_masses;
  rq.frm_mc.frm_mc_len = n_masses;

  // no deep copy either
  rq.bsn_md.bsn_md_val = bsn_remez_arg_md;
  rq.bsn_mc.bsn_mc_val = bsn_remez_arg_mc;
  rq.frm_md.frm_md_val = frm_remez_arg_md;
  rq.frm_mc.frm_mc_val = frm_remez_arg_mc;

  return rq.Encode(rat_quo_arg->rat_poles_file, "rq");
}

CPS_END_NAMESPACE
