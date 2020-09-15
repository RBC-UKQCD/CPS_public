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
#include<util/lattice/fbfm.h>
#include<alg/alg_int.h>
#include<alg/alg_remez.h>
#include <util/timer.h>
#include <util/lattice/fbfm.h>
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
//    skip_force(false)
{

  cname = "AlgActionRationalQuotient";
  const char *fname = "AlgActionRationalQuotient()";

  if(!UniqueID()){ printf("AlgActionRationalQuotient constructor started\n"); fflush(stdout); }

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
#ifdef USE_BFM
      if (rat_quo_arg->bi_arg.fermion == F_CLASS_BFM) {
	  // AlgActionBilinear does not set fermion field size correctly for Fbfm
	  int Ls = Fbfm::arg_map.at(rat_quo_arg->bsn_mass.bsn_mass_val[0]).Ls;

	  VRB.Result(cname, fname, "Recalculating fermion field size for Fbfm based on Ls = %d\n", Ls);

	  //!< Number of Floats in a Vector array
	  f_size = GJP.VolNodeSites() * Ls * (2 * 3 * 4) / 2; // (reim * color * spin) / Ncheckerboard
	  //!< Number of Vectors in a Vector array
	  f_vec_count = f_size / (2 * 3);
	  //!< Number of lattice sites
	  f_sites = f_size / (2 * 3 * 4);

	  VRB.Result(cname, fname, "Allocating phi fields\n");
	  for (int i = 0; i < n_masses; i++) {
	      phi[i] = (Vector *)smalloc(f_size*sizeof(Float), "phi[i]", fname, cname);
	  }
      }
#endif


    bsn_mass = (Float*) smalloc(n_masses * sizeof(Float), "bsn_mass", fname, cname);
    frm_mass = (Float*) smalloc(n_masses * sizeof(Float), "frm_mass", fname, cname);

    for(int i=0; i<n_masses; i++) {
      bsn_mass[i] = rat_quo_arg->bsn_mass.bsn_mass_val[i];
      frm_mass[i] = rat_quo_arg->frm_mass.frm_mass_val[i];

#ifdef USE_BFM
      if (rat_quo_arg->bi_arg.fermion == F_CLASS_BFM) {
	  // Make sure all quotients have the same Ls
	  int Ls = Fbfm::arg_map.at(bsn_mass[0]).Ls;
	  if (Fbfm::arg_map.at(bsn_mass[i]).Ls != Ls) {
	      ERR.General(cname, fname, "Boson mass #%d doesn't have the same Ls as boson mass #0!\n", i);
	  }
	  if (Fbfm::arg_map.at(frm_mass[i]).Ls != Ls) {
	      ERR.General(cname, fname, "Fermion mass #%d doesn't have the same Ls as boson mass #0!\n", i);
	  }
      }
#endif
    }

    //CK: added for twisted mass fermions
    bsn_mass_epsilon = (Float*) smalloc(n_masses * sizeof(Float), "bsn_mass_epsilon", fname, cname);
    frm_mass_epsilon = (Float*) smalloc(n_masses * sizeof(Float), "frm_mass_epsilon", fname, cname);
    
    for(int i=0; i<n_masses; i++) {
      bsn_mass_epsilon[i] = rat_quo_arg->bsn_mass_epsilon.bsn_mass_epsilon_val[i];
      frm_mass_epsilon[i] = rat_quo_arg->frm_mass_epsilon.frm_mass_epsilon_val[i];
    }
    
    //!< construct approximation if necessary
    if(!loadPoles()) {
      //CK: Note that the Remez approximation does not actually depend upon the quark mass. The mass
      //    is only used to bound the eigenvalues when RATIONAL_BOUNDS_AUTOMATIC is chosen (currently
      //    this is only implemented for ASQTAD, Staggered and P4 actions)
      //    This means that we do not have to modify the Remez part to accomodate the twist term for 
      //    twisted Wilson fermions
      if(!UniqueID()){ printf("Generating rational approximation\n"); fflush(stdout); }

      generateApprox(frm_mass,&frm_remez_arg_md,&frm_remez_arg_mc,
                     rat_quo_arg->fermions.fermions_val);
      generateApprox(bsn_mass,&bsn_remez_arg_md,&bsn_remez_arg_mc,
                     rat_quo_arg->bosons.bosons_val);
      savePoles();
      if(!UniqueID()){ printf("Finished generating rational approximation\n"); fflush(stdout); }
    }
    generateCgArg(frm_mass,frm_mass_epsilon,&frm_cg_arg_fg, &frm_cg_arg_md,&frm_cg_arg_mc,"frm_cg_arg",
		  rat_quo_arg->fermions.fermions_val);
    generateCgArg(bsn_mass,bsn_mass_epsilon,&bsn_cg_arg_fg, &bsn_cg_arg_md,&bsn_cg_arg_mc,"bsn_cg_arg",
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
    
    VRB.Result(cname, fname, "allocating fermion fields of size %d Floats\n", f_size);

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

  if(!UniqueID()){ printf("AlgActionRationalQuotient constructor finished\n"); fflush(stdout); }
}

AlgActionRationalQuotient::~AlgActionRationalQuotient() {

  const char *fname = "~AlgActionRationalQuotient()" ;
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

    //CK: Added for twisted mass fermions
    sfree(bsn_mass_epsilon,"bsn_mass_epsilon",fname,cname);
    sfree(frm_mass_epsilon,"frm_mass_epsilon",fname,cname);

    if (rat_quo_arg->eigen.eigen_measure == EIGEN_MEASURE_YES) destroyEigArg();
  }

}

//!< Heat Bath for the pseudo-fermions (phi)
void AlgActionRationalQuotient::reweight(Float *rw_fac,Float *norm) {

  const char *fname = "reweight()";


//    if (n_masses>1)
//      ERR.General(cname,fname,"not implemented for n_mass>1(%d)\n",n_masses);

    //!< Before energy is measured, do we want to check bounds?
    if (rat_quo_arg->eigen.eigen_measure == EIGEN_MEASURE_YES) {
      checkApprox(bsn_mass, bsn_mass_epsilon, bsn_remez_arg_mc, rat_quo_arg->eigen);
      checkApprox(frm_mass, frm_mass_epsilon, frm_remez_arg_mc, rat_quo_arg->eigen);
    }

    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  
    
    for(int i=0; i<n_masses; i++){
#ifdef USE_BFM
	if (rat_quo_arg->bi_arg.fermion == F_CLASS_BFM) {
	    // Fbfm needs current_key_mass set before calling RandGaussVector
	    Fbfm::current_key_mass = bsn_mass[i];
	    VRB.Result(cname, fname, "Setting Fbfm::current_key_mass = %e before calling RandGaussVector\n", Fbfm::current_key_mass);
	}
#endif
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
  char fname[20+strlen(force_label)];
  sprintf(fname, "heatbath()[%s]",force_label);
    
  //const char *fname = "heatbath()";
  static Timer timer(cname, fname);
  timer.start(true);
  Float dtime = -dclock();

  //!< Only evaluate heatbath if necessary
  if (!heatbathEval) {

    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  
    h_init = 0.0;
    
    for(int i=0; i<n_masses; i++){
#ifdef USE_BFM
	if (rat_quo_arg->bi_arg.fermion == F_CLASS_BFM) {
	    // Fbfm needs current_key_mass set before calling RandGaussVector
	    Fbfm::current_key_mass = bsn_mass[i];
	    VRB.Result(cname, fname, "Setting Fbfm::current_key_mass = %e before calling RandGaussVector\n", Fbfm::current_key_mass);
	}
#endif
      lat.RandGaussVector(phi[i], 0.5, Ncb);
      VRB.Result(cname,fname,
	"phi vector for mass %d: %.9e %.9e %.9e .....\n",force_label,i, ((Float*)phi[i])[0],((Float*)phi[i])[1], ((Float*)phi[i])[2]);

      Float h_i = lat.FhamiltonNode(phi[i],phi[i]);
      h_init += h_i;
      VRB.Result(cname,fname,"h_init=%0.16e\n",h_init);

      Float total_h_i = h_i;
      glb_sum(&total_h_i);
      VRB.Result(cname, fname, "heatbath: mass ratio %0.4f/%0.4f initial ham = %0.16e\n", frm_cg_arg_mc[i][0]->mass, bsn_cg_arg_mc[i][0]->mass, total_h_i);

      if(GJP.Gparity1fX() && GJP.Gparity1fY()){
      	if(!UniqueID()){ printf("Putting minus sign on fermion source in UR quadrant\n"); fflush(stdout); }
	//make source on upper-right quadrant negative (RNGs should be correct)
      	for(int s=0;s<GJP.SnodeSites();s++){
      	  for(int t=0;t<GJP.TnodeSites();t++){
      	    for(int z=0;z<GJP.ZnodeSites();z++){
      	      for(int y=0;y<GJP.YnodeSites();y++){
      		for(int x=0;x<GJP.XnodeSites();x++){
      		  if( (x+y+z+t+s)%2 == 0) continue; //ferm vect is odd parity only

      		  int gx = x+GJP.XnodeCoor()*GJP.XnodeSites();
      		  int gy = y+GJP.YnodeCoor()*GJP.YnodeSites();

      		  if(gx>=GJP.Xnodes()*GJP.XnodeSites()/2 && gy>=GJP.Ynodes()*GJP.YnodeSites()/2){
      		    int pos[5] = {x,y,z,t,s};
      		    int f_off = lat.FsiteOffsetChkb(pos) * lat.SpinComponents();

      		    for(int spn=0;spn<lat.SpinComponents();spn++) *(frmn[0]+f_off+spn) *=-1;
      		  }
      		}
      	      }
      	    }
      	  }
      	}
      }

      VRB.Result(cname,fname, "phi vector for mass %d: %.9e %.9e %.9e .....\n",force_label,i, ((Float*)phi[i])[0],((Float*)phi[i])[1], ((Float*)phi[i])[2]);
      
      Float delta_h = lat.FhamiltonNode(phi[i],phi[i]);
      {
	Float gsum_h(delta_h);
	glb_sum(&gsum_h);
      VRB.Result(cname,fname, "delta H for mass %d:  %.9e\n",i,gsum_h);
      }
     // h_init += delta_h;

      //!< First apply the fermion rational
      frmn[0] -> 
	VecEqualsVecTimesEquFloat(phi[i],frm_remez_arg_mc[i].norm_inv,f_size);
      
      VRB.Result(cname, fname, "Calling lat.FmatEvlMInv for fermion mass %d\n", i);
      cg_iter = lat.FmatEvlMInv(frmn, phi[i], frm_remez_arg_mc[i].pole_inv, 
				frm_remez_arg_mc[i].degree, 0, 
				frm_cg_arg_mc[i], CNV_FRM_NO, SINGLE,
				frm_remez_arg_mc[i].residue_inv);
     {
        Float * tmp_f = (Float *)frmn[0];
        VRB.Result(cname,fname,"frm: phi=%g\n",*tmp_f);
      }

      VRB.Result(cname, fname, "fermion mass %e multishift inversion cg_iter = %d\n", frm_cg_arg_mc[i][0]->mass, cg_iter);

      //!< Now apply the boson rational
      phi[i] -> 
	VecEqualsVecTimesEquFloat(frmn[0],bsn_remez_arg_mc[i].norm_inv,f_size);

      VRB.Result(cname, fname, "Calling lat.FmatEvlMInv for boson mass %d\n", i);
      cg_iter = lat.FmatEvlMInv(phi+i, frmn[0], bsn_remez_arg_mc[i].pole_inv, 
				bsn_remez_arg_mc[i].degree, 0, 
				bsn_cg_arg_mc[i], CNV_FRM_NO, SINGLE,
				bsn_remez_arg_mc[i].residue_inv);
      VRB.Result(cname, fname, "boson mass %e multishift inversion cg_iter = %d\n", bsn_cg_arg_mc[i][0]->mass, cg_iter);
      {
        Float * tmp_f = (Float *)phi[i];
        VRB.Result(cname,fname,"bsn: phi=%g\n",*tmp_f);
      }


      updateCgStats(frm_cg_arg_mc[i][0]);
      updateCgStats(bsn_cg_arg_mc[i][0]);
    }

    LatticeFactory::Destroy();

    evolved = 0;
    heatbathEval = 1;
    energyEval = 0;
//    traj++;
  }
  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
  timer.stop(true);
}

// Calculate rhmc fermion contribution to the Hamiltonian
Float AlgActionRationalQuotient::energy() {

  //const char *fname="energy()";
  char fname[20+strlen(force_label)];
  sprintf(fname, "energy()[%s]",force_label);
   VRB.Result(cname,fname,"energeEval evolved h_init=%d %d %g\n",
        energyEval,evolved,h_init);


  if (energyEval) {
    return 0.0;
  } else if (!evolved) {
    energyEval = 1;
    {
    Float glb_h = h_init;
    glb_sum(&glb_h);
//    VRB.Result(cname, fname, "glb_h = %0.15e\n", glb_h);
      if(UniqueID()==0)   printf("AlgActionRationalQuotient::energy() [%s] %.9e\n",force_label,glb_h);
    }
    return h_init;
  } else {
    static Timer timer(cname, fname);
    timer.start(true);
    Float dtime = -dclock();
    int shift = 0;
    Float h = 0.0;

    //!< Before energy is measured, do we want to check bounds?
    if (rat_quo_arg->eigen.eigen_measure == EIGEN_MEASURE_YES) {
      VRB.Result(cname, fname, "Checking boson eigenvalue bounds\n");
      if(!UniqueID()) printf("AlgActionRationalQuotient::energy() [%s] checking bosonic eigenvalue bounds\n",force_label);
      checkApprox(bsn_mass, bsn_mass_epsilon, bsn_remez_arg_mc, rat_quo_arg->eigen);
      VRB.Result(cname, fname, "Checking fermion eigenvalue bounds\n");
      if(!UniqueID()) printf("AlgActionRationalQuotient::energy() [%s] checking fermionic eigenvalue bounds\n",force_label);
      checkApprox(frm_mass, frm_mass_epsilon, frm_remez_arg_mc, rat_quo_arg->eigen);
    }

    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  

    Float total_h_i;
    for (int i=0; i<n_masses; i++) {

      //!< First apply boson rational
      frmn[0] -> VecEqualsVecTimesEquFloat(phi[i],bsn_remez_arg_mc[i].norm,
					   f_size);
      
      VRB.Result(cname, fname, "Calling lat.FmatEvlMInv for boson mass %d\n", i);
      cg_iter = lat.FmatEvlMInv(frmn, phi[i], bsn_remez_arg_mc[i].pole, 
				bsn_remez_arg_mc[i].degree, 0, 
				bsn_cg_arg_mc[i], CNV_FRM_NO, SINGLE, 
				bsn_remez_arg_mc[i].residue);
      VRB.Result(cname, fname, "boson mass %e multishift inversion cg_iter = %d\n", bsn_cg_arg_mc[i][0]->mass, cg_iter);
      
      //!< Now apply fermion rational
      frmn[1] -> VecEqualsVecTimesEquFloat(frmn[0],frm_remez_arg_mc[i].norm,
					   f_size);
      
      VRB.Result(cname, fname, "Calling lat.FmatEvlMInv for fermion mass %d\n", i);
      cg_iter = lat.FmatEvlMInv(frmn+1, frmn[0], frm_remez_arg_mc[i].pole, 
				frm_remez_arg_mc[i].degree, 0, 
				frm_cg_arg_mc[i], CNV_FRM_NO, SINGLE, 
				frm_remez_arg_mc[i].residue);
      VRB.Result(cname, fname, "fermion mass %e multishift inversion cg_iter = %d\n", frm_cg_arg_mc[i][0]->mass, cg_iter);

      updateCgStats(bsn_cg_arg_mc[i][0]);
      updateCgStats(frm_cg_arg_mc[i][0]);

      // shift this evaluation into minvcg?
      Float h_i = lat.FhamiltonNode(frmn[1], frmn[1]);
      h += h_i;
      total_h_i = h_i;
      glb_sum(&total_h_i);
      VRB.Result(cname, fname, "energy: mass ratio %0.4f/%0.4f final ham = %0.16e\n", frm_cg_arg_mc[i][0]->mass, bsn_cg_arg_mc[i][0]->mass, total_h_i);
    }

    LatticeFactory::Destroy();

    energyEval = 1;

    dtime += dclock();
    print_flops(cname, fname, 0, dtime);
    timer.stop(true);
    {
      Float gsum_h(h);
      glb_sum(&gsum_h);
      if(UniqueID()==0)   printf("AlgActionRationalQuotient::energy() [%s] %.16e\n",force_label,gsum_h);
    }
    total_h_i = h-h_init;
    glb_sum(&total_h_i);
    VRB.Result(cname, fname, "energy: delta_h = %0.16e\n", total_h_i);

    return h;
  }
}

void AlgActionRationalQuotient::prepare_fg(Matrix * force, Float dt_ratio)
{
  //const char * fname = "prepare_fg(M*,F)";

  char fname[30+strlen(force_label)];
  sprintf(fname, "prepare_fg(M*,F)[%s]",force_label);
  static Timer timer(cname, fname);
  timer.start(true);

  Float dtime = -dclock();
  Float dtime_cg = 0.;
  Float dtime_force = 0.;
  
  if (skip_force) {
    VRB.Result(cname, fname, "WARNING! skipping prepare_fg() because AlgActionRationalQuotient::skip_force is true!\n");
    evolved = 1;
    timer.stop(true);
    return;
  }

  if(!UniqueID()){    
    Float pvals[4];
    for(int ii=0;ii<4;ii++){
      int off = 18 * ii + 2;
      pvals[ii] = ((Float*)force)[off];
    }
    VRB.Debug(cname,fname,"start, input temp conj mom Px(0) = %.9e, Py(0) = %.9e, Pz(0) = %.9e, Pt(0) = %.9e\n",force_label,pvals[0],pvals[1],pvals[2],pvals[3]);
  }  

  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  

  for(int i=0; i<n_masses; i++){	
    int bsn_deg = bsn_remez_arg_md[i].degree;
    int frm_deg = frm_remez_arg_md[i].degree;

    //! First apply boson rational
    dtime_cg -= dclock();
    int shift = 0;
    VRB.Result(cname, fname, "Calling first lat.FmatEvlMInv for boson mass %d\n", i);
    cg_iter = lat.FmatEvlMInv(frmn, phi[i],
                              bsn_remez_arg_md[i].pole,
                              bsn_deg, 0,
                              bsn_cg_arg_fg[i], CNV_FRM_NO,
                              frmn_d+shift);
    dtime_cg += dclock();
    VRB.Result(cname, fname, "boson mass %e first multishift inversion cg_iter = %d\n", bsn_cg_arg_fg[i][0]->mass, cg_iter);

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
    VRB.Result(cname, fname, "Calling lat.FmatEvlMInv for fermion mass %d\n", i);
    cg_iter = lat.FmatEvlMInv(frmn+shift, eta[0], 
                              frm_remez_arg_md[i].pole, 
                              frm_deg, 0,
                              frm_cg_arg_fg[i], CNV_FRM_NO, 
                              frmn_d+shift);
    dtime_cg += dclock();
    VRB.Result(cname, fname, "fermion mass %e multishift inversion cg_iter = %d\n", frm_cg_arg_fg[i][0]->mass, cg_iter);

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
    VRB.Result(cname, fname, "Calling second lat.FmatEvlMInv for boson mass %d\n", i);
    cg_iter = lat.FmatEvlMInv(frmn+shift, eta[1], 
                              bsn_remez_arg_md[i].pole, 
                              bsn_deg, 0, 
                              bsn_cg_arg_fg[i], CNV_FRM_NO, 
                              frmn_d+shift);	
    dtime_cg += dclock();
    VRB.Result(cname, fname, "boson mass %e second multishift inversion cg_iter = %d\n", bsn_cg_arg_fg[i][0]->mass, cg_iter);

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
    //CK: Modified for twisted mass fermions
    if(fermion == F_CLASS_WILSON_TM) 
      Fdt = dynamic_cast<FwilsonTm&>(lat).RHMC_EvolveMomFforce(mom_tmp, frmn_tmp, 3*bsn_deg, 0,
							       all_res, bsn_mass[i], bsn_mass_epsilon[i], dt_ratio, frmn_d, 
							       force_measure);
    else
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
    //CK: Modified for twisted mass fermions
    if(fermion == F_CLASS_WILSON_TM) 
      Fdt = dynamic_cast<FwilsonTm&>(lat).RHMC_EvolveMomFforce(mom_tmp, frmn+bsn_deg, frm_deg, 0,
							       frm_remez_arg_md[i].residue, frm_mass[i], frm_mass_epsilon[i],
							       dt_ratio, frmn_d, force_measure);
    else
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

  char fname_cg[30+strlen(force_label)];
  sprintf(fname_cg,"prepare_fg::cg()[%s]",force_label);
  char fname_force[30+strlen(force_label)];
  sprintf(fname_force,"prepare_fg::force()[%s]",force_label);

  print_flops(cname, fname, 0, dtime);
  print_flops(cname, fname_cg, 0, dtime_cg);
  print_flops(cname, fname_force, 0, dtime_force);
  timer.stop(true);
  //print_flops(cname, "prepare_fg::cg() [%s]", 0, dtime_cg);
  //print_flops(cname, "prepare_fg::force() [%s]", 0, dtime_force);

  if(!UniqueID()){    
    Float pvals[4];
    for(int ii=0;ii<4;ii++){
      int off = 18 * ii + 2;
      pvals[ii] = ((Float*)force)[off];
    }
    if(UniqueID()==0) printf("AlgActionRationalQuotient::prepare_fg() [%s] end, output temp conj mom Px(0) = %.9e, Py(0) = %.9e, Pz(0) = %.9e, Pt(0) = %.9e\n",force_label,pvals[0],pvals[1],pvals[2],pvals[3]);
  } 

}

//!< run method evolves the integrator
void AlgActionRationalQuotient::evolve(Float dt, int nsteps)
{
  //const char * fname = "evolve(Float, int)";
  char fname[30+strlen(force_label)];
  sprintf(fname, "evolve(F,i)[%s]",force_label);
  static Timer timer(cname, fname);
  timer.start(true);

  Float dtime = -dclock();
  Float dtime_cg = 0.;
  Float dtime_force = 0.;

  if (skip_force) {
    VRB.Result(cname, fname, "WARNING! skipping evolve() because AlgActionRationalQuotient::skip_force is true!\n");
    evolved = 1;
    timer.stop(true);
    return;
  }

  //!< Create an appropriate lattice
  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  
    
  {
    Float pvals[4];
    for(int ii=0;ii<4;ii++){
      int off = 18 * ii + 2;
      pvals[ii] = ((Float*)mom)[off];
    }
    if(UniqueID()==0) printf("AlgActionRationalQuotient::evolve() [%s] start, conj mom Px(0) = %.9e, Py(0) = %.9e, Pz(0) = %.9e, Pt(0) = %.9e\n",force_label,pvals[0],pvals[1],pvals[2],pvals[3]);
  }  


  for(int steps = 0; steps<nsteps; steps++) {
    for(int i=0; i<n_masses; i++){
      if(UniqueID()==0) printf("AlgActionRationalQuotient::evolve() [%s] step %d mass %d\n",force_label,steps,i);

      int bsn_deg = bsn_remez_arg_md[i].degree;
      int frm_deg = frm_remez_arg_md[i].degree;

      //! First apply boson rational
      int shift = 0;
      dtime_cg -= dclock();
//      VRB.Result(cname, fname, "Calling first lat.FmatEvlMInv for boson mass %d\n", i);
      cg_iter = lat.FmatEvlMInv(frmn, phi[i], 
                                bsn_remez_arg_md[i].pole, 
                                bsn_deg, 0, 
                                bsn_cg_arg_md[i], CNV_FRM_NO, 
                                frmn_d+shift);	
      dtime_cg += dclock();
      VRB.Result(cname, fname, "boson mass %e first multishift inversion cg_iter = %d\n", bsn_cg_arg_md[i][0]->mass, cg_iter);

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
      VRB.Result(cname, fname, "Calling lat.FmatEvlMInv for fermion mass %d\n", i);
      cg_iter = lat.FmatEvlMInv(frmn+shift, eta[0], 
                                frm_remez_arg_md[i].pole, 
                                frm_deg, 0, 
                                frm_cg_arg_md[i], CNV_FRM_NO, 
                                frmn_d+shift);
      dtime_cg += dclock();
      VRB.Result(cname, fname, "fermion mass %e multishift inversion cg_iter = %d\n", frm_cg_arg_md[i][0]->mass, cg_iter);

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
      VRB.Result(cname, fname, "Calling second lat.FmatEvlMInv for boson mass %d\n", i);
      cg_iter = lat.FmatEvlMInv(frmn+shift, eta[1], 
                                bsn_remez_arg_md[i].pole, 
                                bsn_deg, 0, 
                                bsn_cg_arg_md[i], CNV_FRM_NO, 
                                frmn_d+shift);	
      dtime_cg += dclock();
      VRB.Result(cname, fname, "boson mass %e second multishift inversion cg_iter = %d\n", bsn_cg_arg_md[i][0]->mass, cg_iter);

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
      //CK: Modified for twisted mass fermions
      if(fermion == F_CLASS_WILSON_TM) 
	Fdt = dynamic_cast<FwilsonTm&>(lat).RHMC_EvolveMomFforce(mom_tmp, frmn_tmp, 3*bsn_deg, 0,
								 all_res, bsn_mass[i], bsn_mass_epsilon[i], dt, frmn_d, 
								 force_measure);
      else
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
      if(fermion == F_CLASS_WILSON_TM) 
	Fdt = dynamic_cast<FwilsonTm&>(lat).RHMC_EvolveMomFforce(mom_tmp, frmn+bsn_deg, frm_deg, 0,
								 frm_remez_arg_md[i].residue, frm_mass[i], frm_mass_epsilon[i],
								 dt, frmn_d, force_measure);
      else
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
      {
	Float pvals[4];
	for(int ii=0;ii<4;ii++){
	  int off = 18 * ii + 2;
	  pvals[ii] = ((Float*)mom)[off];
	}
	if(UniqueID()==0) printf("AlgActionRationalQuotient::evolve() [%s] end of step %d, mass %d: conj mom Px(0) = %.9e, Py(0) = %.9e, Pz(0) = %.9e, Pt(0) = %.9e\n",
				 force_label,steps,i,pvals[0],pvals[1],pvals[2],pvals[3]);
      }  


    }

    evolved = 1;
    heatbathEval = 0;
    energyEval = 0;
    md_steps++;
  }

  {
    Float pvals[4];
    for(int ii=0;ii<4;ii++){
      int off = 18 * ii + 2;
      pvals[ii] = ((Float*)mom)[off];
    }
    if(UniqueID()==0) printf("AlgActionRationalQuotient::evolve() [%s] end, conj mom Px(0) = %e, Py(0) = %e, Pz(0) = %e, Pt(0) = %e\n",force_label,pvals[0],pvals[1],pvals[2],pvals[3]);
  }  
    
  LatticeFactory::Destroy();

  char fname_cg[30+strlen(force_label)];
  sprintf(fname_cg,"evolve::cg()[%s]",force_label);
  char fname_force[30+strlen(force_label)];
  sprintf(fname_force,"evolve::force()[%s]",force_label);

  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
  print_flops(cname, fname_cg, 0, dtime_cg);
  print_flops(cname, fname_force, 0, dtime_force);
  timer.stop(true);
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
  
  if(md.field_type != r.field_type) { printf("a\n"); return false; }
  if(mc.field_type != r.field_type) { printf("b\n"); return false; }

  if(mc.power_num != r.power_num) { printf("c\n"); return false; }
  if(mc.power_den != r.power_den * 2) { printf("d\n"); return false; }

  if(md.power_num != r.power_num) { printf("e\n"); return false; }
  if(r.field_type == BOSON) {
    if(md.power_den != r.power_den * 2) { printf("f\n"); return false; }
  } else {
    if(md.power_den != r.power_den) { printf("g\n"); return false; }
  }

  if(md.degree != r.md_approx.stop_rsd.stop_rsd_len) { printf("h\n"); return false; }
  if(mc.degree != r.mc_approx.stop_rsd.stop_rsd_len) { printf("i\n"); return false; }

  if(fabs(md.lambda_low  - r.md_approx.lambda_low ) > 1e-3 * fabs(r.md_approx.lambda_low ) ) { printf("j\n"); return false; }
  if(fabs(md.lambda_high - r.md_approx.lambda_high) > 1e-3 * fabs(r.md_approx.lambda_high) ) { printf("k\n"); return false; }
  if(fabs(mc.lambda_low  - r.mc_approx.lambda_low ) > 1e-3 * fabs(r.mc_approx.lambda_low ) ) { printf("l\n"); return false; }
  if(fabs(mc.lambda_high - r.mc_approx.lambda_high) > 1e-3 * fabs(r.mc_approx.lambda_high) ) { printf("m\n"); return false; }
  if(fabs(mc.lambda_high - r.mc_approx.lambda_high) > 1e-3 * fabs(r.mc_approx.lambda_high) ) return false;

  return true;
}

void PrintPFE(const char* name, Float norm, const Float* residue, const Float* pole, int degree)
{
  const char* cname = "AlgActionRationalQuotient";
  const char* fname = "PrintPFE";

  char pfe[10000];
  sprintf(pfe, "%.15e", norm);

  for(int i = 0; i < degree; i++) {
    char term[1024];
    sprintf(term, " + %.15e/(x + %.15e)", residue[i], pole[i]);
    strcat(pfe, term);
  }

  if(!UniqueID()) printf("PFE of %s:    %s\n", name, pfe);
}

void PrintRemezArg(const char* name, const RemezArg& remez_arg)
{
  char pfe_name[1024];

  sprintf(pfe_name, "%s normal", name);
  PrintPFE(pfe_name, remez_arg.norm, remez_arg.residue, remez_arg.pole, remez_arg.degree);
  sprintf(pfe_name, "%s inv", name);
  PrintPFE(pfe_name, remez_arg.norm_inv, remez_arg.residue_inv, remez_arg.pole_inv, remez_arg.degree);
}

// return true if we successfully loaded from a file.
bool AlgActionRationalQuotient::loadPoles(void)
{
  const char* fname = "loadPoles()";

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
  
    PrintRemezArg("frm_remez_arg_md", frm_remez_arg_md[i]);
    PrintRemezArg("frm_remez_arg_mc", frm_remez_arg_mc[i]);
    PrintRemezArg("bsn_remez_arg_md", bsn_remez_arg_md[i]);
    PrintRemezArg("bsn_remez_arg_mc", bsn_remez_arg_mc[i]);
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
