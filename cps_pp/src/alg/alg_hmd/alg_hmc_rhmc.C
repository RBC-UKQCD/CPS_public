#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
/*!\file
  \brief Definitions of the AlgHmcRHMC methods.

  $Id: alg_hmc_rhmc.C,v 1.28 2007-06-25 15:49:20 chulwoo Exp $
*/
//--------------------------------------------------------------------
/*
  $Author: chulwoo $
  $Date: 2007-06-25 15:49:20 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_hmd/alg_hmc_rhmc.C,v 1.28 2007-06-25 15:49:20 chulwoo Exp $
  $Id: alg_hmc_rhmc.C,v 1.28 2007-06-25 15:49:20 chulwoo Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: alg_hmc_rhmc.C,v $
  $Revision: 1.28 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_hmd/alg_hmc_rhmc.C,v $
  $State: Exp $
*/
//--------------------------------------------------------------------


//------------------------------------------------------------------
//
// alg_hmc_rhmc.C
//
// AlgHmcRHMC is derived from Alg and is relevant to the Rational
// Hybrid Monte Carlo Algorithm.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/data_shift.h>
#include<util/checksum.h>
#include<util/qcdio.h>
#include<math.h>
#include<alg/alg_hmd.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<comms/glb.h>
#include<alg/alg_remez.h>
#include<util/data_shift.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
CPS_START_NAMESPACE

//------------------------------------------------------------------
/*!
  \param latt The lattice on which the HMC algorithm runs.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
*/
//------------------------------------------------------------------

AlgHmcRHMC::AlgHmcRHMC(Lattice& latt, 
		       CommonArg *c_arg,
		       HmdArg *arg) : 
  AlgHmd(latt, c_arg, arg) 
{
  cname = "AlgHmcRHMC";
  char *fname = "AlgHmcRHMC(L&,CommonArg*,HmdArg*)";
  VRB.Func(cname,fname);

  // Approx type must be constant if this constructor is called
  //----------------------------------------------------------------
  hmd_arg->approx_type = CONSTANT;

  // Currently RHMC force term is not implemented for Asqtad with any
  // dimension less than 4
  //----------------------------------------------------------------
  if( latt.Fclass() == F_CLASS_ASQTAD && 
      (GJP.XnodeSites()==2 ||
       GJP.YnodeSites()==2 ||
       GJP.ZnodeSites()==2 ||
       GJP.TnodeSites()==2 ) )
    ERR.General(cname,fname," RHMC force term is not implemented for Fasqtad with any dimension less than 4\n");

  // construct approximation if necessary
  generateApprox(arg);

  init();


}

//------------------------------------------------------------------
/*!

\param latt The lattice on which the HMC algorithm runs.
\param c_arg The common argument structure for all algorithms.
\param arg The algorithm parameters.
\param e_arg ?
*/
//------------------------------------------------------------------

AlgHmcRHMC::AlgHmcRHMC(Lattice& latt, CommonArg *c_arg, HmdArg *arg, 
		       EigArg *e_arg) : AlgHmd(latt, c_arg, arg)

{
  cname = "AlgHmcRHMC";
  char *fname = "AlgHmcRHMC(L&,CommonArg*,HmdArg*)";
  VRB.Func(cname,fname);

  // construct approximation if necessary
  generateApprox(arg);

  init();

  eig_arg = e_arg;
}

void AlgHmcRHMC::init()
{
  int i,j;
  cname = "AlgHmcRHMC";
  char *fname = "init()";
  VRB.Func(cname,fname);

  // Number of lattice sites
  f_sites = GJP.SnodeSites()*GJP.VolNodeSites() / (AlgLattice().FchkbEvl()+1);
  // Number of Vectors in a Vector array
  f_vec_count = f_sites * AlgLattice().SpinComponents();
  // Number of Floats in a Vector array
  f_size = f_vec_count * AlgLattice().Colors() * 2;

  // Initialize the number of dynamical fermion masses
  //----------------------------------------------------------------
  n_frm_masses = hmd_arg->n_frm_masses;
  if(n_frm_masses > MAX_HMD_MASSES){
    ERR.General(cname,fname,
		"hmd_arg->n_frm_masses = %d is larger than MAX_HMD_MASSES = %d\n",
		n_frm_masses, MAX_HMD_MASSES);
  }

  // Initialize the number of dynamical boson masses
  //----------------------------------------------------------------
  n_bsn_masses = hmd_arg->n_bsn_masses;
  if(n_bsn_masses > MAX_HMD_MASSES){
    ERR.General(cname,fname,
		"hmd_arg->n_bsn_masses = %d is larger than MAX_HMD_MASSES = %d\n",
		n_bsn_masses, MAX_HMD_MASSES);
  }

  // Allocate memory for the fermion CG arguments.
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    frm_cg_arg = (CgArg **) smalloc (cname,fname, "frm_cg_arg", 
				     n_frm_masses * sizeof(CgArg*));

    for(i=0; i<n_frm_masses; i++) {
      frm_cg_arg[i] = (CgArg *) 
	smalloc(cname,fname, "frm_cg_arg[i]",sizeof(CgArg));
    }

    // Initialize the fermion CG arguments
    //----------------------------------------------------------------
    for(i=0; i<n_frm_masses; i++){
      frm_cg_arg[i]->mass = hmd_arg->frm_mass[i];
      frm_cg_arg[i]->max_num_iter = hmd_arg->max_num_iter[i];
      frm_cg_arg[i]->stop_rsd = hmd_arg->stop_rsd_mc[i];
    }
  }
  // Allocate memory for the boson CG arguments.
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    bsn_cg_arg = (CgArg **) smalloc(n_bsn_masses * sizeof(CgArg*), 
				    cname, fname, "bsn_cg_arg");
    
    for(i=0; i<n_bsn_masses; i++){
      bsn_cg_arg[i] = (CgArg *) 
	smalloc(sizeof(CgArg), cname, fname, "bsn_cg_arg[i]");
    } 

  }

  // Initialize the boson CG arguments
  //----------------------------------------------------------------
  //??? Complete this
  for(i=0; i<n_bsn_masses; i++){
    bsn_cg_arg[i]->mass = hmd_arg->bsn_mass[i];
    bsn_cg_arg[i]->max_num_iter = hmd_arg->max_num_iter[i];
    bsn_cg_arg[i]->stop_rsd = hmd_arg->stop_rsd_mc[i];
  }

  // Allocate memory for the phi pseudo fermion field.
  //----------------------------------------------------------------
  phi = (Vector **) smalloc(n_frm_masses * sizeof(Vector*),
			    cname,fname, "phi");

  for(i=0; i<n_frm_masses; i++) {
    phi[i] = (Vector *) smalloc(f_size*sizeof(Float),cname,fname, "phi[i]");
  }      
  
  // Allocate memory for the frmn solution fermion fields.
  //----------------------------------------------------------------
  total_size = 0;
  for (i=0; i<n_frm_masses; i++) total_size += hmd_arg->FRatDeg[i];
  
  // array holding coefficents used in dwf force, folded into MInv
  alpha = (Float**) smalloc(hmd_arg->n_frm_masses*sizeof(Float*),
			    cname,fname,"alpha");
  for (i=0; i<n_frm_masses; i++)
    alpha[i] = (Float*) smalloc(f_size*hmd_arg->FRatDeg[i]*sizeof(Float),
				cname,fname,"alpha[i]");    

  if (AlgLattice().Fclass() == F_CLASS_DWF) {
    // For dwf we need the solution vector contiguous in memory
    frmn = (Vector**) smalloc(total_size*sizeof(Vector*), cname, fname, "frmn");
    
    frmn[0] = (Vector*) smalloc(f_size*total_size*sizeof(Float), 
				cname, fname, "frmn[0]");
    
    for (i=1; i<total_size; i++) frmn[i] = frmn[0] + i*f_vec_count;

    frmn_d = 0;
  } else {
    // For asqtad we need them checkerboarded with dslash applied
    frmn = (Vector**) smalloc(total_size*sizeof(Vector*), cname, fname, "frmn");    
    frmn_d = (Vector**) smalloc(total_size*sizeof(Vector*), cname, fname, "frmn_d");
    
    for (i=0; i<total_size; i++) {
      frmn[i] = (Vector*) smalloc(2*f_size*sizeof(Float));
      frmn_d[i] = frmn[i] + f_vec_count;
    }

  }	

  // Allocate memory for the boson field bsn.
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    bsn = (Vector **) smalloc(n_bsn_masses * sizeof(Vector*),
			      cname,fname, "bsn");
    for(i=0; i<n_bsn_masses; i++) {
      bsn[i] = (Vector*)smalloc(f_size*sizeof(Float), cname, fname, "bsn[i]");
    }

  }

  // Allocate memory for the initial gauge field.
  //----------------------------------------------------------------
  gauge_field_init = (Matrix *) smalloc(g_size * sizeof(Float), 
					cname,fname, "gauge_field_init");

  if (hmd_arg->reproduce != 1 && hmd_arg->reproduce != 0)
    hmd_arg->reproduce = 0;

  if (hmd_arg->reproduce) {
    
    // Allocate memory for the final gauge field.
    //----------------------------------------------------------------

    gauge_field_final = (Matrix *) smalloc(g_size * sizeof(Float), 
					   cname,fname, "gauge_field_final");
        
    // Allocate memory for the initial rng state
    //----------------------------------------------------------------
    rng4d_init = (unsigned int**) smalloc(LRG.NStates(FOUR_D)*sizeof(unsigned int*),
					  cname, fname, "rng4d_init");
    rng5d_init = (unsigned int**) smalloc(LRG.NStates()*sizeof(unsigned int*),
					  cname, fname, "rng5d_init");
    for (int i=0; i<LRG.NStates(FOUR_D); i++)
      rng4d_init[i] = (unsigned int*) smalloc(LRG.StateSize()*sizeof(unsigned int), 
					      cname, fname, "rng4d_init[i]");
    for (int i=0; i<LRG.NStates(); i++)
      rng5d_init[i] = (unsigned int*) smalloc(LRG.StateSize()*sizeof(unsigned int), 
					      cname, fname, "rng5d_init[i]");

  }

  // Used for storing residue coefficients for staggered optimisation
  //----------------------------------------------------------------
  if ((AlgLattice().Fclass() == F_CLASS_ASQTAD || AlgLattice().Fclass() == F_CLASS_P4) && total_size > 0) {
    all_res = (Float *)smalloc(total_size * sizeof(Float),
			       cname,fname, "all_res");
    Float *res = all_res;
    for (i=0; i<n_frm_masses; i++) 
      for (j=0; j<hmd_arg->FRatDeg[i]; j++)
	*(res++) = hmd_arg->FRatRes[i][j];
  }

} 

//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgHmcRHMC::~AlgHmcRHMC() {

  int i;
  char *fname = "~AlgHmcRHMC()" ;
  VRB.Func(cname,fname);

  // Free memory for the residue coefficients
  //----------------------------------------------------------------
  if ((AlgLattice().Fclass() == F_CLASS_ASQTAD || AlgLattice().Fclass() == F_CLASS_P4) && total_size > 0) 
    sfree(all_res, cname,fname, "all_res");

  // Free memory for the initial gauge field.
  //----------------------------------------------------------------
  sfree(gauge_field_init, cname,fname, "gauge_field_init");

  if (hmd_arg->reproduce) {

    // Free memory for the initial rng state
    //----------------------------------------------------------------
    for (int i=0; i<LRG.NStates(); i++)
      sfree(rng5d_init[i], cname, fname, "rng5d_init[i]");
    for (int i=0; i<LRG.NStates(FOUR_D); i++)
      sfree(rng4d_init[i], cname, fname, "rng4d_init[i]");
    
    sfree(rng5d_init, cname, fname, "rng5d_init");
    sfree(rng4d_init, cname, fname, "rng4d_init");
    // Free memory for the final gauge field.
    //----------------------------------------------------------------
    sfree(gauge_field_final, cname, fname, "gauge_field_final");

  }

  // Free memory for the boson field bsn.
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    for(i=0; i<n_bsn_masses; i++)
      sfree(bsn[i], cname,fname, "bsn[i]");
    sfree(bsn, cname,fname, "bsn");
  }

  // Free memory for the phi (pseudo fermion) fermion field.
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    for(i=0; i<n_frm_masses; i++)
      sfree(phi[i], cname,fname, "phi[i]");
    sfree(phi, cname,fname, "phi");

    // Free memory for the frmn (pseudo fermion) solution fields.
    //----------------------------------------------------------------
    if (AlgLattice().Fclass() == F_CLASS_DWF) {
      sfree(frmn[0], cname, fname, "frmn[0]");
    } else {
      sfree(frmn_d, cname,fname, "frmn_d");    
      for (i=0; i<total_size; i++) sfree(frmn[i], cname, fname, "frmn[i]");
    }
    sfree(frmn, cname, fname, "frmn");
  }

  // Free memory for alpha coefficients
  //----------------------------------------------------------------
  for (i=0; i<n_frm_masses; i++)
    sfree(alpha[i],cname,fname,"alpha[i]");
  sfree(alpha,cname,fname,"alpha");

  // Free memory for the boson CG arguments
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    for(i=0; i<n_bsn_masses; i++) {
      sfree(bsn_cg_arg[i], cname,fname, "bsn_cg_arg[i]");
    }
    sfree(bsn_cg_arg, cname,fname, "bsn_cg_arg");
  }

  // Free memory for the fermion CG arguments
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    for(i=0; i<n_frm_masses; i++){
      sfree(frm_cg_arg[i], cname,fname, "frm_cg_arg[i]");
    }
    sfree(frm_cg_arg, cname,fname, "frm_cg_arg");
  }
  
}


//------------------------------------------------------------------
//
// run(): The Rational Hybrid Monte Carlo algorithm.
/*!
  \return The probability used in the metropolis step.
  \post The following results are written to the file specified in the
  CommonArg structure:
  -# Number of molecular dynamics steps + 1
  -# Change in the hamiltonian
  -# 1/0 if accepted/rejected
  -# A measure of the change in the gauge field due to reunitarisation
  -# Another measure of the change in the gauge field due to reunitarisation
  -# Average number of solver iterations.
  -# Minimum number of solver iterations.
  -# Maximum number of solver iterations.
  -# Average solver residue
  -# Minimum solver residue
  -# Maximum solver residue
*/
//------------------------------------------------------------------
Float AlgHmcRHMC::run(void)
{
#if TARGET==cpsMPI
  using MPISCU::fprintf;
#endif
  int step;                            // Trajectory step
  Float h_init;                        // Initial Hamiltonian
  Float h_final;                       // Final Hamiltonian
  Float delta_h;                       // Final-Initial Hamiltonian
  int accept;

  int cg_iter;
  Float cg_iter_av=0;
  int cg_iter_min=0;
  int cg_iter_max=0;
  Float true_res_av=0;
  Float true_res_min=0;
  Float true_res_max=0;
  int cg_calls=0;
  int shift;
  int i, j;
  char *fname = "run()";
  char *md_time_str = "MD_time/step_size = ";
  FILE *fp;
  VRB.Func(cname,fname);
  unsigned long step_cnt = 0;
  CSM.SaveComment(step_cnt);

  Float dev = 0.0;
  Float max_diff = 0.0;

  Float trueMass=0;
  Float acceptance;                            // The acceptance probability
  Float efficiency;

 
  // Get the Lattice object
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // Set the microcanonical time step
  //----------------------------------------------------------------
  Float dt = hmd_arg->step_size;

  if(hmd_arg->metropolis) {
    // Save initial gauge field configuration
    //--------------------------------------------------------------
    lat.CopyGaugeField(gauge_field_init);
  }

  Float delta_h0;                       // save first run dH
  int Ntests;
  LRGState lrg_state;

  if (hmd_arg->reproduce) {
    // Save initial rng state
    //--------------------------------------------------------------
//    LRG.GetStates(rng5d_init, FIVE_D);
//    LRG.GetStates(rng4d_init, FOUR_D);
    lrg_state.GetStates();
    
    
    // Need to save initial lattice if have not already done so
    if (!hmd_arg->metropolis)
      lat.CopyGaugeField(gauge_field_init);


    if (hmd_arg->reproduce_attempt_limit < 1 ||
	hmd_arg->reproduce_attempt_limit > 5)
      hmd_arg->reproduce_attempt_limit = 3;
    Ntests = 2;
  } else {
    hmd_arg->reproduce_attempt_limit = 1;
    Ntests = 1;
  }

  if (lat.Fclass() == F_CLASS_DWF || lat.Fclass() == F_CLASS_ASQTAD || lat.Fclass() == F_CLASS_P4) {
    for (int i=0; i<n_frm_masses; i++)
      for (int j=0; j<hmd_arg->FRatDeg[i]; j++)
	alpha[i][j] = sqrt(hmd_arg->FRatRes[i][j]);
  }

  GDS.Set(0,0,0,0);
  GDS.SetOrigin(0,0,0,0);
  // Try attempt_limit times to generate the same final gauge config 
  // consecutively
  for (int attempt = 0; attempt < hmd_arg->reproduce_attempt_limit; attempt++) {

    for (int test = 0; test < Ntests; test++) {
      VRB.Result(cname,fname,"Running test %d of %d, attempt %d of %d\n", 
		 test+1, Ntests, attempt+1, hmd_arg->reproduce_attempt_limit);

      // Restore the initial gauge field and rng state

      if ( !(test == 0 && attempt ==0) ) {
	lat.GaugeField(gauge_field_init);
//	LRG.SetStates(rng4d_init, FOUR_D);
//	LRG.SetStates(rng5d_init, FIVE_D);
        lrg_state.SetStates();
        GDS.Set(1,1,1,0);
        LRG.Shift();
        lat.Shift();
        GDS.SetOrigin(1,1,1,0);
      }

      // Initialize Hamiltonian variables
      //----------------------------------------------------------------
      h_init=0;                        // Initial Hamiltonian
      h_final=0;                       // Final Hamiltonian
      delta_h=0;                       // Final-Initial Hamiltonian
      
      // Initialize monitor variables
      //----------------------------------------------------------------
      cg_iter_av   = 0.0;
      cg_iter_min  = 1000000;
      cg_iter_max  = 0;
      true_res_av  = 0.0;
      true_res_min = 3.4e+38;
      true_res_max = 0.0;
      cg_calls     = 0;
      
      // reset MD time in Lattice
      //----------------------------------------------------------------
      lat.MdTime(0.0);
      VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
      
      
      // Heat Bath for the conjugate momenta
      //----------------------------------------------------------------
      lat.RandGaussAntiHermMatrix(mom, 1.0);
      
      // Heat Bath for the boson field bsn
      //----------------------------------------------------------------
      for(i=0; i<n_bsn_masses; i++){
	lat.RandGaussVector(frmn[0], 0.5, Ncb);
	lat.RandGaussVector(bsn[i], 0.5, Ncb);
	lat.SetPhi(frmn[1], frmn[0], bsn[i], hmd_arg->bsn_mass[i]);
	lat.RandGaussVector(bsn[i], 0.5, Ncb);
	lat.FmatEvlInv(bsn[i], frmn[1], bsn_cg_arg[i], CNV_FRM_NO);
      }
      
      // Heat Bath for the pseudo-fermions (phi)
      // Use the SI approximation since need the inverse of the action
      //----------------------------------------------------------------
      
      h_init = lat.GhamiltonNode() + lat.MomHamiltonNode(mom);
      for(i=0; i<n_frm_masses; i++){
    
	lat.RandGaussVector(frmn[i], 0.5, Ncb);
	h_init += lat.FhamiltonNode(frmn[i],frmn[i]);
	
	massRenormalise(&(frm_cg_arg[i]->mass), &trueMass, hmd_arg->SRatDeg[i], 
			hmd_arg->SIRatPole[i], RENORM_FORWARDS);

	phi[i] -> CopyVec(frmn[i], f_size);
	phi[i] -> VecTimesEquFloat(hmd_arg->SIRatNorm[i], f_size);
	cg_iter = lat.FmatEvlMInv(phi+i, frmn[i], hmd_arg->SIRatPole[i], hmd_arg->SRatDeg[i],
				  hmd_arg->isz, frm_cg_arg[i], CNV_FRM_NO, SINGLE, 
				  hmd_arg->SIRatRes[i]);
	
	massRenormalise(&(frm_cg_arg[i]->mass), &trueMass, hmd_arg->SRatDeg[i], 
			hmd_arg->SIRatPole[i], RENORM_BACKWARDS);
	
	cg_iter_av += cg_iter;
	if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
	if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
	true_res_av += frm_cg_arg[i]->true_rsd;
	if(frm_cg_arg[i]->true_rsd < true_res_min) true_res_min = frm_cg_arg[i]->true_rsd;
	if(frm_cg_arg[i]->true_rsd > true_res_max) true_res_max = frm_cg_arg[i]->true_rsd;
	cg_calls++;      

      }

      // Calculate initial boson contribution to the Hamiltonian
      //---------------------------------------------------------------
      for(i=0; i<n_bsn_masses; i++)
	h_init += lat.BhamiltonNode(bsn[i], hmd_arg->bsn_mass[i]);

      //----------------------------------------------------------------
      // Molecular Dynamics Trajectory
      //----------------------------------------------------------------

      // Reset the residual error for the Molecular Dynamics
      //--------------------------------------------------------------
      for (i=0; i<n_frm_masses; i++)
	frm_cg_arg[i]->stop_rsd = hmd_arg->stop_rsd_md[i];

      // Perform initial QPQ integration
      //--------------------------------------------------------------
      if (hmd_arg->steps_per_traj > 0) {
	lat.EvolveGfield(mom, 0.25*dt/(Float)hmd_arg->sw);
	for (i=0; i<hmd_arg->sw; i++) {
	  for(j=0; j<n_bsn_masses; j++)
	    lat.EvolveMomFforce(mom, bsn[j], hmd_arg->bsn_mass[j], 
				-0.5*dt/(Float)hmd_arg->sw);
	  lat.EvolveMomGforce(mom, 0.5*dt/(Float)hmd_arg->sw);
	  if (i < hmd_arg->sw-1) lat.EvolveGfield(mom, 0.5*dt/(Float)hmd_arg->sw);
	  else lat.EvolveGfield(mom, 0.25*dt/(Float)hmd_arg->sw);
	}
      }

      // First leap frog step has occurred. Increment MD Time clock in
      // Lattice one half time step.
      //----------------------------------------------------------------
      lat.MdTimeInc(0.5) ;
      VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));

      // Run through the trajectory
      //----------------------------------------------------------------

      for(step=0; step < hmd_arg->steps_per_traj; step++){
	CSM.SaveComment(++step_cnt);


	// Evolve momenta by one step using the fermion force
	//--------------------------------------------------------------

	shift = 0;
	for(i=0; i<n_frm_masses; i++){

	  massRenormalise(&(frm_cg_arg[i]->mass), &trueMass, hmd_arg->FRatDeg[i], 
			  hmd_arg->FRatPole[i], RENORM_FORWARDS);
      
	  for (j=0; j<hmd_arg->FRatDeg[i]; j++) 
	    bzero((char *)frmn[j+shift],f_size*sizeof(Float));

      	  cg_iter = lat.FmatEvlMInv(frmn + shift, phi[i], hmd_arg->FRatPole[i], 
				    hmd_arg->FRatDeg[i], hmd_arg->isz, frm_cg_arg[i], 
				    CNV_FRM_NO, frmn_d + shift);

	  massRenormalise(&(frm_cg_arg[i]->mass), &trueMass, hmd_arg->FRatDeg[i], 
			  hmd_arg->FRatPole[i], RENORM_BACKWARDS);
      
	  cg_iter_av += cg_iter;
	  if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
	  if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
	  true_res_av += frm_cg_arg[i]->true_rsd;
	  if(frm_cg_arg[i]->true_rsd < true_res_min) 
	    true_res_min = frm_cg_arg[i]->true_rsd;
	  if(frm_cg_arg[i]->true_rsd > true_res_max) 
	    true_res_max = frm_cg_arg[i]->true_rsd;
	  cg_calls++;      

	  if ((lat.Fclass() != F_CLASS_ASQTAD && lat.Fclass() != F_CLASS_P4))
	    lat.RHMC_EvolveMomFforce(mom, frmn+shift, hmd_arg->FRatDeg[i], 0,
				     hmd_arg->FRatRes[i], hmd_arg->frm_mass[i],
				     dt, frmn_d+shift, FORCE_MEASURE_NO);

	  shift += hmd_arg->FRatDeg[i];
	}

	// Only for the case of asqtad fermions do we perform this optimisation
	//--------------------------------------------------------------      
	if ((lat.Fclass() == F_CLASS_ASQTAD || lat.Fclass() == F_CLASS_P4) && n_frm_masses != 0)
	  {
	    VRB.Flow(cname,fname,"start EvolveMomFforce\n");
#if 1
#if 1
	    lat.RHMC_EvolveMomFforce(mom, frmn, total_size, 0, all_res, 0.0, dt, frmn_d, FORCE_MEASURE_NO);
#else
	    for(int jj = 0; jj< total_size; jj++)
	      {
		VRB.Flow(cname,fname,"before shift = %d, **(frmn+jj) = %e, all_res[jj]=%e\n",jj,*((IFloat *) *(frmn+jj)), all_res[jj]);
		lat.RHMC_EvolveMomFforce(mom, frmn+jj, 1, 0, all_res+jj, 0.0, dt, frmn_d+jj, FORCE_MEASURE_NO);
		VRB.Flow(cname,fname,"after shift = %d, **(frmn+jj) = %e, all_res[jj]=%e\n",jj,*((IFloat *) *(frmn+jj)), all_res[jj]);
	      }
#endif
#else
	    for(int jj = 0; jj < total_size; jj++)
	      {
		//VRB.Flow(cname,fname,"before shift = %d, **(frmn+jj) = %e, all_res[jj]=%e\n",jj,*((IFloat *) *(frmn+jj)), all_res[jj]);
		frmn[jj]->VecTimesEquFloat(all_res[jj], GJP.VolNodeSites()*lat.FsiteSize()/2);
		//VRB.Flow(cname,fname,"shift = %d, **(frmn+jj) = %e\n",jj,*((IFloat *) *(frmn+jj)));
		lat.EvolveMomFforce(mom, *(frmn+jj), 0.0, dt);
		//VRB.Flow(cname,fname,"after shift = %d, **(frmn+jj) = %e, all_res[jj]=%e\n",jj,*((IFloat *) *(frmn+jj)), all_res[jj]);
	      }
#endif
	    VRB.Flow(cname,fname,"*((IFloat *)mom+777) = %e\n", *((IFloat *) mom+777));
	  }

	// Increment MD Time clock in Lattice by one half time step.
	//--------------------------------------------------------------
	lat.MdTimeInc(0.5);
	VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));


	if (step < hmd_arg->steps_per_traj-1) 
	  {
	    // Perform QPQ integration (pure gauge, boson and momenta)
	    //----------------------------------------------------------------
	    lat.EvolveGfield(mom, 0.5*dt/(Float)hmd_arg->sw);
	    for (i=0; i<hmd_arg->sw; i++) {
	      for(j=0; j<n_bsn_masses; j++)
		lat.EvolveMomFforce(mom, bsn[j], hmd_arg->bsn_mass[j], 
				    -dt/(Float)hmd_arg->sw);
	      lat.EvolveMomGforce(mom, dt/(Float)hmd_arg->sw);
	      if (i < hmd_arg->sw-1) lat.EvolveGfield(mom, dt/(Float)hmd_arg->sw);
	      else lat.EvolveGfield(mom, 0.5*dt/(Float)hmd_arg->sw);
	    }

	    // Increment MD Time clock in Lattice by one half time step.
	    //--------------------------------------------------------------
	    lat.MdTimeInc(0.5);
	    VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
	  } 
	else 
	  {
	    // Perform final QPQ integration
	    //----------------------------------------------------------------
	    lat.EvolveGfield(mom, 0.25*dt/(Float)hmd_arg->sw);
	    for (i=0; i<hmd_arg->sw; i++) {
	      for(j=0; j<n_bsn_masses; j++)
		lat.EvolveMomFforce(mom, bsn[j], hmd_arg->bsn_mass[j], 
				    -0.5*dt/(Float)hmd_arg->sw);
	      lat.EvolveMomGforce(mom, 0.5*dt/(Float)hmd_arg->sw);
	      if (i < hmd_arg->sw-1) lat.EvolveGfield(mom, 0.5*dt/(Float)hmd_arg->sw);
	      else lat.EvolveGfield(mom, 0.25*dt/(Float)hmd_arg->sw);
	    }

	  }
      }

      CSM.SaveComment(++step_cnt);
      // Reunitarize
      //----------------------------------------------------------------
      if(hmd_arg->reunitarize == REUNITARIZE_YES){
	lat.Reunitarize(dev, max_diff);
      }
      h_final = lat.GhamiltonNode();
      IFloat h_gauge = h_final;

      // Measure bounds and generate new approximations
      if (hmd_arg->approx_type == DYNAMIC) {
#ifndef GMP
	ERR.General(cname,fname,"Dynamical rational generation not possible without gmp\n");
#else
	dynamicalApprox();
#endif
      }

      // Calculate final fermion contribution to the Hamiltonian
      //---------------------------------------------------------------
      shift = 0;
      for (i=0; i<n_frm_masses; i++) {
	massRenormalise(&(frm_cg_arg[i]->mass), &trueMass, hmd_arg->SRatDeg[i], 
			hmd_arg->SRatPole[i], RENORM_FORWARDS);

	// Reset the residual error to the mc error
	frm_cg_arg[i]->stop_rsd = hmd_arg->stop_rsd_mc[i];

	frmn[shift] -> CopyVec(phi[i],f_size);
	frmn[shift] -> VecTimesEquFloat(hmd_arg->SRatNorm[i], f_size);

	cg_iter = lat.FmatEvlMInv(frmn+shift, phi[i], hmd_arg->SRatPole[i], 
				  hmd_arg->SRatDeg[i], hmd_arg->isz, frm_cg_arg[i], 
				  CNV_FRM_NO, SINGLE, hmd_arg->SRatRes[i]);

	massRenormalise(&(frm_cg_arg[i]->mass), &trueMass, hmd_arg->SRatDeg[i], 
			hmd_arg->SRatPole[i], RENORM_BACKWARDS);

	cg_iter_av += cg_iter;
	if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
	if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
	true_res_av += frm_cg_arg[i]->true_rsd;
	if(frm_cg_arg[i]->true_rsd < true_res_min) true_res_min = frm_cg_arg[i]->true_rsd;
	if(frm_cg_arg[i]->true_rsd > true_res_max) true_res_max = frm_cg_arg[i]->true_rsd;
	cg_calls++;      

	h_final += lat.FhamiltonNode(frmn[shift], frmn[shift]);
	shift += hmd_arg->FRatDeg[i];

      }
  
      IFloat h_fermion = h_final - h_gauge;
      h_final += lat.MomHamiltonNode(mom);
      IFloat h_mom = h_final - h_gauge - h_fermion;

      // Calculate final boson contribution to the Hamiltonian
      //---------------------------------------------------------------
      for(i=0; i<n_bsn_masses; i++)
	h_final += lat.BhamiltonNode(bsn[i], hmd_arg->bsn_mass[i]);

      // Calculate Final-Initial Hamiltonian 
      //---------------------------------------------------------------
      delta_h = h_final - h_init;
      glb_sum(&h_gauge);
      glb_sum(&h_fermion);
      glb_sum(&h_mom);
      glb_sum(&h_init);
      glb_sum(&h_final);
      glb_sum(&delta_h);
      VRB.Result(cname,fname,"hamilton_gauge = %e\n", h_gauge);
      VRB.Result(cname,fname,"hamilton_fermion = %e\n", h_fermion);
      VRB.Result(cname,fname,"hamilton_mom = %e\n", h_mom);
      VRB.Result(cname,fname,"hamilton_final = %e\n", h_final);
      VRB.Result(cname,fname,"hamilton_initial = %e\n", h_init);
      VRB.Result(cname,fname,"delta_hamilton = %e\n", delta_h);

      // Check that delta_h is the same accross all s-slices 
      // (relevant only if GJP.Snodes() != 1)
      //----------------------------------------------------------------
      if(GJP.Snodes() != 1) {
	VRB.Flow(cname,fname, "Checking Delta H across s-slices\n");
	lat.SoCheck(delta_h);
      }

      if (hmd_arg->reproduce) {
	// Save final gauge field configuration
        if ( test != 0) {
          GDS.Set(-1,-1,-1,0);
          LRG.Shift();
          lat.Shift();
          GDS.SetOrigin(0,0,0,0);
        }
        else {
	  lat.CopyGaugeField(gauge_field_final);
	  delta_h0 = delta_h;
        }
      }
    } // end test loop
  
    if (hmd_arg->reproduce) {  
      // Compare the final gauge configs generated
      // if passed continue, else try again
      
      if (lat.CompareGaugeField(gauge_field_final)) {
	VRB.Result(cname,fname,"Passed reproducibility test\n");
	break;
      } else {
	VRB.Result(cname,fname, 
		   "Failed reproducibility dH0 = %18.16e, dH1 = %18.16e, delta dH = %18.16e\n", 
		   delta_h0, delta_h, delta_h-delta_h0);
      }

      VRB.Result(cname,fname,"Failed reproducibility test %d\n",attempt);
      if (attempt == hmd_arg->reproduce_attempt_limit-1) 
	ERR.General(cname,fname,"Failed to reproduce\n");
    } else {
      break;
    }

  } // end attempt loop

  // Metropolis step
  //---------------------------------------------------------------
  if(hmd_arg->metropolis){
    // Metropolis accept-reject step
    //--------------------------------------------------------------
    accept = lat.MetropolisAccept(delta_h,&acceptance);
    if( !(accept) ){
      // Trajectory rejected
      //------------------------------------------------------------
      lat.GaugeField(gauge_field_init);
      VRB.Result(cname,fname,"Metropolis step -> Rejected\n");
    }
    else {
      // Trajectory accepted. 
      // Increment the Gauge Update counter in Lattice.
      //-------------------------------------------------------------
      lat.GupdCntInc(1);
      VRB.Result(cname,fname,"Metropolis step -> Accepted\n");
    }
  } 
  else {
    accept = 1;
    acceptance = 1.0;
    lat.GupdCntInc(1);
    VRB.Result(cname,fname,"No Metropolis step -> Accepted\n");
  }    


  // Calculate average of monitor variables
  //---------------------------------------------------------------
  
  efficiency = Float(acceptance) / Float(cg_iter_av);
  cg_iter_av = Float(cg_iter_av) / Float(cg_calls);
  true_res_av = Float(true_res_av) / Float(cg_calls);


  // If GJP.Snodes() !=1  the gauge field is spread out
  // across s-slices of processors. It must be identical
  // on each slice. Check to make sure and exit if it
  // is not identical. A case where this is relevant
  // is the DWF spread-out case.
  //----------------------------------------------------------------
  if(GJP.Snodes() != 1) {
    VRB.Flow(cname,fname, "Checking gauge field across s-slices\n");
    lat.GsoCheck();
  }

  // Print out monitor info
  //---------------------------------------------------------------
  if(common_arg->results != 0){
    if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    Fprintf(fp,"%d %e %e %d %e %e %e %d %d %e %e %e\n",
	    hmd_arg->steps_per_traj,
	    IFloat(delta_h), 
	    acceptance,
	    accept, 
	    IFloat(dev),
	    IFloat(max_diff),
	    IFloat(cg_iter_av),
	    cg_iter_min,
	    cg_iter_max,
	    IFloat(true_res_av),
	    IFloat(true_res_min),
	    IFloat(true_res_max));
    Fclose(fp);
  }

  VRB.Result(cname,fname,
	     "Hmc steps = %d, Delta_hamilton = %e, accept = %d, dev = %e, max_diff = %e\n",
	     hmd_arg->steps_per_traj,
	     IFloat(delta_h), 
	     accept,
	     IFloat(dev),
	     IFloat(max_diff));
  VRB.Result(cname,fname,
	     "CG iterations: average = %e, min = %d, max = %d\n",
	     IFloat(cg_iter_av),
	     cg_iter_min,
	     cg_iter_max);
  VRB.Result(cname,fname,
	     "True Residual / |source|: average = %e, min = %e, max = %e\n",
	     IFloat(true_res_av),
	     IFloat(true_res_min),
	     IFloat(true_res_max));

  VRB.Result(cname,fname,
	     "Efficiency (acceptance per cg iteration) = %e\n",
	     efficiency);

  VRB.Result(cname,fname,"Configuration number = %d\n", lat.GupdCnt());

  for (int i=0; i<n_frm_masses; i++)
    VRB.Result(cname,fname,"%d psi_bar psi = %e\n", 
	       i, lat.FhamiltonNode(phi[i],phi[i]));
  

  // Reset Molecular Dynamics time counter
  //----------------------------------------------------------------
  lat.MdTime(0.0);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));

  //  ERR.HdwCheck(cname,fname);

  return acceptance;
}

/*!
  Renormalise the smallest shift into the mass parameters - reduces linear algebra
  in the multi-mass solver.
*/
void AlgHmcRHMC::massRenormalise(Float *mass, Float *trueMass, int degree, 
				  Float *shift, MassRenormaliseDir direction) {

    // Can only renormalise mass for staggered or asqtad cases
  if (AlgLattice().Fclass() == F_CLASS_ASQTAD || AlgLattice().Fclass() == F_CLASS_STAG || AlgLattice().Fclass() == F_CLASS_P4) {
    if (direction == RENORM_FORWARDS) {
      *trueMass = *mass;
      *mass = sqrt((*trueMass)*(*trueMass) + shift[hmd_arg->isz]/4.0);
      Float zeroPole = shift[hmd_arg->isz];
      for (int j=0; j<degree; j++) shift[j] -= zeroPole;
    } else if (direction == RENORM_BACKWARDS) {
      Float zeroPole = 4.0*((*mass)*(*mass) - (*trueMass)*(*trueMass));
      for (int j=0; j<degree; j++) shift[j] += zeroPole;
      *mass = *trueMass;
    }
  }

}

/*!
  Generate the optimal rational approximation for the RHMC simulation
*/
void AlgHmcRHMC::generateApprox(HmdArg *hmd_arg)
{

  // Construct approximations
  for (int i=0; i<hmd_arg->n_frm_masses; i++) {
    for (int j=0; j<i; j++) {
      // no need to recalculate approximation if same mass
      if (hmd_arg->frm_mass[j] == hmd_arg->frm_mass[i] &&
	  hmd_arg->field_type[j] == hmd_arg->field_type[i] &&
	  hmd_arg->frm_power_num[j] ==hmd_arg-> frm_power_num[i] &&
	  hmd_arg->frm_power_den[j] ==hmd_arg-> frm_power_den[i] ) {
	hmd_arg->FRatDeg[i] = hmd_arg->FRatDeg[j];
	hmd_arg->FRatNorm[i] = hmd_arg->FRatNorm[j];
	for (int k=0; k<hmd_arg->FRatDeg[i]; k++) {
	  hmd_arg->FRatRes[i][k] = hmd_arg->FRatRes[j][k];
	  hmd_arg->FRatPole[i][k] = hmd_arg->FRatPole[j][k];
	}
	hmd_arg->SRatDeg[i] = hmd_arg->SRatDeg[j];
	hmd_arg->SRatNorm[i] = hmd_arg->SRatNorm[j];
	hmd_arg->SIRatNorm[i] = hmd_arg->SIRatNorm[j];
	for (int k=0; k<hmd_arg->SRatDeg[i]; k++) {
	  hmd_arg->SRatRes[i][k] = hmd_arg->SRatRes[j][k];
	  hmd_arg->SRatPole[i][k] = hmd_arg->SRatPole[j][k];
	  hmd_arg->SIRatRes[i][k] = hmd_arg->SIRatRes[j][k];
	  hmd_arg->SIRatPole[i][k] = hmd_arg->SIRatPole[j][k];
	}
	hmd_arg->valid_approx[i] = 1;
      }
    }

// CJ: backward compatibility

    if (!hmd_arg->valid_approx[i]) {
//      AlgRemez remez(hmd_arg->lambda_low[i],hmd_arg->lambda_high[i],hmd_arg->precision);
//      hmd_arg->FRatError[i] = remez.generateApprox(hmd_arg->FRatDeg[i],hmd_arg->frm_power_num[i], hmd_arg->frm_power_den[i]);

        RemezArg remez_arg;
	//-----------------------------------------------
	//Ugly hack to make things work with this deprecated style of RHMC
	//michaelc 03/01/06
	remez_arg.approx_type = RATIONAL_APPROX_POWER;
	//-----------------------------------------------
        remez_arg.degree = hmd_arg->FRatDeg[i];
        remez_arg.field_type = hmd_arg->field_type[i];
        remez_arg.lambda_low = hmd_arg->lambda_low[i];
        remez_arg.lambda_high = hmd_arg->lambda_high[i];
        remez_arg.power_num = hmd_arg->frm_power_num[i];
        remez_arg.power_den = hmd_arg->frm_power_den[i];
        remez_arg.precision = hmd_arg->precision;
	remez_arg.delta_m = 0.0;
        {
            AlgRemez remez(remez_arg);
            remez.generateApprox();
            hmd_arg->FRatError[i] = remez_arg.error;
            if (hmd_arg->field_type[i] == BOSON) {
                remez.getPFE(hmd_arg->FRatRes[i],hmd_arg->FRatPole[i],&hmd_arg->FRatNorm[i]);
             } else {
                remez.getIPFE(hmd_arg->FRatRes[i],hmd_arg->FRatPole[i],&hmd_arg->FRatNorm[i]);
             }
        }

   //   hmd_arg->SRatError[i] = remez.generateApprox(hmd_arg->SRatDeg[i],hmd_arg->frm_power_num[i], 2*hmd_arg->frm_power_den[i]);
        remez_arg.degree = hmd_arg->SRatDeg[i];
        remez_arg.power_num = hmd_arg->frm_power_num[i];
        remez_arg.power_den = 2*hmd_arg->frm_power_den[i];
        {
            AlgRemez remez(remez_arg);
            remez.generateApprox();
            hmd_arg->SRatError[i] = remez_arg.error;

            if (hmd_arg->field_type[i] == BOSON) {
                remez.getPFE(hmd_arg->SRatRes[i],hmd_arg->SRatPole[i],&hmd_arg->SRatNorm[i]);
                remez.getIPFE(hmd_arg->SIRatRes[i],hmd_arg->SIRatPole[i],&hmd_arg->SIRatNorm[i]);
            }else {
                remez.getIPFE(hmd_arg->SRatRes[i],hmd_arg->SRatPole[i],&hmd_arg->SRatNorm[i]);
                remez.getPFE(hmd_arg->SIRatRes[i],hmd_arg->SIRatPole[i],&hmd_arg->SIRatNorm[i]);
            }
        }
        hmd_arg->valid_approx[i] = 1;
    }      
    
  }
  
}

/*!
  Measures the eigenvalue bounds and regenerates the approximation
  if necessary.
*/
void AlgHmcRHMC::dynamicalApprox()
{

  char *fname = "dynamicalApprox()";
  VRB.Func(cname,fname);

  // Get the Lattice object
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // First measure highest and lowest eigenvalues
  eig_arg->N_eig = 1;

  Vector **psi = (Vector**) smalloc(sizeof(Vector*),
				    cname,fname, "psi");
  Float *lambda = (Float*) smalloc(2*sizeof(Float), 
				   cname,fname, "lambda");
  Float *chirality = (Float*) smalloc(sizeof(Float),
				      cname,fname, "chirality");
  int *valid_eig = (int*) smalloc(sizeof(int), 
				  cname,fname, "valid_eig");
  Float **hsum = (Float**) smalloc(sizeof(Float*), 
				   cname,fname, "hsum");

  size_t hsum_size=0;
  switch(eig_arg->hsum_dir == 0){
  case DIR_X:
    hsum_size = GJP.Xnodes()*GJP.XnodeSites() * sizeof(Float);
    break;
  case DIR_Y:
    hsum_size = GJP.Ynodes()*GJP.YnodeSites() * sizeof(Float);
    break;
  case DIR_Z:
    hsum_size = GJP.Znodes()*GJP.ZnodeSites() * sizeof(Float);
    break;
  case DIR_T:
    hsum_size = GJP.Tnodes()*GJP.TnodeSites() * sizeof(Float);
    break;
  default:
    ERR.General(cname, fname, "Bad value %d for eig_arg.hsum_dir\n",
		eig_arg->hsum_dir);
  }
  hsum[0] = (Float*) smalloc(hsum_size, cname,fname, "hsum[0]");
  psi[0] = (Vector*) smalloc(f_size*sizeof(Float), cname,fname, "psi[0]");

  for (int i=0; i<n_frm_masses; i++) {
    // Reset the required degree of approximation
    if (hmd_arg->FRatDegNew[i] == 0)
      hmd_arg->FRatDegNew[i] = hmd_arg->FRatDeg[i];

    if (hmd_arg->SRatDegNew[i] == 0)
      hmd_arg->SRatDegNew[i] = hmd_arg->SRatDeg[i];

    eig_arg->mass = hmd_arg->frm_mass[i];

    // Measure lowest e-value
    eig_arg->RitzMatOper = MATPCDAG_MATPC;
    lat.RandGaussVector(psi[0], 0.5, 1);
    lat.FeigSolv(psi, lambda, chirality, valid_eig, hsum, eig_arg, CNV_FRM_NO);
      
    // Measure highest e-value
    eig_arg->RitzMatOper = NEG_MATPCDAG_MATPC;
    lat.RandGaussVector(psi[0], 0.5, 1);
    lat.FeigSolv(psi, lambda+1, chirality, valid_eig, hsum, eig_arg,CNV_FRM_NO);
    lambda[1] *= -1;

    VRB.Result(cname,fname, "Old Spectral Bounds are [%f,%f]\n",
	       hmd_arg->lambda_low[i],hmd_arg->lambda_high[i]);
    VRB.Result(cname,fname, "New Spectral Bounds are [%f,%f]\n",
	       lambda[0],lambda[1]);

    // Need to reconstruct approximations if either the spectrum leaves the 
    // interval or if it moves much inside the interval
    
    Float delta = eig_arg->Rsdlam + hmd_arg->spread;
    
    if (lambda[0] < hmd_arg->lambda_low[i]  || // below interval
	lambda[1] > hmd_arg->lambda_high[i] || // above interval
	lambda[0] > hmd_arg->lambda_low[i]  * (1.0+delta) || //inside from below
	lambda[1] < hmd_arg->lambda_high[i] * (1.0-delta) ||
	hmd_arg->FRatDeg[i] != hmd_arg->FRatDegNew[i] ||
	hmd_arg->SRatDeg[i] != hmd_arg->SRatDegNew[i] ) { // inside from above
      VRB.Result(cname,fname,"Reconstructing approximation\n");
      
      if (hmd_arg->SRatDeg[i] != hmd_arg->SRatDegNew[i])
	hmd_arg->SRatDeg[i] = hmd_arg->SRatDegNew[i];
      
      // If bounded, remain within bounds
      if (lambda[0]*(1-delta) > hmd_arg->lambda_min[i]) 
	hmd_arg->lambda_low[i] = lambda[0]*(1-delta);
      else
	hmd_arg->lambda_low[i] = hmd_arg->lambda_min[i];
    
      if (lambda[1]*(1+delta) < hmd_arg->lambda_max[i]) 
	hmd_arg->lambda_high[i] = lambda[1]*(1+delta);
      else
	hmd_arg->lambda_high[i] = hmd_arg->lambda_max[i];

//      AlgRemez remez(hmd_arg->lambda_low[i], hmd_arg->lambda_high[i], hmd_arg->precision);
      
//      remez.generateApprox(hmd_arg->FRatDeg[i],hmd_arg->frm_power_num[i], hmd_arg->frm_power_den[i]);
      RemezArg remez_arg;
      remez_arg.degree = hmd_arg->FRatDeg[i];
      remez_arg.field_type = hmd_arg->field_type[i];
      remez_arg.lambda_low = hmd_arg->lambda_low[i];
      remez_arg.lambda_high = hmd_arg->lambda_high[i];
      remez_arg.power_num = hmd_arg->frm_power_num[i];
      remez_arg.power_den = hmd_arg->frm_power_den[i];
      remez_arg.precision = hmd_arg->precision;
      {
        AlgRemez remez(remez_arg);
        remez.generateApprox();
        hmd_arg->FRatError[i] = remez_arg.error;
        if (hmd_arg->field_type[i] == BOSON) {
          remez.getPFE(hmd_arg->FRatRes[i], hmd_arg->FRatPole[i], &hmd_arg->FRatNorm[i]);
        } else {
          remez.getIPFE(hmd_arg->FRatRes[i], hmd_arg->FRatPole[i], &hmd_arg->FRatNorm[i]);
        }
      }

//      remez.generateApprox(hmd_arg->SRatDeg[i],hmd_arg->frm_power_num[i], 2*hmd_arg->frm_power_den[i]);

      remez_arg.degree = hmd_arg->SRatDeg[i];
      remez_arg.power_num = hmd_arg->frm_power_num[i];
      remez_arg.power_den = 2*hmd_arg->frm_power_den[i];
      {
        AlgRemez remez(remez_arg);
        remez.generateApprox();
        hmd_arg->SRatError[i] = remez_arg.error;
        if (hmd_arg->field_type[i] == BOSON) {
          remez.getPFE(hmd_arg->SRatRes[i],hmd_arg->SRatPole[i], &hmd_arg->SRatNorm[i]);
          remez.getIPFE(hmd_arg->SIRatRes[i],hmd_arg->SIRatPole[i], &hmd_arg->SIRatNorm[i]);
        } else {
          remez.getIPFE(hmd_arg->SRatRes[i],hmd_arg->SRatPole[i], &hmd_arg->SRatNorm[i]);
          remez.getPFE(hmd_arg->SIRatRes[i],hmd_arg->SIRatPole[i], &hmd_arg->SIRatNorm[i]);
        }
      }
    } else {
      VRB.Result(cname,fname,"Reconstruction not necessary\n");
    }
    
  }

  sfree(hsum[0], cname,fname, "hsum[0]");
  sfree(hsum, cname,fname, "hsum");
  sfree(psi[0], cname,fname, "psi[0]");
  sfree(psi, cname,fname, "psi");
  sfree(lambda, cname,fname, "lambda");
  sfree(chirality, cname,fname, "chirality");
  sfree(valid_eig, cname,fname, "valid_eig");

}


CPS_END_NAMESPACE
