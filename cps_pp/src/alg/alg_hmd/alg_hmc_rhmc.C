#include<config.h>
CPS_START_NAMESPACE 

//------------------------------------------------------------------
//
// alg_hmc_rhmc.C
//
// AlgHmcRHMC is derived from Alg and is relevant to the Rational
// Hybrid Monte Carlo Algorithm.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<stdio.h>
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
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
AlgHmcRHMC::AlgHmcRHMC(Lattice& latt, 
		     CommonArg *c_arg,
		     HmdArg *arg) : 
		     AlgHmd(latt, c_arg, arg) 
{
  cname = "AlgHmcRHMC";
  char *fname = "AlgHmcRHMC(L&,CommonArg*,HmdArg*)";
  VRB.Func(cname,fname);

  // Calculate the fermion field size.
  //----------------------------------------------------------------
  f_size = GJP.VolNodeSites() * latt.FsiteSize() / (latt.FchkbEvl()+1);


  // Approx type must be constant if this constructor is called
  //----------------------------------------------------------------
  hmd_arg->approx_type = CONSTANT;

  init();
}

//------------------------------------------------------------------
// Constructor for when dynamical approximations are used
//------------------------------------------------------------------

AlgHmcRHMC::AlgHmcRHMC(Lattice& latt, CommonArg *c_arg, HmdArg *arg, 
		       EigArg *e_arg) : AlgHmd(latt, c_arg, arg)

{
  cname = "AlgHmcRHMC";
  char *fname = "AlgHmcRHMC(L&,CommonArg*,HmdArg*)";
  VRB.Func(cname,fname);

  // Calculate the fermion field size.
  //----------------------------------------------------------------
  f_size = GJP.VolNodeSites() * latt.FsiteSize() / (latt.FchkbEvl()+1);

  init();
  eig_arg = e_arg;
}

void AlgHmcRHMC::init()
{
  int i,j;
  cname = "AlgHmcRHMC";
  char *fname = "init()";
  VRB.Func(cname,fname);
  int n_masses;

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
    frm_cg_arg = (CgArg **) smalloc(n_frm_masses * sizeof(CgArg*));
    if(frm_cg_arg == 0)
      ERR.Pointer(cname,fname, "frm_cg_arg");
    VRB.Smalloc(cname,fname,
		"frm_cg_arg",frm_cg_arg, n_frm_masses * sizeof(CgArg*));
    
    for(i=0; i<n_frm_masses; i++){
      frm_cg_arg[i] = (CgArg *) smalloc(sizeof(CgArg));
      if(frm_cg_arg[i] == 0)
	ERR.Pointer(cname,fname, "frm_cg_arg[i]");
      VRB.Smalloc(cname,fname,
		  "frm_cg_arg[i]", frm_cg_arg[i], sizeof(CgArg));
    } 
  }

  // Initialize the fermion CG arguments
  //----------------------------------------------------------------
  for(i=0; i<n_frm_masses; i++){
    frm_cg_arg[i]->mass = hmd_arg->frm_mass[i];
    frm_cg_arg[i]->max_num_iter = hmd_arg->max_num_iter[i];
    frm_cg_arg[i]->stop_rsd = hmd_arg->stop_rsd[i];
  }

  // Allocate memory for the boson CG arguments.
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    bsn_cg_arg = (CgArg **) smalloc(n_bsn_masses * sizeof(CgArg*));
    if(bsn_cg_arg == 0)
      ERR.Pointer(cname,fname, "bsn_cg_arg");
    VRB.Smalloc(cname,fname,
		"bsn_cg_arg",bsn_cg_arg, n_bsn_masses * sizeof(CgArg*));
    
    for(i=0; i<n_bsn_masses; i++){
      bsn_cg_arg[i] = (CgArg *) smalloc(sizeof(CgArg));
      if(bsn_cg_arg[i] == 0)
	ERR.Pointer(cname,fname, "bsn_cg_arg[i]");
      VRB.Smalloc(cname,fname,
		  "bsn_cg_arg[i]", bsn_cg_arg[i], sizeof(CgArg));
    } 
  }

  // Initialize the boson CG arguments
  //----------------------------------------------------------------
  //??? Complete this
  for(i=0; i<n_bsn_masses; i++){
    bsn_cg_arg[i]->mass = hmd_arg->bsn_mass[i];
    bsn_cg_arg[i]->max_num_iter = hmd_arg->max_num_iter[i];
    bsn_cg_arg[i]->stop_rsd = hmd_arg->stop_rsd[i];
  }

  // Allocate memory for the phi pseudo fermion field.
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    phi = (Vector **) smalloc(n_frm_masses * sizeof(Vector*));
    VRB.Smalloc(cname,fname, "phi",phi, n_frm_masses * sizeof(Vector*));
    for(i=0; i<n_frm_masses; i++){
      phi[i] = (Vector *) smalloc(f_size * sizeof(Float));
      VRB.Smalloc(cname,fname, "phi[i]", phi[i], f_size * sizeof(Float));
    }  
  }

  // Allocate memory for the frmn solution fermion fields.
  //----------------------------------------------------------------
  total_size = 0;
  for (i=0; i<n_frm_masses; i++) total_size += hmd_arg->FRatDeg[i];

  frmn = (Vector **) smalloc(total_size * sizeof(Vector*));
  VRB.Smalloc(cname,fname, "frmn",frmn, n_frm_masses * sizeof(Vector*));
  frmn_d = (Vector **)smalloc(total_size * sizeof(Vector*));
  VRB.Smalloc(cname,fname, "frmn_d",frmn_d, total_size * sizeof(Vector*));
  for (i=0; i<total_size; i++) {
    frmn[i] = (Vector *) smalloc(f_size *sizeof(Float));
    VRB.Smalloc(cname,fname, "frmn[i]", frmn[i], f_size *sizeof(Float));
    frmn_d[i] = (Vector*)smalloc(f_size * sizeof(Float));
    VRB.Smalloc(cname,fname, "frmn_d[i]", frmn_d[i], f_size * sizeof(Float));
  }

  // Allocate memory for the boson field bsn.
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    bsn = (Vector **) smalloc(n_bsn_masses * sizeof(Vector*));
    if(bsn == 0)
      ERR.Pointer(cname,fname, "bsn");
    VRB.Smalloc(cname,fname, "bsn",bsn, n_bsn_masses * sizeof(Vector*));
    for(i=0; i<n_bsn_masses; i++){
      bsn[i] = (Vector *) smalloc(f_size * sizeof(Float));
      if(bsn[i] == 0)
	ERR.Pointer(cname,fname, "bsn[i]");
      VRB.Smalloc(cname,fname, "bsn[i]", bsn[i], f_size * sizeof(Float));
    }  
  }

  // Allocate memory for the initial gauge field.
  //----------------------------------------------------------------
  gauge_field_init = (Matrix *) smalloc(g_size * sizeof(Matrix));
  if(gauge_field_init == 0)
    ERR.Pointer(cname,fname, "gauge_field_init"); 
  VRB.Smalloc(cname,fname,
	      "gauge_field_init",gauge_field_init, 
	      g_size * sizeof(Matrix));

  // Used for storing residue coefficients for staggered optimisation
  //----------------------------------------------------------------
  if (AlgLattice().Fclass() == F_CLASS_ASQTAD) {
    all_res = (Float *)smalloc(total_size * sizeof(Float));
    VRB.Smalloc(cname,fname, "all_res", all_res, total_size * sizeof(Float));
    Float *res = all_res;
    for (i=0; i<n_frm_masses; i++) {
      for (j=0; j<hmd_arg->FRatDeg[i]; j++)
	*(res++) = hmd_arg->FRatRes[i][j];
    }
  }
}

//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgHmcRHMC::~AlgHmcRHMC() {
  int i,j;
  char *fname = "~AlgHmcRHMC()" ;
  VRB.Func(cname,fname);

  // Free memory for the residue coefficients
  //----------------------------------------------------------------
  if (AlgLattice().Fclass() == F_CLASS_ASQTAD) {
    VRB.Sfree(cname,fname, "all_res",all_res);
    sfree(all_res);
  }

  // Free memory for the initial gauge field.
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "gauge_field_init",gauge_field_init);
  sfree(gauge_field_init);

  // Free memory for the boson field bsn.
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    for(i=0; i<n_bsn_masses; i++){
      VRB.Sfree(cname,fname, "bsn[i]",bsn[i]);
      sfree(bsn[i]);
    }
    VRB.Sfree(cname,fname, "bsn",bsn);
    sfree(bsn);
  }

  // Free memory for the phi (pseudo fermion) fermion field.
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    for(i=0; i<n_frm_masses; i++){
      VRB.Sfree(cname,fname, "phi[i]",phi[i]);
      sfree(phi[i]);
    }
    VRB.Sfree(cname,fname, "phi",phi);
    sfree(phi);
  }

  // Free memory for the frmn (pseudo fermion) solution fields.
  //----------------------------------------------------------------
  for(i=0; i<total_size; i++){
    VRB.Sfree(cname,fname, "frmn[i]",frmn[i]);
    sfree(frmn[i]);
    VRB.Sfree(cname,fname, "frmn_d[i]",frmn_d[i]);
    sfree(frmn_d[i]);
  }
  VRB.Sfree(cname,fname, "frmn",frmn);
  sfree(frmn);
  VRB.Sfree(cname,fname, "frmn_d",frmn_d);
  sfree(frmn_d);

  // Free memory for the boson CG arguments
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    for(i=0; i<n_bsn_masses; i++){
      VRB.Sfree(cname,fname, "bsn_cg_arg[i]",bsn_cg_arg[i]);
      sfree(bsn_cg_arg[i]);
    }
    VRB.Sfree(cname,fname, "bsn_cg_arg",bsn_cg_arg);
    sfree(bsn_cg_arg);
  }

  // Free memory for the fermion CG arguments
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    for(i=0; i<n_frm_masses; i++){
      VRB.Sfree(cname,fname, "frm_cg_arg[i]",frm_cg_arg[i]);
      sfree(frm_cg_arg[i]);
    }
    VRB.Sfree(cname,fname, "frm_cg_arg",frm_cg_arg);
    sfree(frm_cg_arg);
  }
  
}


//------------------------------------------------------------------
//
// run(): The Rational Hybrid Monte Carlo algorithm.
//
//------------------------------------------------------------------
Float AlgHmcRHMC::run(void)
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  int step;                            // Trajectory step
  Float h_init=0;                        // Initial Hamiltonian
  Float h_final=0;                       // Final Hamiltonian
  Float delta_h=0;                       // Final-Initial Hamiltonian
  int accept;
  int cg_iter=0;
  Float cg_iter_av=0;
  int cg_iter_min=0;
  int cg_iter_max=0;
  Float true_res=0;
  Float true_res_av=0;
  Float true_res_min=0;
  Float true_res_max=0;
  int cg_calls=0;
  int shift;
  int i=0, j=0;
  char *fname = "run()";
  char *md_time_str = "MD_time/step_size = ";
  FILE *fp;
  VRB.Func(cname,fname);

  Float trueMass=0;
  Float zeroPole;
  Float acceptance;                            // The acceptance probability
  Float efficiency;
 
  // Get the Lattice object
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // Set the microcanonical time step
  //----------------------------------------------------------------
  Float dt = hmd_arg->step_size;

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
  
  
  if(hmd_arg->metropolis){
    // Save initial gauge field configuration
    //--------------------------------------------------------------
    lat.CopyGaugeField(gauge_field_init);
  }


  // Heat Bath for the conjugate momenta
  //----------------------------------------------------------------
  lat.RandGaussAntiHermMatrix(mom, 1.0);
  
  // Heat Bath for the boson field bsn
  //----------------------------------------------------------------
  for(i=0; i<n_bsn_masses; i++){
    lat.RandGaussVector(frmn[0], 0.5, Ncb);
    lat.RandGaussVector(bsn[i], 0.5, Ncb);
    lat.SetPhi(phi[i], frmn[0], bsn[i], hmd_arg->bsn_mass[i]);
    lat.RandGaussVector(bsn[i], 0.5, Ncb);
    lat.FmatEvlInv(bsn[i], phi[i], bsn_cg_arg[i], CNV_FRM_NO);
  }

  // Heat Bath for the pseudo-fermions (phi)
  // Use the SI approximation since need the inverse of the action
  //----------------------------------------------------------------

  h_init = lat.GhamiltonNode() + lat.MomHamiltonNode(mom);
  for(i=0; i<n_frm_masses; i++){
    
    lat.RandGaussVector(frmn[0], 0.5, Ncb);
    h_init += lat.FhamiltonNode(frmn[0],frmn[0]);

    phi[i] -> CopyVec(frmn[0],f_size);
    phi[i] -> VecTimesEquFloat(hmd_arg->SIRatNorm[i], f_size);

    // Can only renormalise mass for staggered or asqtad cases
    if (lat.Fclass() == F_CLASS_ASQTAD || lat.Fclass() == F_CLASS_STAG) {
      trueMass = frm_cg_arg[i] -> mass;
      frm_cg_arg[i] -> mass = sqrt(trueMass*trueMass + hmd_arg->SIRatPole[i][0]/4.0);
      zeroPole = hmd_arg->SIRatPole[i][0];
      for (j=0; j<hmd_arg->SRatDeg[i]; j++) hmd_arg->SIRatPole[i][j] -= zeroPole;
    }

    cg_iter = lat.FmatEvlMInv(phi+i, frmn[0], hmd_arg->SIRatPole[i], hmd_arg->SRatDeg[i],
    			      hmd_arg->isz, frm_cg_arg[i], CNV_FRM_NO, SINGLE, hmd_arg->SIRatRes[i]);

    if (lat.Fclass() == F_CLASS_ASQTAD || lat.Fclass() == F_CLASS_STAG) {
      for (j=0; j<hmd_arg->SRatDeg[i]; j++) hmd_arg->SIRatPole[i][j] += zeroPole;
      frm_cg_arg[i] -> mass = trueMass;
    }

    cg_iter_av += cg_iter;
    if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
    if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
    true_res_av += true_res;
    if(true_res < true_res_min) true_res_min = true_res;
    if(true_res > true_res_max) true_res_max = true_res;
    cg_calls++;      

  }

  // Calculate initial boson contribution to the Hamiltonian
  //---------------------------------------------------------------
  for(i=0; i<n_bsn_masses; i++)
    h_init += lat.BhamiltonNode(bsn[i], hmd_arg->bsn_mass[i]);

  //----------------------------------------------------------------
  // Molecular Dynamics Trajectory
  //----------------------------------------------------------------

  // Perform initial UQPQ integration
  //--------------------------------------------------------------
  lat.EvolveGfield(mom, 0.25*dt/(Float)hmd_arg->sw);
  for (i=0; i<hmd_arg->sw; i++) {
    lat.EvolveMomGforce(mom, 0.5*dt/(Float)hmd_arg->sw);
    if (i < hmd_arg->sw-1) lat.EvolveGfield(mom, 0.5*dt/(Float)hmd_arg->sw);
    else lat.EvolveGfield(mom, 0.25*dt/(Float)hmd_arg->sw);
  }

  // First leap frog step has occurred. Increment MD Time clock in
  // Lattice one half time step.
  //----------------------------------------------------------------
  lat.MdTimeInc(0.5) ;
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));


  // Run through the trajectory
  //----------------------------------------------------------------

  for(step=0; step < hmd_arg->steps_per_traj; step++){

    // Evolve momenta by one step using the fermion force
    //--------------------------------------------------------------

    shift = 0;
    for(i=0; i<n_frm_masses; i++){
      // Can only renormalise mass for staggered or asqtad cases
      if (lat.Fclass() == F_CLASS_ASQTAD || lat.Fclass() == F_CLASS_STAG) {
	trueMass = frm_cg_arg[i] -> mass;
	frm_cg_arg[i] -> mass = sqrt(trueMass*trueMass + hmd_arg->FRatPole[i][0]/4.0);
	zeroPole = hmd_arg->FRatPole[i][0];
	for (j=0; j<hmd_arg->FRatDeg[i]; j++) hmd_arg->FRatPole[i][j] -= zeroPole;
      }

      for (j=0; j<hmd_arg->FRatDeg[i]; j++) frmn[j] -> VecTimesEquFloat(0.0,f_size);

      cg_iter = lat.FmatEvlMInv(frmn+shift, phi[i], hmd_arg->FRatPole[i], hmd_arg->FRatDeg[i],
				hmd_arg->isz, frm_cg_arg[i], CNV_FRM_NO, frmn_d+shift);
      shift += hmd_arg->FRatDeg[i];

      if (lat.Fclass() == F_CLASS_ASQTAD || lat.Fclass() == F_CLASS_STAG) {
	for (j=0; j<hmd_arg->FRatDeg[i]; j++) hmd_arg->FRatPole[i][j] += zeroPole;
	frm_cg_arg[i] -> mass = trueMass;
      }

      cg_iter_av += cg_iter;
      if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
      if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
      true_res_av += true_res;
      if(true_res < true_res_min) true_res_min = true_res;
      if(true_res > true_res_max) true_res_max = true_res;
      cg_calls++;      
      
      if (lat.Fclass() != F_CLASS_ASQTAD)
	lat.RHMC_EvolveMomFforce(mom, frmn, hmd_arg->FRatDeg[i], hmd_arg->FRatRes[i], 
				 hmd_arg->frm_mass[i], dt, frmn_d);
      //for(j=0; j<hmd_arg->FRatDeg[i]; j++)
      //	lat.EvolveMomFforce(mom, frmn[j], hmd_arg->frm_mass[i], hmd_arg->FRatRes[i][j]*dt);
    }

    // Only for the case of asqtad fermions do we perform this optimisation
    //--------------------------------------------------------------
    if (lat.Fclass() == F_CLASS_ASQTAD)
      lat.RHMC_EvolveMomFforce(mom, frmn, total_size, all_res, 0.0, dt, frmn_d);

    // Evolve momenta by one step using the boson force
    //--------------------------------------------------------------
    for(i=0; i<n_bsn_masses; i++)
      lat.EvolveMomFforce(mom, bsn[i], hmd_arg->bsn_mass[i], -dt);

    // Increment MD Time clock in Lattice by one half time step.
    //--------------------------------------------------------------
    lat.MdTimeInc(0.5);
    VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));


    if (step < hmd_arg->steps_per_traj-1) 
      {
	// Perform UQPQ integration (pure gauge and momenta)
	//----------------------------------------------------------------
	lat.EvolveGfield(mom, 0.5*dt/(Float)hmd_arg->sw);
	for (i=0; i<hmd_arg->sw; i++) {
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
	// Perform final UQPQ integration
	//----------------------------------------------------------------
	lat.EvolveGfield(mom, 0.25*dt/(Float)hmd_arg->sw);
	for (i=0; i<hmd_arg->sw; i++) {
	  lat.EvolveMomGforce(mom, 0.5*dt/(Float)hmd_arg->sw);
	  if (i < hmd_arg->sw-1) lat.EvolveGfield(mom, 0.5*dt/(Float)hmd_arg->sw);
	  else lat.EvolveGfield(mom, 0.25*dt/(Float)hmd_arg->sw);
	}

      }

  }

  // Reunitarize
  //----------------------------------------------------------------
  Float dev = 0.0;
  Float max_diff = 0.0;
  if(hmd_arg->reunitarize == REUNITARIZE_YES){
    lat.Reunitarize(dev, max_diff);
  }
  h_final = lat.GhamiltonNode();

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
    frmn[0] -> CopyVec(phi[i],f_size);
    frmn[0] -> VecTimesEquFloat(hmd_arg->SRatNorm[i], f_size);

    // Can only renormalise mass for staggered or asqtad cases
    if (lat.Fclass() == F_CLASS_ASQTAD || lat.Fclass() == F_CLASS_STAG) {
      trueMass = frm_cg_arg[i] -> mass;
      frm_cg_arg[i] -> mass = sqrt(trueMass*trueMass + hmd_arg->SRatPole[i][0]/4.0);
      zeroPole = hmd_arg->SRatPole[i][0];
      for (j=0; j<hmd_arg->SRatDeg[i]; j++) hmd_arg->SRatPole[i][j] -= zeroPole;
    }

    cg_iter = lat.FmatEvlMInv(frmn+shift, phi[i], hmd_arg->SRatPole[i], hmd_arg->SRatDeg[i],
    			      hmd_arg->isz, frm_cg_arg[i], CNV_FRM_NO, SINGLE, hmd_arg->SRatRes[i]);
    shift += hmd_arg->FRatDeg[i];

    if (lat.Fclass() == F_CLASS_ASQTAD || lat.Fclass() == F_CLASS_STAG) {
      for (j=0; j<hmd_arg->SRatDeg[i]; j++) hmd_arg->SRatPole[i][j] += zeroPole;
      frm_cg_arg[i] -> mass = trueMass;
    }

    cg_iter_av += cg_iter;
    if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
    if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
    true_res_av += true_res;
    if(true_res < true_res_min) true_res_min = true_res;
    if(true_res > true_res_max) true_res_max = true_res;
    cg_calls++;      

    h_final += lat.FhamiltonNode(frmn[0], frmn[0]);
  }
  
  h_final = h_final + lat.MomHamiltonNode(mom);

  // Calculate final boson contribution to the Hamiltonian
  //---------------------------------------------------------------
  for(i=0; i<n_bsn_masses; i++)
    h_final = h_final + lat.BhamiltonNode(bsn[i], 
				      hmd_arg->bsn_mass[i]);

  // Calculate Final-Initial Hamiltonian 
  //---------------------------------------------------------------
  delta_h = h_final - h_init;
  glb_sum(&delta_h);


  // Check that delta_h is the same accross all s-slices 
  // (relevant only if GJP.Snodes() != 1)
  //----------------------------------------------------------------
  if(GJP.Snodes() != 1) {
    VRB.Flow(cname,fname, "Checking Delta H across s-slices\n");
    lat.SoCheck(delta_h);
  }


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
  // accross s-slices of processors. It must be identical
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
    if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    fprintf(fp,"%d %e %d %e %e %e %d %d %e %e %e\n",
	    hmd_arg->steps_per_traj+2,
	    IFloat(delta_h), 
	    accept, 
	    IFloat(dev),
	    IFloat(max_diff),
	    IFloat(cg_iter_av),
	    cg_iter_min,
	    cg_iter_max,
	    IFloat(true_res_av),
	    IFloat(true_res_min),
	    IFloat(true_res_max));
    fclose(fp);
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
  

  // Reset Molecular Dynamics time counter
  //----------------------------------------------------------------
  lat.MdTime(0.0);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));

  return acceptance;
}

// Measure Eigenvalue bounds and regenerate approximation if necessary
void AlgHmcRHMC::dynamicalApprox()
{

  char *fname = "dynamicalApprox()";
  VRB.Func(cname,fname);

  // Get the Lattice object
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // First measure highest and lowest eigenvalues
  eig_arg->N_eig = 1;

  Vector **psi = (Vector**) smalloc(sizeof(Vector*));
  if(psi == 0) ERR.Pointer(cname,fname, "psi");
  VRB.Smalloc(cname,fname, "psi",psi,sizeof(Vector*));
  Float *lambda = (Float*) smalloc(2*sizeof(Float));
  Float *chirality = (Float*) smalloc(sizeof(Float));
  int *valid_eig = (int*) smalloc(sizeof(int));
  Float **hsum = (Float**) smalloc(sizeof(Float*));
  if (eig_arg->hsum_dir == 0)
    *hsum = (Float*) smalloc(GJP.Xnodes()*GJP.XnodeSites() * sizeof(Float));
  else if (eig_arg->hsum_dir == 1)
    *hsum = (Float*) smalloc(GJP.Ynodes()*GJP.YnodeSites() * sizeof(Float));
  else if (eig_arg->hsum_dir == 2)
    *hsum = (Float*) smalloc(GJP.Ynodes()*GJP.YnodeSites() * sizeof(Float));
  else if (eig_arg->hsum_dir == 3)
    *hsum = (Float*) smalloc(GJP.Ynodes()*GJP.YnodeSites() * sizeof(Float));

  psi[0] = (Vector *) smalloc(f_size * sizeof(Float));
  if(psi[0] == 0) ERR.Pointer(cname,fname, "psi[0]");
  
  for (int i=0; i<n_frm_masses; i++) {
    eig_arg->mass = hmd_arg->frm_mass[i];
    
    // Measure lowest e-value
    eig_arg->RitzMatOper = MATPCDAG_MATPC;
    lat.RandGaussVector(psi[0], 0.5, 1);
    lat.FeigSolv(psi, lambda, chirality, valid_eig, hsum, eig_arg, CNV_FRM_NO);
    
    // Measure highest e-value
    eig_arg->RitzMatOper = NEG_MATPCDAG_MATPC;
    lat.RandGaussVector(psi[0], 0.5, 1);
    lat.FeigSolv(psi, lambda+1, chirality, valid_eig, hsum, eig_arg, CNV_FRM_NO);
    lambda[1] *= -1;

    VRB.Flow(cname,fname, "Old Spectral Bounds are [%f,%f]\n",
	     hmd_arg->lambda_low[i],hmd_arg->lambda_high[i]);
    VRB.Flow(cname,fname, "New Spectral Bounds are [%f,%f]\n",lambda[0],lambda[1]);

    // Need to reconstruct approximations if either the spectrum leaves the interval 
    // or if it moves much inside the interval
    
    Float delta = eig_arg->Rsdlam + hmd_arg->spread;

    // Reset the required degree of approximation
    if (hmd_arg->FRatDegNew[i] != 0 && hmd_arg->SRatDegNew[i] != 0) {
      hmd_arg->FRatDeg[i] = hmd_arg->FRatDegNew[i];
      hmd_arg->SRatDeg[i] = hmd_arg->SRatDegNew[i];
    }

    if (lambda[0] < hmd_arg->lambda_low[i]  || // below interval
	lambda[1] > hmd_arg->lambda_high[i] || // above interval
	lambda[0] > hmd_arg->lambda_low[i]  * (1.0+delta) || // inside from below
	lambda[1] < hmd_arg->lambda_high[i] * (1.0-delta) ||
	hmd_arg->FRatDeg[i] != hmd_arg->FRatDegNew[i] ||
	hmd_arg->SRatDeg[i] != hmd_arg->SRatDegNew[i] ) { // inside from above
      VRB.Flow(cname,fname,"Reconstructing approximation\n");

      // Free and Allocate memory for force approximation if necessary
       if (hmd_arg->FRatDeg[i] != hmd_arg->FRatDegNew[i]) {

	for(i=0; i<total_size; i++){
	  VRB.Sfree(cname,fname, "frmn[i]",frmn[i]);
	  sfree(frmn[i]);
	  VRB.Sfree(cname,fname, "frmn_d[i]",frmn_d[i]);
	  sfree(frmn_d[i]);
	}
	VRB.Sfree(cname,fname, "frmn",frmn);
	sfree(frmn);
	VRB.Sfree(cname,fname, "frmn_d",frmn_d);
	sfree(frmn_d);

	total_size = 0;
	for (i=0; i<n_frm_masses; i++) total_size += hmd_arg->FRatDeg[i];
	
	frmn = (Vector **) smalloc(total_size * sizeof(Vector*));
	VRB.Smalloc(cname,fname, "frmn",frmn, n_frm_masses * sizeof(Vector*));
	frmn_d = (Vector **)smalloc(total_size * sizeof(Vector*));
	VRB.Smalloc(cname,fname, "frmn_d",frmn_d, total_size * sizeof(Vector*));
	for (i=0; i<total_size; i++) {
	  frmn[i] = (Vector *) smalloc(f_size *sizeof(Float));
	  VRB.Smalloc(cname,fname, "frmn[i]", frmn[i], f_size *sizeof(Float));
	  frmn_d[i] = (Vector*)smalloc(f_size * sizeof(Float));
	  VRB.Smalloc(cname,fname, "frmn_d[i]", frmn_d[i], f_size * sizeof(Float));
	}
	hmd_arg->FRatDeg[i] = hmd_arg->FRatDegNew[i];
      }

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

      AlgRemez remez(hmd_arg->lambda_low[i], hmd_arg->lambda_high[i], hmd_arg->precision);
      
      remez.generateApprox(hmd_arg->FRatDeg[i],hmd_arg->frm_power_num[i], hmd_arg->frm_power_den[i]);
      remez.getIPFE(hmd_arg->FRatRes[i], hmd_arg->FRatPole[i], hmd_arg->FRatNorm+i);
      remez.generateApprox(hmd_arg->SRatDeg[i],hmd_arg->frm_power_num[i],2*hmd_arg->frm_power_den[i]);
      remez.getIPFE(hmd_arg->SRatRes[i],hmd_arg->SRatPole[i], hmd_arg->SRatNorm+i);
      remez.getPFE(hmd_arg->SIRatRes[i],hmd_arg->SIRatPole[i], hmd_arg->SIRatNorm+i);
      
    } else {
      VRB.Flow(cname,fname,"Reconstruction not necessary\n");
    }
    
  }

  sfree(hsum[0]);
  sfree(psi[0]);
  sfree(psi);
  sfree(hsum);
  sfree(lambda);
  sfree(chirality);
  sfree(valid_eig);

}


CPS_END_NAMESPACE
