#define MINV 1

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
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
AlgHmcRHMC::AlgHmcRHMC(Lattice& latt, 
		     CommonArg *c_arg,
		     HmdArg *arg) : 
		     AlgHmd(latt, c_arg, arg) 
{
  int i,j;
  cname = "AlgHmcRHMC";
  char *fname = "AlgHmcRHMC(L&,CommonArg*,HmdArg*)";
  VRB.Func(cname,fname);
  int n_masses;

  isz = hmd_arg->isz;
  sw = hmd_arg->sw;

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

  // Calculate the fermion field size.
  //----------------------------------------------------------------
  f_size = GJP.VolNodeSites() * latt.FsiteSize() / (latt.FchkbEvl()+1);

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
    if(phi == 0)
      ERR.Pointer(cname,fname, "phi");
    VRB.Smalloc(cname,fname, "phi",phi, n_frm_masses * sizeof(Vector*));
    for(i=0; i<n_frm_masses; i++){
      phi[i] = (Vector *) smalloc(f_size * sizeof(Vector));
      if(phi[i] == 0)
	ERR.Pointer(cname,fname, "phi[i]");
      VRB.Smalloc(cname,fname, "phi[i]", phi[i], f_size * sizeof(Vector));
    }  
  }

  // Allocate memory for the boson field bsn.
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    bsn = (Vector **) smalloc(n_bsn_masses * sizeof(Vector*));
    if(bsn == 0)
      ERR.Pointer(cname,fname, "bsn");
    VRB.Smalloc(cname,fname, "bsn",bsn, n_bsn_masses * sizeof(Vector*));
    for(i=0; i<n_bsn_masses; i++){
      bsn[i] = (Vector *) smalloc(f_size * sizeof(Vector));
      if(bsn[i] == 0)
	ERR.Pointer(cname,fname, "bsn[i]");
      VRB.Smalloc(cname,fname, "bsn[i]", bsn[i], f_size * sizeof(Vector));
    }  
  }

  // Allocate memory for the initial gauge field.
  //----------------------------------------------------------------
  gauge_field_init = (Matrix *) smalloc(g_size * sizeof(Matrix));
  if(gauge_field_init == 0)
    ERR.Pointer(cname,fname, "gauge_field_init"); 
  VRB.Smalloc(cname,fname,
	      "gauge_field_init",gauge_field_init, 
	      g_size * sizeof(Float));

  // Allocate memory for the rational functions and copy over.
  //----------------------------------------------------------

  FRatDeg = (int*)smalloc(n_frm_masses * sizeof(int));
  FRatNorm = (Float*)smalloc(n_frm_masses * sizeof(Float));
  FRatPole = (Float**)smalloc(n_frm_masses * sizeof(Float*));
  FRatRes = (Float**)smalloc(n_frm_masses * sizeof(Float*));
  SRatDeg = (int*)smalloc(n_frm_masses * sizeof(int));
  SRatNorm = (Float*)smalloc(n_frm_masses * sizeof(Float));
  SRatPole = (Float**)smalloc(n_frm_masses * sizeof(Float*));
  SRatRes = (Float**)smalloc(n_frm_masses * sizeof(Float*));
  SIRatNorm = (Float*)smalloc(n_frm_masses * sizeof(Float));
  SIRatPole = (Float**)smalloc(n_frm_masses * sizeof(Float*));
  SIRatRes = (Float**)smalloc(n_frm_masses * sizeof(Float*));

  for (i=0; i<n_frm_masses; i++) {
    FRatDeg[i] = hmd_arg -> FRatDeg[i];
    FRatNorm[i] = hmd_arg -> FRatNorm[i];
    FRatRes[i] = (Float*)smalloc(FRatDeg[i] * sizeof(Float));
    FRatPole[i] = (Float*)smalloc(FRatDeg[i] * sizeof(Float));
    for (j=0; j<FRatDeg[i]; j++) {
      FRatRes[i][j] = hmd_arg -> FRatRes[i][j];
      FRatPole[i][j] = hmd_arg -> FRatPole[i][j];
    }
    SRatDeg[i] = hmd_arg -> SRatDeg[i];
    SRatNorm[i] = hmd_arg -> SRatNorm[i];
    SRatRes[i] = (Float*)smalloc(SRatDeg[i] * sizeof(Float));
    SRatPole[i] = (Float*)smalloc(SRatDeg[i] * sizeof(Float));
    SIRatNorm[i] = hmd_arg -> SIRatNorm[i];
    SIRatRes[i] = (Float*)smalloc(SRatDeg[i] * sizeof(Float));
    SIRatPole[i] = (Float*)smalloc(SRatDeg[i] * sizeof(Float));
    for (j=0; j<SRatDeg[i]; j++) {
      SRatRes[i][j] = hmd_arg -> SRatRes[i][j];
      SRatPole[i][j] = hmd_arg -> SRatPole[i][j];
      SIRatRes[i][j] = hmd_arg -> SIRatRes[i][j];
      SIRatPole[i][j] = hmd_arg -> SIRatPole[i][j];
    }
  }

  // Allocate memory for fermion/boson field arrays.
  //----------------------------------------------------------------
  n_masses = n_frm_masses;
  if(n_bsn_masses > n_frm_masses)
    n_masses = n_bsn_masses;

  if(n_masses != 0){
    frmn = (Vector ***) smalloc(n_masses * sizeof(Vector**));
    if(frmn == 0) ERR.Pointer(cname,fname, "frmn");
    VRB.Smalloc(cname,fname, "frmn",frmn, n_masses * sizeof(Vector**));

    for(i=0; i<n_masses; i++){
      // Ensure the solution vector has dimension of the largest approximation
      int rat_size = (SRatDeg[i] > FRatDeg[i]) ? SRatDeg[i] : FRatDeg[i];
      frmn[i] = (Vector**) smalloc(rat_size*sizeof(Vector*));
      if(frmn[i] == 0) ERR.Pointer(cname,fname, "frmn[i]");
      VRB.Smalloc(cname,fname, "frmn[i]", frmn[i], rat_size * sizeof(Vector*));

      for (j=0; j<rat_size; j++) {
	frmn[i][j] = (Vector*) smalloc(f_size*sizeof(Float));
	if(frmn[i][j] == 0) ERR.Pointer(cname,fname, "frmn[i][j]");
	VRB.Smalloc(cname,fname, "frmn[i][j]", frmn[i][j], f_size * sizeof(Float));
      }

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
  int n_masses;

  // Free memory for fermion solution field
  //----------------------------------------------------------------
  n_masses = n_frm_masses;
  if(n_bsn_masses > n_frm_masses)
    n_masses = n_bsn_masses;

  if(n_masses != 0){
    for(i=0; i<n_masses; i++){
      int rat_size = (SRatDeg[i] > FRatDeg[i]) ? SRatDeg[i] : FRatDeg[i];
      for (j=0; j<rat_size; j++) {
	VRB.Sfree(cname,fname, "frmn[i][j]",frmn[i][j]);
	sfree(frmn[i][j]);
      }
      VRB.Sfree(cname,fname, "frmn[i]",frmn[i]);
      sfree(frmn[i]);
    }
    VRB.Sfree(cname,fname, "frmn",frmn);
    sfree(frmn);
  }
  
  // Free memory for the rational functions
  //----------------------------------------------------------------
  for (i=0; i<n_frm_masses; i++) {
    sfree(SIRatPole[i]);
    sfree(SIRatRes[i]);
    sfree(SRatPole[i]);
    sfree(SRatRes[i]);
    sfree(FRatPole[i]);
    sfree(FRatRes[i]);
  }
  sfree(SIRatPole);
  sfree(SIRatRes);
  sfree(SIRatNorm);
  sfree(SRatPole);
  sfree(SRatRes);
  sfree(SRatNorm);
  sfree(SRatDeg);
  sfree(FRatPole);
  sfree(FRatRes);
  sfree(FRatNorm);
  sfree(FRatDeg);

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
  int i=0, j=0;
  char *fname = "run()";
  char *md_time_str = "MD_time/step_size = ";
  FILE *fp;
  VRB.Func(cname,fname);

  Float trueMass=0;
  Float acceptance;                            // The acceptance probability
 
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
    lat.RandGaussVector(frmn[i][0], 0.5, Ncb);
    lat.RandGaussVector(bsn[i], 0.5, Ncb);
    lat.SetPhi(phi[i], frmn[i][0], bsn[i], hmd_arg->bsn_mass[i]);
    lat.RandGaussVector(bsn[i], 0.5, Ncb);
    lat.FmatEvlInv(bsn[i], phi[i], bsn_cg_arg[i], CNV_FRM_NO);
  }
  

  // Heat Bath for the pseudo-fermions (phi)
  // Use the SI approximation since need the inverse of the action
  //----------------------------------------------------------------

  h_init = lat.GhamiltonNode() + lat.MomHamiltonNode(mom);
  for(i=0; i<n_frm_masses; i++){
    lat.RandGaussVector(phi[i], 0.5, Ncb);
    h_init += lat.FhamiltonNode(phi[i],phi[i]);

#ifdef MINV
    for (j=0; j<SRatDeg[i]; j++) frmn[i][j] -> VecTimesEquFloat(0.0,f_size);
    cg_iter = lat.FmatEvlMInv(frmn[i], phi[i], SIRatPole[i], SRatDeg[i],
			      isz, frm_cg_arg[i], CNV_FRM_NO);
#else
    trueMass = frm_cg_arg[i] -> mass;
    cg_iter = 0;
    for (j=0; j<SRatDeg[i]; j++) {
      frm_cg_arg[i] -> mass = sqrt(trueMass*trueMass + SIRatPole[i][j]/4.0);
      cg_iter += lat.FmatEvlInv(frmn[i][j], phi[i], frm_cg_arg[i], &true_res, CNV_FRM_NO);
    }
    frm_cg_arg[i] -> mass = trueMass;
#endif

    cg_iter_av += cg_iter;
    if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
    if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
    true_res_av += true_res;
    if(true_res < true_res_min) true_res_min = true_res;
    if(true_res > true_res_max) true_res_max = true_res;
    cg_calls++;      

    phi[i] -> VecTimesEquFloat(SIRatNorm[i], f_size);
    for (j=0; j<SRatDeg[i]; j++)
      phi[i] -> FTimesV1PlusV2(SIRatRes[i][j],frmn[i][j],phi[i],f_size);
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
  lat.EvolveGfield(mom, 0.25*dt/(Float)sw);
  for (i=0; i<sw; i++) {
    lat.EvolveMomGforce(mom, 0.5*dt/(Float)sw);
    if (i < sw-1) lat.EvolveGfield(mom, 0.5*dt/(Float)sw);
    else lat.EvolveGfield(mom, 0.25*dt/(Float)sw);
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

    for(i=0; i<n_frm_masses; i++){
#ifdef MINV
      for (j=0; j<FRatDeg[i]; j++) frmn[i][j] -> VecTimesEquFloat(0.0, f_size);
      cg_iter = lat.FmatEvlMInv(frmn[i], phi[i], FRatPole[i], FRatDeg[i],
				isz, frm_cg_arg[i], CNV_FRM_NO);
#else
      cg_iter = 0;
      trueMass = frm_cg_arg[i] -> mass;
      for (j=0; j<FRatDeg[i]; j++) {
	frm_cg_arg[i] -> mass = sqrt(trueMass*trueMass + FRatPole[i][j]/4.0);
	cg_iter += lat.FmatEvlInv(frmn[i][j],phi[i],frm_cg_arg[i],&true_res,CNV_FRM_NO);
      }
      frm_cg_arg[i] -> mass = trueMass;
#endif

      cg_iter_av += cg_iter;
      if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
      if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
      true_res_av += true_res;
      if(true_res < true_res_min) true_res_min = true_res;
      if(true_res > true_res_max) true_res_max = true_res;
      cg_calls++;      
      
      for (j=0; j<FRatDeg[i]; j++)
	lat.EvolveMomFforce(mom,frmn[i][j],hmd_arg->frm_mass[i],dt*FRatRes[i][j]);

    }

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
	lat.EvolveGfield(mom, 0.5*dt/(Float)sw);
	for (i=0; i<sw; i++) {
	  lat.EvolveMomGforce(mom, dt/(Float)sw);
	  if (i < sw-1) lat.EvolveGfield(mom, dt/(Float)sw);
	  else lat.EvolveGfield(mom, 0.5*dt/(Float)sw);
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
	lat.EvolveGfield(mom, 0.25*dt/(Float)sw);
	for (i=0; i<sw; i++) {
	  lat.EvolveMomGforce(mom, 0.5*dt/(Float)sw);
	  if (i < sw-1) lat.EvolveGfield(mom, 0.5*dt/(Float)sw);
	  else lat.EvolveGfield(mom, 0.25*dt/(Float)sw);
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


  // Calculate final fermion contribution to the Hamiltonian
  //---------------------------------------------------------------
  for (i=0; i<n_frm_masses; i++) {

#ifdef MINV
    for (j=0; j<SRatDeg[i]; j++) frmn[i][j] -> VecTimesEquFloat(0.0, f_size);
    cg_iter = lat.FmatEvlMInv(frmn[i], phi[i], SRatPole[i], SRatDeg[i],
			      isz, frm_cg_arg[i], CNV_FRM_NO);        
#else
    cg_iter = 0;
    trueMass = frm_cg_arg[i] -> mass;
    for (j=0; j<SRatDeg[i]; j++) {
      frm_cg_arg[i] -> mass = sqrt(trueMass*trueMass + SRatPole[i][j]/4.0);
      cg_iter += lat.FmatEvlInv(frmn[i][j], phi[i], frm_cg_arg[i], &true_res, CNV_FRM_NO);
    }
    frm_cg_arg[i] -> mass = trueMass;
#endif

    cg_iter_av += cg_iter;
    if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
    if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
    true_res_av += true_res;
    if(true_res < true_res_min) true_res_min = true_res;
    if(true_res > true_res_max) true_res_max = true_res;
    cg_calls++;      

    phi[i] -> VecTimesEquFloat(SRatNorm[i],f_size);
    for (j=0; j<SRatDeg[i]; j++)
      phi[i] -> FTimesV1PlusV2(SRatRes[i][j],frmn[i][j],phi[i],f_size);
    h_final += lat.FhamiltonNode(phi[i], phi[i]);
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

  // Calculate average of monitor variables
  //---------------------------------------------------------------
  cg_iter_av = Float(cg_iter_av) / Float(cg_calls);
  true_res_av = Float(true_res_av) / Float(cg_calls);

  
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

  VRB.Result(cname,fname,"Configuration number = %d\n", lat.GupdCnt());
  

  // Reset Molecular Dynamics time counter
  //----------------------------------------------------------------
  lat.MdTime(0.0);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));

  return acceptance;
}
CPS_END_NAMESPACE
