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
#include<string.h>
#include<stdio.h>
#include<math.h>
#include<alg/alg_hmd.h>
#include<alg/common_arg.h>
#include<alg/hmd_arg.h>
#include<alg/cg_arg.h>
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

#ifdef RHMC_DEBUG 
  printf("Started RHMC Constructor\n");
#endif

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
    frm_cg_arg = (CgArg **) smalloc(n_frm_masses * sizeof(int));
    if(frm_cg_arg == 0)
      ERR.Pointer(cname,fname, "frm_cg_arg");
    VRB.Smalloc(cname,fname,
		"frm_cg_arg",frm_cg_arg, n_frm_masses * sizeof(int));
    
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
  //??? Complete this
  for(i=0; i<n_frm_masses; i++){
    frm_cg_arg[i]->mass = hmd_arg->frm_mass[i];
    frm_cg_arg[i]->max_num_iter = hmd_arg->max_num_iter[i];
    frm_cg_arg[i]->stop_rsd = hmd_arg->stop_rsd[i];
  }


  // Allocate memory for the boson CG arguments.
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    bsn_cg_arg = (CgArg **) smalloc(n_bsn_masses * sizeof(int));
    if(bsn_cg_arg == 0)
      ERR.Pointer(cname,fname, "bsn_cg_arg");
    VRB.Smalloc(cname,fname,
		"bsn_cg_arg",bsn_cg_arg, n_bsn_masses * sizeof(int));
    
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
    phi = (Vector **) smalloc(n_frm_masses * sizeof(int));
    if(phi == 0)
      ERR.Pointer(cname,fname, "phi");
    VRB.Smalloc(cname,fname, "phi",phi, n_frm_masses * sizeof(int));
    for(i=0; i<n_frm_masses; i++){
      phi[i] = (Vector *) smalloc(f_size * sizeof(Float));
      if(phi[i] == 0)
	ERR.Pointer(cname,fname, "phi[i]");
      VRB.Smalloc(cname,fname, "phi[i]", phi[i], f_size * sizeof(Float));
    }  
  }


  // Allocate memory for the boson field bsn.
  //----------------------------------------------------------------
  if(n_bsn_masses != 0){
    bsn = (Vector **) smalloc(n_bsn_masses * sizeof(int));
    if(bsn == 0)
      ERR.Pointer(cname,fname, "bsn");
    VRB.Smalloc(cname,fname, "bsn",bsn, n_bsn_masses * sizeof(int));
    for(i=0; i<n_bsn_masses; i++){
      bsn[i] = (Vector *) smalloc(f_size * sizeof(Float));
      if(bsn[i] == 0)
	ERR.Pointer(cname,fname, "bsn[i]");
      VRB.Smalloc(cname,fname, "bsn[i]", bsn[i], f_size * sizeof(Float));
    }  
  }


  // Allocate memory for the initial gauge field.
  //----------------------------------------------------------------
  gauge_field_init = (Matrix *) smalloc(g_size * sizeof(Float));
  if(gauge_field_init == 0)
    ERR.Pointer(cname,fname, "gauge_field_init"); 
  VRB.Smalloc(cname,fname,
	      "gauge_field_init",gauge_field_init, 
	      g_size * sizeof(Float));


  // Allocate memory for 2 general purpose fermion/boson field 
  // arrays (frm1,frm2).
  //----------------------------------------------------------------
  n_masses = n_frm_masses;
  if(n_bsn_masses > n_frm_masses)
    n_masses = n_bsn_masses;

  if(n_masses != 0){
    frm1 = (Vector **) smalloc(n_masses * sizeof(int));
    if(frm1 == 0)
      ERR.Pointer(cname,fname, "frm1");
    VRB.Smalloc(cname,fname, "frm1",frm1, n_masses * sizeof(int));
    frm2 = (Vector **) smalloc(n_masses * sizeof(int));
    if(frm2 == 0)
      ERR.Pointer(cname,fname, "frm2");
    VRB.Smalloc(cname,fname, "frm2",frm2, n_masses * sizeof(int));
    for(i=0; i<n_masses; i++){
      frm1[i] = (Vector *) smalloc(f_size * sizeof(Float));
      if(frm1[i] == 0)
	ERR.Pointer(cname,fname, "frm1[i]");
      VRB.Smalloc(cname,fname, "frm1[i]", frm1[i], f_size * sizeof(Float));
      frm2[i] = (Vector *) smalloc(f_size * sizeof(Float));
      if(frm2[i] == 0)
	ERR.Pointer(cname,fname, "frm2[i]");
      VRB.Smalloc(cname,fname, "frm2[i]", frm2[i], f_size * sizeof(Float));
    }
  }

  // Allocate memory for the rational functions and copy over.
  //----------------------------------------------------------

  FRatDeg = (int*)smalloc(n_frm_masses * sizeof(int));
  FRatNorm = (Float*)smalloc(n_frm_masses * sizeof(Float));
  FRatConst = (Float**)smalloc(n_frm_masses * sizeof(Float*));
  FRatPole = (Float**)smalloc(n_frm_masses * sizeof(Float*));
  HBRatDeg = (int*)smalloc(n_frm_masses * sizeof(int));
  HBRatNorm = (Float*)smalloc(n_frm_masses * sizeof(Float));
  HBRatPole = (Float**)smalloc(n_frm_masses * sizeof(Float*));
  HBRatConst = (Float**)smalloc(n_frm_masses * sizeof(Float*));

  for (i=0; i<n_frm_masses; i++) {
    FRatDeg[i] = hmd_arg -> FRatDeg[i];
    FRatNorm[i] = hmd_arg -> FRatNorm[i];
    FRatConst[i] = (Float*)smalloc(FRatDeg[i] * sizeof(Float));
    FRatPole[i] = (Float*)smalloc(FRatDeg[i] * sizeof(Float));
    for (j=0; j<FRatDeg[i]; j++) {
      FRatConst[i][j] = hmd_arg -> FRatConst[i][j];
      FRatPole[i][j] = hmd_arg -> FRatPole[i][j];
    }
    HBRatDeg[i] = hmd_arg -> HBRatDeg[i];
    HBRatNorm[i] = hmd_arg -> HBRatNorm[i];
    HBRatConst[i] = (Float*)smalloc(HBRatDeg[i] * sizeof(Float));
    HBRatPole[i] = (Float*)smalloc(HBRatDeg[i] * sizeof(Float));
    for (j=0; j<HBRatDeg[i]; j++) {
      HBRatConst[i][j] = hmd_arg -> HBRatConst[i][j];
      HBRatPole[i][j] = hmd_arg -> HBRatPole[i][j];
    }
  }

#ifdef RHMC_DEBUG 
  printf("Finished RHMC Constructor\n");
#endif

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgHmcRHMC::~AlgHmcRHMC() {
  int i;
  char *fname = "~AlgHmcRHMC()" ;
  VRB.Func(cname,fname);
  int n_masses;

#ifdef RHMC_DEBUG 
  printf("Start of deconstructor\n");
  fflush(stdout);
#endif  

  // Free memory for the rational functions
  //----------------------------------------------------------------
  for (i=0; i<n_frm_masses; i++) {
    sfree(HBRatPole[i]);
    sfree(HBRatConst[i]);
    sfree(FRatPole[i]);
    sfree(FRatConst[i]);
  }
  sfree(HBRatPole);
  sfree(HBRatConst);
  sfree(HBRatNorm);
  sfree(HBRatDeg);
  sfree(FRatPole);
  sfree(FRatConst);
  sfree(FRatNorm);
  sfree(FRatDeg);

#ifdef RHMC_DEBUG 
  printf("Freed rational functions\n");
  fflush(stdout);
#endif  

  // Free memory for 2 general purpose fermion/boson field 
  // arrays (frm1,frm2).
  //----------------------------------------------------------------
  n_masses = n_frm_masses;
  if(n_bsn_masses > n_frm_masses)
    n_masses = n_bsn_masses;

  if(n_masses != 0){
    for(i=0; i<n_masses; i++){
      VRB.Sfree(cname,fname, "frm1[i]",frm1[i]);
      sfree(frm1[i]);
      VRB.Sfree(cname,fname, "frm2[i]",frm2[i]);
      sfree(frm2[i]);
    }
    VRB.Sfree(cname,fname, "frm1",frm1);
    sfree(frm1);
    VRB.Sfree(cname,fname, "frm2",frm2);
    sfree(frm2);
  }

#ifdef RHMC_DEBUG 
  printf("Freed 2 fermion fields\n");
  fflush(stdout);
#endif  

  // Free memory for the initial gauge field.
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "gauge_field_init",gauge_field_init);
  sfree(gauge_field_init);

#ifdef RHMC_DEBUG 
  printf("Freed initial gauge field\n");
  fflush(stdout);
#endif  

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

#ifdef RHMC_DEBUG 
  printf("Freed boson field\n");
  fflush(stdout);
#endif  

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

#ifdef RHMC_DEBUG 
  printf("Freed phi fermion field\n");
  fflush(stdout);
#endif  

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

#ifdef RHMC_DEBUG 
  printf("Freed boson CG\n");
  fflush(stdout);
#endif

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

#ifdef RHMC_DEBUG 
  printf("Freed fermion CG\n");
  fflush(stdout);
#endif  
  
}


//------------------------------------------------------------------
//
// run(): The Rational Hybrid Monte Carlo algorithm.
//
//------------------------------------------------------------------
void AlgHmcRHMC::run(void)
{
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
  Vector *cg_sol;
  int i=0, j=0;
  char *fname = "run()";
  char *md_time_str = "MD_time/step_size = ";
  FILE *fp;
  VRB.Func(cname,fname);

  Float trueMass=0;

  
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
    lat.RandGaussVector(frm1[i], 0.5, Ncb);
    lat.RandGaussVector(frm2[i], 0.5, Ncb);
    lat.SetPhi(phi[i], frm1[i], frm2[i], hmd_arg->bsn_mass[i]);
    lat.RandGaussVector(bsn[i], 0.5, Ncb);
    lat.FmatEvlInv(bsn[i], phi[i], bsn_cg_arg[i], CNV_FRM_NO);
  }
  
  
  // Heat Bath for the pseudo-fermions (phi)
  //----------------------------------------------------------------

#ifdef RHMC_DEBUG 
  printf("Performing heatbath for fermions\n");
#endif
  for(i=0; i<n_frm_masses; i++){
    trueMass = frm_cg_arg[i] -> mass;
    lat.RandGaussVector(frm1[i], 0.5, Ncb);
    *phi[i] *= (Float)0.0;
    for(j=HBRatDeg[i]-1; j>0; j--) {
      frm_cg_arg[i] -> mass = sqrt(trueMass*trueMass + HBRatPole[i][j]);
      //*frm2[i] *= 0.0;
      cg_iter = lat.FmatEvlInv(frm2[i], frm1[i], frm_cg_arg[i], &true_res, CNV_FRM_NO);
      *frm2[i] *= HBRatConst[i][j];
      *phi[i] += *frm2[i];
#ifdef RHMC_DEBUG 
      printf("%d mass, %d Pole, %d cg interations\n",i,j,cg_iter);fflush(stdout);
#endif
    }
    *frm1[i] *= HBRatNorm[i];
    *phi[i] += *frm1[i];
    frm_cg_arg[i] -> mass = trueMass;
  }
#ifdef RHMC_DEBUG 
  printf("Finished heatbath for fermions\n");
#endif

  //----------------------------------------------------------------
  // Molecular Dynamics Trajectory
  //----------------------------------------------------------------

  h_init = lat.GhamiltonNode() + lat.MomHamiltonNode(mom);


  // Evolve momenta by half a step using the pure gauge force
  //----------------------------------------------------------------
  lat.EvolveMomGforce(mom, 0.5*dt);


  // Evolve momenta by half a step using the fermion force
  //----------------------------------------------------------------
#ifdef RHMC_DEBUG 
  printf("Starting initial fermion force calculation\n");
#endif
  for(i=0; i<n_frm_masses; i++){
    lat.RandGaussVector(frm1[i], 0.5, Ncb);
    true_res = 1.0; // Set to non-zero so that FmatEvlInv 
                    // will return the true residual.

    trueMass = frm_cg_arg[i] -> mass;
    for(j=0; j<FRatDeg[i]; j++) {
      frm_cg_arg[i] -> mass = sqrt(trueMass*trueMass + FRatPole[i][j]);
      cg_iter = lat.FmatEvlInv(frm2[i], phi[i], frm_cg_arg[i], &true_res, CNV_FRM_NO);
      cg_iter_av += cg_iter;
      if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
      if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
      true_res_av += true_res;
      if(true_res < true_res_min) true_res_min = true_res;
      if(true_res > true_res_max) true_res_max = true_res;
      cg_calls++;      

      h_init += FRatConst[i][j] * lat.FhamiltonNode(phi[i], frm2[i]);
      *frm2[i] *= sqrt(FRatConst[i][j]);
      lat.EvolveMomFforce(mom,frm2[i],hmd_arg->frm_mass[i], 0.5*dt);
    }
    h_init += FRatNorm[i] * lat.FhamiltonNode(phi[i], phi[i]);
    frm_cg_arg[i] -> mass = trueMass;
  }
#ifdef RHMC_DEBUG 
  printf("Finished initial fermion force calculation\n");
#endif


  // Evolve momenta by half a step using the boson force
  //----------------------------------------------------------------
  for(i=0; i<n_bsn_masses; i++){
    h_init = h_init + lat.BhamiltonNode(bsn[i], 
				    hmd_arg->bsn_mass[i]);
    lat.EvolveMomFforce(mom, bsn[i], 
			hmd_arg->bsn_mass[i], -0.5*dt);
  }


  // First leap frog step has occurred. Increment MD Time clock in
  // Lattice one half time step.
  //----------------------------------------------------------------
  lat.MdTimeInc(0.5) ;
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));


  // Run through the trajectory
  //----------------------------------------------------------------
#ifdef RHMC_DEBUG 
  printf("Starting RHMD trajectory\n");
#endif
  for(step=0; step < hmd_arg->steps_per_traj; step++){


    // Evolve gauge field by one step
    //--------------------------------------------------------------
    lat.EvolveGfield(mom, dt);


    // Increment MD Time clock in Lattice by one half time step.
    //--------------------------------------------------------------
    lat.MdTimeInc(0.5);
    VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));


    // Evolve momenta by one step using the pure gauge force
    //--------------------------------------------------------------
    lat.EvolveMomGforce(mom, dt);


    // Evolve momenta by one step using the fermion force
    //--------------------------------------------------------------

#ifdef RHMC_DEBUG 
  printf("Starting %d fermion force calculation\n",step); fflush(stdout);
#endif
    for(i=0; i<n_frm_masses; i++){
      trueMass = frm_cg_arg[i] -> mass;
      for(j=0; j<FRatDeg[i]; j++) {
	frm_cg_arg[i] -> mass = sqrt(trueMass*trueMass + FRatPole[i][j]);
	cg_iter = lat.FmatEvlInv(frm2[i], phi[i], frm_cg_arg[i], &true_res, CNV_FRM_NO);
	cg_iter_av += cg_iter;
	if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
	if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
	true_res_av += true_res;
	if(true_res < true_res_min) true_res_min = true_res;
	if(true_res > true_res_max) true_res_max = true_res;
	cg_calls++;      
	
	*frm2[i] *= sqrt(FRatConst[i][j]);
	lat.EvolveMomFforce(mom,frm2[i],hmd_arg->frm_mass[i], dt);	
      }      
      frm_cg_arg[i] -> mass = trueMass;
    }
#ifdef RHMC_DEBUG 
    printf("Finished %d fermion force calculation\n",step); fflush(stdout);
#endif

    // Evolve momenta by one step using the boson force
    //--------------------------------------------------------------
    for(i=0; i<n_bsn_masses; i++){
      lat.EvolveMomFforce(mom, bsn[i], 
			  hmd_arg->bsn_mass[i], -dt);
    }

    // Increment MD Time clock in Lattice by one half time step.
    //--------------------------------------------------------------
    lat.MdTimeInc(0.5);
    VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));

  }
#ifdef RHMC_DEBUG 
  printf("Finished RHMD trajectory\n"); fflush(stdout);
#endif


  // Evolve gauge field by one last step
  //----------------------------------------------------------------
  lat.EvolveGfield(mom, dt);


  // Reunitarize
  //----------------------------------------------------------------
  Float dev = 0.0;
  Float max_diff = 0.0;
  if(hmd_arg->reunitarize == REUNITARIZE_YES){
    lat.Reunitarize(dev, max_diff);
  }
  h_final = lat.GhamiltonNode();


  // Increment MD Time clock in Lattice by one half time step.
  //----------------------------------------------------------------
  lat.MdTimeInc(0.5);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));


  // Evolve momenta by a last half step using the pure gauge force
  //----------------------------------------------------------------
  lat.EvolveMomGforce(mom, 0.5*dt);


#ifdef RHMC_DEBUG 
  printf("Starting final fermion force calculation\n",step);
#endif
  // Evolve momenta by a last half step using the fermion force
  //----------------------------------------------------------------
  for(i=0; i<n_frm_masses; i++){
      trueMass = frm_cg_arg[i] -> mass;
      for(j=0; j<FRatDeg[i]; j++) {
	frm_cg_arg[i] -> mass = sqrt(trueMass*trueMass + FRatPole[i][j]);
	cg_iter = lat.FmatEvlInv(frm2[i], phi[i], frm_cg_arg[i], &true_res, CNV_FRM_NO);
	cg_iter_av += cg_iter;
	if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
	if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
	true_res_av += true_res;
	if(true_res < true_res_min) true_res_min = true_res;
	if(true_res > true_res_max) true_res_max = true_res;
	cg_calls++;      
	
	h_final += FRatConst[i][j] * lat.FhamiltonNode(phi[i], frm2[i]);
	*frm2[i] *= sqrt(FRatConst[i][j]);
	lat.EvolveMomFforce(mom,frm2[i],hmd_arg->frm_mass[i], 0.5*dt);	
      }
      h_final += FRatNorm[i] * lat.FhamiltonNode(phi[i], phi[i]);
      frm_cg_arg[i] -> mass = trueMass;
  }
#ifdef RHMC_DEBUG 
  printf("Finished final fermion force calculation\n",step);
#endif

  // Evolve momenta by a last half step using the boson force
  //---------------------------------------------------------------
  for(i=0; i<n_bsn_masses; i++){
    h_final = h_final + lat.BhamiltonNode(bsn[i], 
				      hmd_arg->bsn_mass[i]);
    lat.EvolveMomFforce(mom, bsn[i], 
			hmd_arg->bsn_mass[i], -0.5*dt);
  }
  h_final = h_final + lat.MomHamiltonNode(mom);


  // Calculate Final-Initial Hemiltonian 
  //---------------------------------------------------------------
  delta_h = h_final - h_init;
  printf("h0 = %e, h1 = %e, dh = %e\n",h_init,h_final,delta_h);
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
    accept = lat.MetropolisAccept(delta_h);
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
	     hmd_arg->steps_per_traj+2,
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


}
CPS_END_NAMESPACE
