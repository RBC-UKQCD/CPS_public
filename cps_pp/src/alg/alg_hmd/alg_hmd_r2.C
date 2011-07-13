#include<config.h>
#include<stdlib.h>
#include<math.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_hmd_r.C
//
// AlgHmdR2 is derived from AlgHmd and is relevant to the Hybrid  
// Molecular Dynamics R2 algorithm. Boson fields are simulated as
// fermion fields with negative flavor number.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdio.h>
#include <alg/alg_hmd.h>
#include <comms/glb.h>
#include <util/lattice.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qcdio.h>
CPS_START_NAMESPACE




//------------------------------------------------------------------
/*!
  \param latt The lattice on which  the HMD algorithm run.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
 */
//------------------------------------------------------------------
AlgHmdR2::AlgHmdR2(Lattice& latt, 
		   CommonArg *c_arg,
		   HmdArg *arg) : 
  AlgHmd(latt, c_arg, arg) 
{
  cname = "AlgHmdR2";
  char *fname = "AlgHmdR2(L&,CommonArg*,HmdArg*)";
  VRB.Func(cname,fname);
  int i;

  if (latt.Fclass() != F_CLASS_STAG && 
      latt.Fclass() != F_CLASS_ASQTAD )
    ERR.General(cname,fname,"Cannot use R2 algorithm with non-staggered quarks\n");

  // Initialize the number of dynamical fermion masses
  //----------------------------------------------------------------
  n_frm_masses = hmd_arg->n_frm_masses;
  if(n_frm_masses != 2){
    ERR.General(cname,fname,"n_frm_masses must equal 2 for R2 algorithm\n");
  }


  // Allocate memory for the flavor time step array
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    flavor_time_step = (Float *) 
      smalloc((n_frm_masses+1)*sizeof(Float), 
	      cname, fname, "flavor_time_step");
    force_coeff = (Float *) 
      smalloc((n_frm_masses)*sizeof(Float), 
	      cname, fname, "force_coeff");
  }


  // Calculate the fermion field size.
  //----------------------------------------------------------------
  // Number of lattice sites
  f_sites = GJP.SnodeSites()*GJP.VolNodeSites() / (AlgLattice().FchkbEvl()+1);
  // Number of Vectors in a Vector array
  f_vec_count = f_sites * AlgLattice().SpinComponents();
  // Number of Floats in a Vector array
  f_size = f_vec_count * AlgLattice().Colors() * 2;

  // Allocate memory for the fermion CG arguments.
  //----------------------------------------------------------------
  frm_cg_arg = (CgArg **) 
    smalloc(sizeof(CgArg*),cname,fname,"frm_cg_arg");
  for (i=0; i<n_frm_masses; i++)
    frm_cg_arg[i] = (CgArg *) 
      smalloc(sizeof(CgArg),cname,fname,"frm_cg_arg[i]");

  light=0;
  heavy=1;
  if (hmd_arg->frm_mass[0] > hmd_arg->frm_mass[1]) {
    Float temp = hmd_arg->frm_mass[0];
    hmd_arg->frm_mass[0] = hmd_arg->frm_mass[1];
    hmd_arg->frm_mass[1] = temp;
  }

  // Initialize the fermion CG argument
  //----------------------------------------------------------------
  for (i=0; i<n_frm_masses; i++) {
    frm_cg_arg[i]->mass = hmd_arg->frm_mass[i];
    frm_cg_arg[i]->max_num_iter = hmd_arg->max_num_iter[i];
    frm_cg_arg[i]->stop_rsd = hmd_arg->stop_rsd[i];
  }

  shift = (Float*) smalloc(n_frm_masses*sizeof(Float), cname, fname, "shift");

  shift[light] = 0.0;
  shift[heavy] 
    = 4.0*(hmd_arg->frm_mass[heavy]*hmd_arg->frm_mass[heavy] - 
	   hmd_arg->frm_mass[light]*hmd_arg->frm_mass[light]);

  // Allocate memory for the phi pseudo fermion field.
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    phi = (Vector **) 
      smalloc(n_frm_masses * sizeof(Vector*), cname,fname, "phi");
    for(i=0; i<n_frm_masses; i++)
      phi[i] = (Vector *) 
	smalloc(f_size * sizeof(Float), cname,fname, "phi[i]");
  }


  // Allocate memory for solution vectors.
  //----------------------------------------------------------------
  frmn = (Vector**) smalloc(n_frm_masses*sizeof(Vector*), cname, fname, "frmn");    
  frmn_d = (Vector**) smalloc(n_frm_masses*sizeof(Vector*), cname, fname, "frmn_d");
  
  for (i=0; i<n_frm_masses; i++) {
    frmn[i] = (Vector*) smalloc(2*f_size*sizeof(Float));
    frmn_d[i] = frmn[i] + f_vec_count;
  }


}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgHmdR2::~AlgHmdR2() {
  int i;
  char *fname = "~AlgHmdR2()" ;
  VRB.Func(cname,fname);

  // Free memory for solution vectors
  //----------------------------------------------------------------
  sfree(frmn_d, cname,fname, "frmn_d");    
  for (i=0; i<n_frm_masses; i++) 
    sfree(frmn[i], cname, fname, "frmn[i]");

  sfree(frmn, cname, fname, "frmn");


  // Free memory for the phi (pseudo fermion) fermion field.
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    for(i=0; i<n_frm_masses; i++){
      sfree(phi[i], cname, fname, "phi[i]");
    }
    sfree(phi, cname, fname, "phi");
  }

  // Free memory for the shifted masses
  //----------------------------------------------------------------
  sfree(shift, cname, fname, "shift");

  // Free memory for the fermion CG arguments
  //----------------------------------------------------------------
  for (i=0; i<n_frm_masses; i++)
    sfree(frm_cg_arg[i], cname,fname, "frm_cg_arg[i]");
  sfree(frm_cg_arg, cname,fname, "frm_cg_arg");

  // Free memory for the flavor time step array
  //----------------------------------------------------------------
  if(n_frm_masses != 0) {
    sfree(flavor_time_step, cname, fname, "flavor_time_step");
    sfree(force_coeff);
  }

}


//------------------------------------------------------------------
/*!
  \post The results are written to the file specified in the common_arg
  structure.
*/
//------------------------------------------------------------------
Float AlgHmdR2::run(void)
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  int i;
  int step;                            // Trajectory step
  Float frm_time_step;
  Float flavor_coeff;
  int   cg_iter;
  Float cg_iter_av;
  int   cg_iter_min;
  int   cg_iter_max;
  Float true_res=0.;
  Float true_res_av;
  Float true_res_min;
  Float true_res_max;
  int cg_calls;
  char *fname = "run()";
  char *md_time_str = "MD_time/step_size = ";

  FILE *fp;
  VRB.Func(cname,fname);

  // Get the Lattice object
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // Set exact flavor coefficient = 1/ ExactFlavors
  //----------------------------------------------------------------
  flavor_coeff = 1.0 / lat.ExactFlavors();

  // Set the microcanonical time step
  //----------------------------------------------------------------
  Float dt = hmd_arg->step_size;

  // Initialize the flavor time step array
  //----------------------------------------------------------------
  int flavor_diff = hmd_arg->frm_flavors[0];
  flavor_time_step[0]  = -0.5 * dt * flavor_diff * flavor_coeff;
  for(i=1; i<n_frm_masses; i++){
    flavor_diff = hmd_arg->frm_flavors[i] - hmd_arg->frm_flavors[i-1];
    flavor_time_step[i]  = -0.5 * dt * flavor_diff * flavor_coeff;
    force_coeff[i] = hmd_arg->frm_flavors[i] * flavor_coeff;
  }
  flavor_diff = hmd_arg->frm_flavors[n_frm_masses-1];
  flavor_time_step[n_frm_masses] = 0.5 * dt * flavor_diff * flavor_coeff;

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


  // generate Gaussian random momentum
  //----------------------------------------------------------------
  lat.RandGaussAntiHermMatrix(mom, 1.0);
  Float mom_sum = lat.MomHamiltonNode(mom);
  glb_sum(&mom_sum);
  VRB.Flow(cname,fname,"mom_sum = %0.14e\n",mom_sum);
  Float *phi_p = (Float *)mom;


  // Evolve gauge field by dt/2
  //--------------------------------------------------------------
  lat.EvolveGfield(mom, 0.5*dt);
  lat.MdTimeInc(0.5);
  VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));


  // Run through the trajectory
  //----------------------------------------------------------------
  step = hmd_arg->steps_per_traj;
  while(step) {

    // Set the field phi for each mass.
    // First evolve gauge field by flavor_time_step.
    // Next generate a Gaussian random vector.
    // Finally calculate phi using the evolved gauge field
    // and gaussian random vector.
    //--------------------------------------------------------------
    for(i=0; i<n_frm_masses; i++){
      lat.EvolveGfield(mom, flavor_time_step[i]);
      lat.MdTimeInc(flavor_time_step[i] / dt);
      VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
      lat.RandGaussVector(frmn[0], 0.5, Ncb);
      lat.RandGaussVector(frmn[1], 0.5, Ncb);	
      lat.SetPhi(phi[i], frmn[0], frmn[1], hmd_arg->frm_mass[i]);
    }

    // Evolve gauge field to the mid-point N*dt + dt/2
    //--------------------------------------------------------------
    lat.EvolveGfield(mom, flavor_time_step[n_frm_masses]);
    lat.MdTimeInc(flavor_time_step[n_frm_masses] / dt);
    VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));

    // Evolve momenta by one step using the pure gauge force
    //--------------------------------------------------------------
    lat.EvolveMomGforce(mom, dt);

    // Evolve momenta by one step using the fermion force
    //--------------------------------------------------------------

    // First set the initial guess for the generalised
    // multi-mass solver
    frmn[light] -> FTimesV1MinusV2(1.0, phi[heavy], phi[light], f_size);
    frmn[light] -> VecTimesEquFloat(1.0/shift[heavy], f_size);
    frmn[heavy] -> CopyVec(frmn[light], f_size);

    cg_iter = lat.FmatEvlMInv(frmn, phi[light], shift, n_frm_masses, light, 
			      frm_cg_arg, CNV_FRM_NO, GENERAL, 
			      frmn_d);

    cg_iter_av = cg_iter_av + cg_iter;
    if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
    if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
    true_res_av = true_res_av + true_res;
    if(true_res < true_res_min) true_res_min = true_res;
    if(true_res > true_res_max) true_res_max = true_res;
    true_res_av = true_res_av + true_res;
    if(true_res < true_res_min) true_res_min = true_res;
    if(true_res > true_res_max) true_res_max = true_res;
    cg_calls++;

    // Now calculate the force (currently uses RHMC force)
    lat.RHMC_EvolveMomFforce(mom, frmn, 2, 0, force_coeff, 0.0, dt, frmn_d,
			     FORCE_MEASURE_NO);
    // decrease the loop counter by 1
    step--;

    // if not last step, move links forward by dt
    // otherwise, exit the loop
    if(step > 0) lat.EvolveGfield(mom, dt);

    // Increment MD Time clock in Lattice by one time step.
    lat.MdTimeInc(1.0);
    VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
  }


  //--------------------------------------------------------------
  // move Links forward half step, finish leap-
  // frog integration
  //--------------------------------------------------------------
  lat.EvolveGfield(mom, 0.5*dt);


  //----------------------------------------------------------------
  // Reunitarize
  //----------------------------------------------------------------
  Float dev = 0.0;
  Float max_diff = 0.0;
  if(hmd_arg->reunitarize == REUNITARIZE_YES){
    lat.Reunitarize(dev, max_diff);
  }

  //----------------------------------------------------------------
  // Update gauge field counter
  //----------------------------------------------------------------
  lat.GupdCntInc(1);

  //----------------------------------------------------------------
  // If GJP.Snodes() !=1  the gauge field is spread out
  // accross s-slices of processors. It must be identical
  // on each slice. Check to make sure and exit if it
  // is not identical. A case where this is relevant
  // is the DWF spread-out case.
  //----------------------------------------------------------------
  lat.GsoCheck();

  // Calculate average of monitor variables
  //---------------------------------------------------------------
  cg_iter_av = Float(cg_iter_av) / Float(cg_calls);
  true_res_av = Float(true_res_av) / Float(cg_calls);


  // Print out monitor info
  //---------------------------------------------------------------
  if(common_arg->results != 0){
    if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    Fprintf(fp,"%d %e %e %e %d %d %e %e %e\n",
	    hmd_arg->steps_per_traj,
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
  "Hmd steps = %d, dev = %e, max_diff = %e\n",
	     hmd_arg->steps_per_traj,
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

  return (Float)1.0;

}






CPS_END_NAMESPACE
