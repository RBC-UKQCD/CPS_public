#include<config.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Definitions of the AlgHmcPhi methods.

  $Id: alg_hmc_phi.C,v 1.24 2006-12-14 17:53:37 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006-12-14 17:53:37 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_hmd/alg_hmc_phi.C,v 1.24 2006-12-14 17:53:37 chulwoo Exp $
//  $Id: alg_hmc_phi.C,v 1.24 2006-12-14 17:53:37 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_hmc_phi.C,v $
//  $Revision: 1.24 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_hmd/alg_hmc_phi.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_hmc_phi.C
//
// AlgHmcPhi is derived from Alg and is relevant to the phi
// Hybrid Monte Carlo Algorithm.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/checksum.h>
#include <util/qcdio.h>
#include <stdlib.h>
#include <alg/alg_hmd.h>
#include <util/lattice.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/glb.h>
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
#include <qcdocos/scu_checksum.h>
#endif
CPS_START_NAMESPACE


//------------------------------------------------------------------
/*!
  \param latt The lattice on which the HMC algorithm runs.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
 */
//------------------------------------------------------------------
AlgHmcPhi::AlgHmcPhi(Lattice& latt, 
		     CommonArg *c_arg,
		     HmdArg *arg) : 
		     AlgHmd(latt, c_arg, arg) 
{
  int i;
  cname = "AlgHmcPhi";
  char *fname = "AlgHmcPhi(L&,CommonArg*,HmdArg*)";
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
      bzero((char *)phi[i],f_size*sizeof(Float));
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


  // Allocate memory for current and previous CG solution pointers.
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    cg_sol_prev = (Vector **) smalloc(n_frm_masses * sizeof(int));
    if(cg_sol_prev == 0)
      ERR.Pointer(cname,fname, "cg_sol_prev");
    VRB.Smalloc(cname,fname, "cg_sol_prev",cg_sol_prev, 
		n_frm_masses * sizeof(int));
    cg_sol_cur = (Vector **) smalloc(n_frm_masses * sizeof(int));
    if(cg_sol_cur == 0)
      ERR.Pointer(cname,fname, "cg_sol_cur");
    VRB.Smalloc(cname,fname, "cg_sol_cur",cg_sol_cur, 
		n_frm_masses * sizeof(int));
  }


}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgHmcPhi::~AlgHmcPhi() {
  int i;
  char *fname = "~AlgHmcPhi()" ;
  VRB.Func(cname,fname);
  int n_masses;


  // Free memory for current and previous CG solution pointers.
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    sfree(cname,fname, "cg_sol_prev",cg_sol_prev);
//    sfree(cg_sol_prev);
    sfree(cname,fname, "cg_sol_cur",cg_sol_cur);
//    sfree(cg_sol_cur);
  }


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
    VRB.FuncEnd(cname,fname);
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
  .
*/
//------------------------------------------------------------------
Float AlgHmcPhi::run(void)
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  int step;                            // Trajectory step
  Float h_init;                        // Initial Hamiltonian
  Float h_final;                       // Final Hamiltonian
  Float delta_h;                       // Final-Initial Hamiltonian
  int accept;
  int   cg_iter;
  Float cg_iter_av;
  int   cg_iter_min;
  int   cg_iter_max;
  Float true_res;
  Float true_res_av;
  Float true_res_min;
  Float true_res_max;
  int cg_calls;
  Vector *cg_sol;
  int i;
  char *fname = "run()";
  char *md_time_str = "MD_time/step_size = ";
  FILE *fp;
  VRB.Func(cname,fname);
  unsigned long step_cnt = 0;
  CSM.SaveComment(step_cnt);
  Float acceptance;

#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
  if(!ScuChecksum::ChecksumsOn())
  ScuChecksum::Initialise(true,true);
#endif
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
  for(i=0; i<n_frm_masses; i++){
    lat.RandGaussVector(frm1[i], 0.5, Ncb);
    lat.RandGaussVector(frm2[i], 0.5, Ncb);
    lat.SetPhi(phi[i], frm1[i], frm2[i], hmd_arg->frm_mass[i]);
  }


  //----------------------------------------------------------------
  // Molecular Dynamics Trajectory
  //----------------------------------------------------------------

  h_init = lat.GhamiltonNode() + lat.MomHamiltonNode(mom);


  // Evolve momenta by half a step using the pure gauge force
  //----------------------------------------------------------------
  lat.EvolveMomGforce(mom, 0.5*dt);

  // Evolve momenta by half a step using the fermion force
  //----------------------------------------------------------------
  for(i=0; i<n_frm_masses; i++){
    lat.RandGaussVector(frm1[i], 0.5, Ncb);
    true_res     = 1.0; // Set to non-zero so that FmatEvlInv 
                        // will return the true residual.
    IFloat *tmp = (IFloat *)frm1[i];
    cg_iter = 
      lat.FmatEvlInv(frm1[i], phi[i], frm_cg_arg[i], &true_res, CNV_FRM_NO);
    tmp = (IFloat *)phi[i];
//    exit(32);
    cg_iter_av = cg_iter_av + cg_iter;
    if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
    if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
    true_res_av = true_res_av + true_res;
    if(true_res < true_res_min) true_res_min = true_res;
    if(true_res > true_res_max) true_res_max = true_res;
    cg_calls++;
    frm2[i]->CopyVec(frm1[i], f_size);
    cg_sol_prev[i] = frm2[i];
    cg_sol_cur[i] = frm1[i];
    h_init = h_init + lat.FhamiltonNode(phi[i], frm1[i]);
    lat.EvolveMomFforce(mom, frm1[i], 
			hmd_arg->frm_mass[i], 0.5*dt);
  }


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
  for(step=0; step < hmd_arg->steps_per_traj; step++){
    CSM.SaveComment(++step_cnt);

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
    for(i=0; i<n_frm_masses; i++){
      cg_sol_prev[i]->FTimesV1MinusV2(2.0, cg_sol_cur[i], 
				      cg_sol_prev[i], f_size);
      cg_sol = cg_sol_prev[i];
      cg_iter = 
	lat.FmatEvlInv(cg_sol, phi[i], frm_cg_arg[i], &true_res, CNV_FRM_NO);
      cg_iter_av = cg_iter_av + cg_iter;
      if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
      if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
      true_res_av = true_res_av + true_res;
      if(true_res < true_res_min) true_res_min = true_res;
      if(true_res > true_res_max) true_res_max = true_res;
      cg_calls++;
      cg_sol_prev[i] = cg_sol_cur[i];
      cg_sol_cur[i] = cg_sol;
      lat.EvolveMomFforce(mom, cg_sol, 
			  hmd_arg->frm_mass[i], dt);
    }


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

  CSM.SaveComment(++step_cnt);
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


  // Evolve momenta by a last half step using the fermion force
  //----------------------------------------------------------------
  for(i=0; i<n_frm_masses; i++){
    cg_sol_prev[i]->FTimesV1MinusV2(2.0, cg_sol_cur[i], 
				    cg_sol_prev[i], f_size);
    cg_sol = cg_sol_prev[i];
    cg_iter = 
      lat.FmatEvlInv(cg_sol, phi[i], frm_cg_arg[i], &true_res, CNV_FRM_NO);
    cg_iter_av = cg_iter_av + cg_iter;
    if(cg_iter < cg_iter_min) cg_iter_min = cg_iter;
    if(cg_iter > cg_iter_max) cg_iter_max = cg_iter;
    true_res_av = true_res_av + true_res;
    if(true_res < true_res_min) true_res_min = true_res;
    if(true_res > true_res_max) true_res_max = true_res;
    cg_calls++;
    h_final = h_final + lat.FhamiltonNode(phi[i], cg_sol);
    lat.EvolveMomFforce(mom, cg_sol, 
			hmd_arg->frm_mass[i], 0.5*dt);
  }


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
  glb_sum(&delta_h);
  VRB.Flow(cname,fname, "delta_h=%e\n",delta_h);


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
    VRB.Flow(cname,fname,"Results ptr is %p\n",common_arg->results);
    if (common_arg->results == NULL) printf("FUCK\n");fflush(stdout);
    if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    Fprintf(fp,"%d %.16e %d %.16e %.16e %.16e %d %d %.16e %.16e %.16e\n",
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
    Fclose(fp);
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

#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
  if ( ! ScuChecksum::CsumSwap() )
    ERR.Hardware(cname,fname, "SCU Checksum mismatch\n");
#endif
  return acceptance;
}

CPS_END_NAMESPACE
