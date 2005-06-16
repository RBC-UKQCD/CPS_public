#include<config.h>
#include<stdlib.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Definitions of the AlgHmdR methods.

  $Id: alg_hmd_r.C,v 1.18 2005-06-16 07:18:55 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-06-16 07:18:55 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_hmd/alg_hmd_r.C,v 1.18 2005-06-16 07:18:55 chulwoo Exp $
//  $Id: alg_hmd_r.C,v 1.18 2005-06-16 07:18:55 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_hmd_r.C,v $
//  $Revision: 1.18 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_hmd/alg_hmd_r.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_hmd_r.C
//
// AlgHmdR is derived from AlgHmd and is relevant to the Hybrid  
// Molecular Dynamics R algorithm. Boson fields are simulated as
// fermion fields with negative flavor number.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/checksum.h>
#include <util/qcdio.h>
#include <alg/alg_hmd.h>
#include <comms/glb.h>
#include <util/lattice.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
#include <qcdocos/scu_checksum.h>
#endif
CPS_START_NAMESPACE


//------------------------------------------------------------------
/*!
  \param latt The lattice on which  the HMD algorithm run.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
 */
//------------------------------------------------------------------
AlgHmdR::AlgHmdR(Lattice& latt, 
		     CommonArg *c_arg,
		     HmdArg *arg) : 
		     AlgHmd(latt, c_arg, arg) 
{
  cname = "AlgHmdR";
  char *fname = "AlgHmdR(L&,CommonArg*,HmdArg*)";
  VRB.Func(cname,fname);
  int i;

  // Initialize the number of dynamical fermion masses
  //----------------------------------------------------------------
  n_frm_masses = hmd_arg->n_frm_masses;
  if(n_frm_masses > MAX_HMD_MASSES){
    ERR.General(cname,fname,
    "hmd_arg->n_frm_masses = %d is larger than MAX_HMD_MASSES = %d\n",
     n_frm_masses, MAX_HMD_MASSES);
  }


  // Allocate memory for the flavor time step array
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    flavor_time_step = 
      (Float *) smalloc( (n_frm_masses+1) * sizeof(Float));
    if(flavor_time_step == 0)
      ERR.Pointer(cname,fname, "flavor_time_step");
    VRB.Smalloc(cname,fname, 
		"flavor_time_step",
		flavor_time_step, 
		(n_frm_masses+1) * sizeof(Float));
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


  // Allocate memory for 2 general purpose fermion fields (frm1,frm2).
  //----------------------------------------------------------------
  frm1 = (Vector *) smalloc(f_size * sizeof(Float));
  if(frm1 == 0)
    ERR.Pointer(cname,fname, "frm1");
  VRB.Smalloc(cname,fname, "frm1", frm1, f_size * sizeof(Float));
#if 1
  frm2 = (Vector *) smalloc(f_size * sizeof(Float));
  if(frm2 == 0)
    ERR.Pointer(cname,fname, "frm2");
  VRB.Smalloc(cname,fname, "frm2", frm2, f_size * sizeof(Float));
#else
  frm2 = &(frm1[GJP.VolNodeSites()/2]);
#endif


}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgHmdR::~AlgHmdR() {
  int i;
  char *fname = "~AlgHmdR()" ;
  VRB.Func(cname,fname);

  // Free memory for 2 general purpose fermion fields (frm1,frm2).
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "frm1",frm1);
  sfree(frm1);
#if 1
  VRB.Sfree(cname,fname, "frm2",frm2);
  sfree(frm2);
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

  // Free memory for the flavor time step array
  //----------------------------------------------------------------
  if(n_frm_masses != 0){
    VRB.Sfree(cname,fname, "flavor_time_step",flavor_time_step);
    sfree(flavor_time_step);
  }

}


//------------------------------------------------------------------
/*!
  \post The following results are written to the file specified in
  the CommonArg structure:
  -# The number of molecular dynamics steps
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
Float AlgHmdR::run(void)
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
  Float true_res;
  Float true_res_av;
  Float true_res_min;
  Float true_res_max;
  int cg_calls;
  char *fname = "run()";
  char *md_time_str = "MD_time/step_size = ";
  FILE *fp;
  VRB.Func(cname,fname);
  unsigned long step_cnt = 0;
  CSM.SaveComment(step_cnt);

#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
  if(!ScuChecksum::ChecksumsOn())
  ScuChecksum::Initialise(true,true);
#endif
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
    CSM.SaveComment(++step_cnt);

    // Set the field phi for each mass.
    // First evolve gauge field by flavor_time_step.
    // Next generate a Gaussian random vector.
    // Finally calculate phi using the evolved gauge field
    // and gaussian random vector.
    //--------------------------------------------------------------
    for(i=0; i<n_frm_masses; i++){
//      if (flavor_time_step[i] > 0.0) 
      lat.EvolveGfield(mom, flavor_time_step[i]);
      lat.MdTimeInc(flavor_time_step[i] / dt);
      VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
      lat.RandGaussVector(frm1, 0.5, Ncb);
      Float phi_sum = (frm1)->NormSqGlbSum(GJP.VolNodeSites()*6/2);
      VRB.Flow(cname,fname,"frm1_sum = %0.14e\n",phi_sum);
      phi_p = (Float *)frm1;

      lat.RandGaussVector(frm2, 0.5, Ncb);
      phi_sum = (frm2)->NormSqGlbSum(GJP.VolNodeSites()/2*6);
      VRB.Flow(cname,fname,"frm2_sum = %0.14e\n",phi_sum);
	
      lat.SetPhi(phi[i], frm1, frm2, hmd_arg->frm_mass[i]);
      phi_sum = (phi[i])->NormSqGlbSum(GJP.VolNodeSites()/2*6);
      phi_p = (Float *)phi[i];
      VRB.Flow(cname,fname,"phi_sum = %0.14e\n",phi_sum);
    }

    // Evolve gauge field to the mid-point N*dt + dt/2
    //--------------------------------------------------------------
    lat.EvolveGfield(mom, flavor_time_step[n_frm_masses]);
    lat.MdTimeInc(flavor_time_step[n_frm_masses] / dt);
    VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));


    // Evolve momenta by one step using the pure gauge force
    //--------------------------------------------------------------
    lat.EvolveMomGforce(mom, dt);
  Float mom_sum = lat.MomHamiltonNode(mom);
  glb_sum(&mom_sum);
  VRB.Flow(cname,fname,"mom_sum = %0.14e\n",mom_sum);
  phi_p = (Float *)mom;

    // Evolve momenta by one step using the fermion force
    //--------------------------------------------------------------
    for(i=0; i<n_frm_masses; i++){
	lat.RandGaussVector(frm1, 0.5, Ncb);
//      lat.Fconvert(frm1,STAG,CANONICAL);
      Float phi_sum = (frm1)->NormSqGlbSum(GJP.VolNodeSites()*6/2);
      VRB.Flow(cname,fname,"frm1_sum = %0.14e\n",phi_sum);
	cg_iter = 
	lat.FmatEvlInv(frm1, phi[i], 
		       frm_cg_arg[i], &true_res, CNV_FRM_NO);
      phi_sum = (phi[i])->NormSqGlbSum(GJP.VolNodeSites()*6/2);
      VRB.Flow(cname,fname,"phi_sum = %0.14e\n",phi_sum);
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
      frm_time_step = 
	hmd_arg->frm_flavors[i] * flavor_coeff * dt ;

      lat.EvolveMomFforce(mom, frm1, 
			  hmd_arg->frm_mass[i], 
			  frm_time_step);
  Float mom_sum = lat.MomHamiltonNode(mom);
  glb_sum(&mom_sum);
  VRB.Flow(cname,fname,"mom_sum = %0.14e\n",mom_sum);

    }

    //--------------------------------------------------------------
    // decrease the loop counter by 1
    //--------------------------------------------------------------
    step--;


    //--------------------------------------------------------------
    // if not last step, move links forward by dt
    // otherwise, exit the loop
    //--------------------------------------------------------------
    if(step > 0) {
      lat.EvolveGfield(mom, dt);
    }


    // Increment MD Time clock in Lattice by one time step.
    //--------------------------------------------------------------
    lat.MdTimeInc(1.0);
    VRB.Flow(cname,fname,"%s%f\n", md_time_str, IFloat(lat.MdTime()));
  }

  CSM.SaveComment(++step_cnt);

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
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
  if ( ! ScuChecksum::CsumSwap() )
    ERR.Hardware(cname,fname, "SCU Checksum mismatch\n");
#endif

  return (Float)1.0;

}






CPS_END_NAMESPACE
