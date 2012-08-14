#include<config.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
/*!\file
  \brief Definitions of the AlgHmc methods.

*/

//------------------------------------------------------------------
//
// alg_hmc.C
//
// AlgHmc is an abstract implementation of the Hybrid Monte Carlo
// Algorithm.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<math.h>
#include<alg/alg_hmc.h>
#include<alg/alg_meas.h>
#include<comms/glb.h>
#include <util/checksum.h>
#include<util/data_shift.h>
#include<util/error.h>
#include<util/gjp.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/smalloc.h>
#include<util/qcdio.h>
#include<util/wilson.h>

#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
#include <qcdocos/scu_checksum.h>
#endif
CPS_START_NAMESPACE

//#if (TARGET==QCDOC) || (TARGET==BGP)
#if (TARGET==QCDOC) 
static const int SHIFT_X = 1;
static const int SHIFT_Y = 1;
static const int SHIFT_Z = 1;
#else
static const int SHIFT_X = 0;
static const int SHIFT_Y = 0;
static const int SHIFT_Z = 0;
#endif

//------------------------------------------------------------------
/*!
  \param Integrator The integrator which defines the HMC type.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
*/
//------------------------------------------------------------------

AlgHmc::AlgHmc(AlgIntAB &Integrator, CommonArg &c_arg, HmcArg &arg)
{
  cname = "AlgHmc";
  char *fname = "AlgHmc(AlgIntAB*,CommonArg*,HmcArg*)";
  VRB.Func(cname,fname);

  integrator = &Integrator;
  hmc_arg = &arg;
  common_arg = &c_arg;

  {
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

    g_size = GJP.VolNodeSites() * lat.GsiteSize();

    //!< Allocate memory for the initial gauge field.
    gauge_field_init = (Matrix *) smalloc(g_size * sizeof(Float), 
					  "gauge_field_init",fname,cname);
    
    if (hmc_arg->reverse == REVERSE_YES) {
      
      //!< Allocate memory for the final gauge field.
      gauge_field_final = (Matrix *) smalloc(g_size * sizeof(Float), 
					     "gauge_field_final",fname,cname);
      
    }

    LatticeFactory::Destroy();
  }

} 

//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgHmc::~AlgHmc() {

//  int i,j;
  char *fname = "~AlgHmc()" ;
  VRB.Func(cname,fname);

  // Free memory for the initial gauge field.
  sfree(gauge_field_init, "gauge_field_init",fname,cname);

  if (hmc_arg->reverse == REVERSE_YES) {
    //!< Free memory for the final gauge field.
    sfree(gauge_field_final, "gauge_field_final",fname,cname);

  }

}


//------------------------------------------------------------------
//
// run(): The Hybrid Monte Carlo algorithm.
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
Float AlgHmc::run(void)
{
  char *fname = "run()";

#if TARGET==cpsMPI
  using MPISCU::fprintf;
#endif

  int accept;

  FILE *fp;
  VRB.Func(cname,fname);

  Float dev = 0.0;
  Float max_diff = 0.0;

  Float acceptance;                            // The acceptance probability
  Float efficiency = 0.0;

#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
  if(!ScuChecksum::ChecksumsOn())
  ScuChecksum::Initialise(true,true);
#endif
 
  // Set the microcanonical time step
//  Float dt = hmc_arg->step_size;

  //!< Save initial lattice and rngs (if necessary)
  int Ntests = saveInitialState();
  VRB.Result(cname,fname,"shifts=%d %d %d\n",SHIFT_X,SHIFT_Y,SHIFT_Z);

  // Try attempt_limit times to generate the same final gauge config 
  // consecutively
  for (int attempt = 0; attempt < hmc_arg->reproduce_attempt_limit; attempt++) {

    for (int test = 0; test < Ntests; test++) {
      VRB.Result(cname,fname,"Running test %d of %d, attempt %d of %d\n", 
		 test+1, Ntests, attempt+1, hmc_arg->reproduce_attempt_limit);

      //!< Must initialise the integrator
      integrator->init();

      //!< Restore state if necessary
      if ( !(test == 0 && attempt ==0) ){
        restoreInitialState();
        shiftStates(SHIFT_X,SHIFT_Y,SHIFT_Z,0);
        GDS.SetOrigin(SHIFT_X,SHIFT_Y,SHIFT_Z,0);
      }

      //!< Evaluate the heatbath
      integrator->heatbath();

      //!< Calculate initial Hamiltonian
      wilson_set_sloppy( false);
      h_init = integrator->energy();
//      Float total_h_init =h_init;
//      glb_sum(&total_h_init);

      // Molecular Dynamics Trajectory
      if(hmc_arg->wfm_md_sloppy) wilson_set_sloppy(true);
      integrator->evolve(hmc_arg->step_size, hmc_arg->steps_per_traj);
      wilson_set_sloppy(false);

      // Reunitarize
      if(hmc_arg->reunitarize == REUNITARIZE_YES){
	Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	lat.Reunitarize(dev, max_diff);
	LatticeFactory::Destroy();
      }

#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
      printf("SCU checksum test\n");
  if ( ! ScuChecksum::CsumSwap() )
    ERR.Hardware(cname,fname, "SCU Checksum mismatch\n");
#endif

      //!< Calculate final Hamiltonian
      h_final = integrator->energy();
//      Float total_h_final =h_final;
//      glb_sum(&total_h_final);

      // Calculate Final-Initial Hamiltonian 
      delta_h = h_final - h_init;
      glb_sum(&delta_h);
//      VRB.Result(cname,fname,"h_init=%0.14e h_final=%0.14e delta_h=%0.14e \n",
//        total_h_init,total_h_final,delta_h);

      // Check that delta_h is the same across all s-slices 
      // (relevant only if GJP.Snodes() != 1)
      if(GJP.Snodes() != 1) {
	Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	VRB.Flow(cname,fname, "Checking Delta H across s-slices\n");
	lat.SoCheck(delta_h);
	LatticeFactory::Destroy();
      }

      if ( !(test == 0 && attempt ==0) ){
        shiftStates(-SHIFT_X,-SHIFT_Y,-SHIFT_Z,0);
        GDS.SetOrigin(0,0,0,0);
      }

      if (hmc_arg->reverse == REVERSE_YES) {
	saveFinalState();

	integrator->reverse();
      if(hmc_arg->wfm_md_sloppy) wilson_set_sloppy(true);
	integrator->evolve(hmc_arg->step_size, hmc_arg->steps_per_traj);
      wilson_set_sloppy(false);

#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
  printf("SCU checksum test\n");
  if ( ! ScuChecksum::CsumSwap() )
    ERR.Hardware(cname,fname, "SCU Checksum mismatch\n");
#endif

	h_delta = h_final - integrator->energy();
	glb_sum(&h_delta);

	reverseTest();
	restoreFinalState();
      }

      Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
      checksum[test] = global_checksum((Float*)lat.GaugeField(),g_size);
      LatticeFactory::Destroy();

    } // end test loop
  
    if (reproduceTest(attempt)) break;

  } // end attempt loop

  {
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    
	VRB.Result(cname,fname,"hmc_arg->metropolis=%d\n",hmc_arg->metropolis);
    //!< Metropolis step
    if(hmc_arg->metropolis == METROPOLIS_YES){
      accept = lat.MetropolisAccept(delta_h,&acceptance);
      if( !(accept) ){
	// Trajectory rejected
	lat.GaugeField(gauge_field_init);
	VRB.Result(cname,fname,"Metropolis step -> Rejected\n");
      }
      else {
	// Trajectory accepted. 
	// Increment the Gauge Update counter in Lattice.
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

    //!< If GJP.Snodes() !=1  the gauge field is spread out
    //!< across s-slices of processors. It must be identical
    //!< on each slice. Check to make sure and exit if it
    //!< is not identical. A case where this is relevant
    //!< is the DWF spread-out case.
    if(GJP.Snodes() != 1) {
      VRB.Flow(cname,fname, "Checking gauge field across s-slices\n");
      lat.GsoCheck();
    }

    config_no = lat.GupdCnt();

    LatticeFactory::Destroy();
  }

  CgStats cg_stats;
  integrator->cost(&cg_stats);

  //!< Calculate average of monitor variables
  if (cg_stats.cg_calls > 0)
    efficiency = Float(acceptance) / Float(cg_stats.cg_iter_av);

  //!< Print out monitor info
  if(common_arg->results != 0){
    if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    Fprintf(fp,"%d %e %e %d %e %e %e %e %d %d %e %e %e\n",
	    hmc_arg->steps_per_traj,
	    IFloat(delta_h), 
	    acceptance,
	    accept, 
	    IFloat(dev),
	    IFloat(max_diff),
	    IFloat(cg_stats.cg_iter_total),
	    IFloat(cg_stats.cg_iter_av),
	    cg_stats.cg_iter_min,
	    cg_stats.cg_iter_max,
	    IFloat(cg_stats.true_rsd_av),
	    IFloat(cg_stats.true_rsd_min),
	    IFloat(cg_stats.true_rsd_max));
    Fclose(fp);
  }

  VRB.Result(cname,fname,
	     "Hmc steps = %d, Delta_hamilton = %.20e, accept = %d, dev = %e, max_diff = %e\n",
	     hmc_arg->steps_per_traj,
	     IFloat(delta_h), 
	     accept,
	     IFloat(dev),
	     IFloat(max_diff));
  VRB.Result(cname,fname,
	     "CG iterations: average = %e, min = %d, max = %d\n",
	     IFloat(cg_stats.cg_iter_av),
	     cg_stats.cg_iter_min,
	     cg_stats.cg_iter_max);
  VRB.Result(cname,fname,
	     "True Residual / |source|: average = %e, min = %e, max = %e\n",
	     IFloat(cg_stats.true_rsd_av),
	     IFloat(cg_stats.true_rsd_min),
	     IFloat(cg_stats.true_rsd_max));

  VRB.Result(cname,fname,
	     "Efficiency (acceptance per cg iteration) = %e\n",
	     efficiency);

  VRB.Result(cname,fname,"Configuration number = %d\n", config_no);


  return acceptance;
}


int AlgHmc::saveInitialState() {

  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

  //!< Save initial gauge field configuration
  if(hmc_arg->metropolis == METROPOLIS_YES ||
     hmc_arg->reproduce == REPRODUCE_YES ||
     hmc_arg->reverse == REVERSE_YES){
    lat.CopyGaugeField(gauge_field_init);
  }

  int Ntests;

  if (hmc_arg->reproduce == REPRODUCE_YES) {
    //!< Save initial rng state
//    LRG.GetStates(rng5d_init, FIVE_D);
//    LRG.GetStates(rng4d_init, FOUR_D);
      lrg_state.GetStates();
    
    if (hmc_arg->reproduce_attempt_limit < 1 ||
	hmc_arg->reproduce_attempt_limit > 5)
      hmc_arg->reproduce_attempt_limit = 3;
    Ntests = 2;
  } else {
    hmc_arg->reproduce_attempt_limit = 1;
    Ntests = 1;
  }

  LatticeFactory::Destroy();

  return Ntests;

}

//!< Save final gauge field configuration
void AlgHmc::saveFinalState() {
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  lat.CopyGaugeField(gauge_field_final);
  LatticeFactory::Destroy();
}

//!< Restore the initial gauge field and rng state
void AlgHmc::restoreInitialState() {
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  lat.GaugeField(gauge_field_init);
//  LRG.SetStates(rng4d_init, FOUR_D);
//  LRG.SetStates(rng5d_init, FIVE_D);
  lrg_state.SetStates();
  LatticeFactory::Destroy();
}

void AlgHmc::shiftStates(int x, int y, int z, int t) {
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  GDS.Set(x,y,z,t);
  LRG.Shift();
  lat.Shift();
  LatticeFactory::Destroy();
}

//!< Restore the final gauge field
void AlgHmc::restoreFinalState() {
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  lat.GaugeField(gauge_field_final);
  LatticeFactory::Destroy();
}

//!< Check for reproducibility
int AlgHmc::reproduceTest(int attempt) {

  char *fname = "reproduceTest(int)" ;

  if (hmc_arg->reproduce == REPRODUCE_YES) {
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    //!< Compare the final gauge configs generated
    if (checksum[0] == checksum[1]) {
      VRB.Result(cname,fname,"Passed reproducibility test\n");
      LatticeFactory::Destroy();
      return 1;
    } else {
      VRB.Result(cname,fname, 
		 "Failed reproducibility first = %p, second = %p\n", 
		 checksum[0], checksum[1]);
      LatticeFactory::Destroy();
    }
    
    VRB.Result(cname,fname,"Failed reproducibility test %d\n",attempt);
    if (attempt == hmc_arg->reproduce_attempt_limit-1) 
      ERR.General(cname,fname,"Failed to reproduce\n");
    return 0;
 } else {
   return 1;
 }

}

void AlgHmc::reverseTest() {

  char *fname = "reverseTest()";

  //!< First calculate delta_delta_h
  Float delta_delta_h = delta_h - h_delta;
  delta_delta_h = sqrt(delta_delta_h*delta_delta_h);

  //!< Now calculate delta_delta_U
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

  ((Vector*)(lat.GaugeField()))->
    VecMinusEquVec((Vector*)gauge_field_init, g_size);
  Float delta_delta_U = ((Vector*)(lat.GaugeField()))->NormSqGlbSum4D(g_size);
  delta_delta_U = sqrt(delta_delta_U);

  LatticeFactory::Destroy();
  
  VRB.Result(cname, fname, 
	     "Reversibility check: delta dH = %e, delta dU = %e\n", 
	     delta_delta_h, delta_delta_U);

}


CPS_END_NAMESPACE
