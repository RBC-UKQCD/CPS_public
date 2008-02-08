/* Quick Wilson Monte Carlo code, which measures the plaquette on each trajectory. */

#include <util/random.h>

#include <stdio.h>
#include <stdlib.h>	// exit()
#include <config.h>

#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qcdio.h>
#include <alg/alg_hmd.h>
#include <alg/alg_s_spect.h>
#include <alg/do_arg.h>
#include <alg/common_arg.h>
#include <alg/cg_arg.h>
#include <alg/hmd_arg.h>
#include <alg/s_spect_arg.h>
#include <alg/fix_gauge_arg.h>
#include <alg/alg_fix_gauge.h>
#include <alg/aots_s.h>
#include <comms/sysfunc_cps.h>
#include <util/vector.h>
#include <util/random.h>
#include <time.h>
#include <sys/time.h>

USING_NAMESPACE_CPS;

// some function prototypes
void setup_do_arg(DoArg& do_arg) ; 
void setup_hmd_arg(HmdArg& hmd_arg) ;

void RerouteStdio(char *directory) ;
int main(int argc,char *argv[])
{
  Start();
  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;
  setup_do_arg(do_arg); 
  GJP.Initialize(do_arg);

  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------
//  VRB.Level(GJP.VerboseLevel());
  VRB.Level(0);
  VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);
  VRB.ActivateLevel(VERBOSE_RESULT_LEVEL);
  char *cname = "wilson_hmc";
  char *fname = "main";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  common_arg.results = (void *)"HMC_clover_output_512.dat";
  HmdArg hmd_arg;
  setup_hmd_arg(hmd_arg) ;


  // parameters for the simulation
  const int no_warmup_sweep = 0 ; 
  const int no_measure_sweep = 1 ; 
  int sweep_counter = 0 ;
  const int total_measure = 6000;
  
  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------


  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  
  char *total_sites_st = "total sites = ";
  VRB.Flow(cname,fname,"%s%f\n",total_sites_st,IFloat(total_sites));
#ifdef PARALLEL
  VRB.Flow(cname,fname,"Nodes = %d %d %d %d\n",SizeX(), SizeY(), SizeZ(), SizeT());
#endif

  FILE *plaquette;
    plaquette = Fopen("plaquette.dat", "w");
    printf("Plaquette fopen returned %x",plaquette);
    if ( plaquette == NULL ) { perror("oops"); exit(-1); }
  printf("Running monte carlo\n");
  //----------------------------------------------------------------
  // Run Wilson Monte Carlo
  //----------------------------------------------------------------
  {
    GwilsonFclover lat;

    {

      //-----------------------------------------------------------------
      // warming up 
      //-----------------------------------------------------------------
      {
	hmd_arg.metropolis = METROPOLIS_NO;
	AlgHmcPhi hmc_r(lat,&common_arg,&hmd_arg);
        for (int n = 0 ; n < no_warmup_sweep ; n++) {
          VRB.Flow(cname,fname,"warmup iteration # = %d/%d\n",n,no_warmup_sweep);
	  Float sum_plaq0 = lat.SumReTrPlaq();
	  Float aver_plaq0 = sum_plaq0/(18.0*total_sites);
	  VRB.Flow(cname,fname,"%d plaquette = %e\n",n,(Float)aver_plaq0); 
	    Fprintf(plaquette,"%d %e\n",n,(Float)aver_plaq0);
	  //  fflush(plaquette);
          unsigned long *binary = (unsigned long *)&aver_plaq0;
          printf("%d %e %x:%x\n",n,(Float)aver_plaq0,binary[0],binary[1]);
	  fflush(stdout);
       	  hmc_r.run();
	  sweep_counter++; 
        }

      }
      for (int i = 0; i < total_measure ; i ++ ) {

        VRB.Flow(cname,fname,"iteration # = %d\n", i);

	//----------------------------------------------------------------
	// calculate action and write it to file
	//----------------------------------------------------------------

	Float sum_plaq0 = lat.SumReTrPlaq();
	Float aver_plaq0 = sum_plaq0/(18.0*total_sites);
	VRB.Flow(cname,fname,"%d plaquette = %e\n",i,(Float)aver_plaq0); 
	printf("%d %e\n",i,(Float)aver_plaq0);
	  Fprintf(plaquette,"%d %e\n",i,(Float)aver_plaq0);
	printf("Written\n");

	//------------------------------------------------------------
	// Run wilson HMC
	//------------------------------------------------------------

	VRB.Flow(cname,fname,"AlgHmcPhi starts....\n");
	{
	  hmd_arg.metropolis = METROPOLIS_NO;
	  AlgHmcPhi hmc_r(lat,&common_arg,&hmd_arg);
          for (int n = 0 ; n < no_measure_sweep; n++)
	    {
	      VRB.Flow(cname,fname,"HMD sweep n= %d/%d\n",n,no_measure_sweep) ;
	      hmc_r.run();
	      sweep_counter++;
	    }

	} // end of HMC scope


      } // end of loop over i


    } 


  }
    Fclose(plaquette);

  return(0);
}

const int lx = 16;
const int ly = 16;
const int lz = 16;
const int lt = 16;

void setup_do_arg(DoArg& do_arg)
{


  do_arg.x_node_sites = lx/SizeX();
  do_arg.y_node_sites = ly/SizeY();
  do_arg.z_node_sites = lz/SizeZ();
  do_arg.t_node_sites = lt/SizeT();
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = SizeS();
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
//  do_arg.colors = 3;
  do_arg.beta = 5.2;
  do_arg.clover_coeff = 2.0171;
//  do_arg.verbose_level = 10;
  do_arg.xi_bare = 1;
  do_arg.xi_v = 1;
  do_arg.xi_v_xi = 1;

}



void setup_hmd_arg(HmdArg& hmd_arg)
{

  hmd_arg.n_frm_masses = 1;
  Float kappa = 0.1350;
  hmd_arg.frm_mass[0] = 1/(2*kappa) - 4.0;
  hmd_arg.frm_flavors[0] = 2;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.max_num_iter[0] = 1000;
  hmd_arg.stop_rsd[0] = 1.0E-8;
  hmd_arg.step_size = 0.01;
  hmd_arg.steps_per_traj = 50;
  hmd_arg.metropolis = METROPOLIS_NO;
  hmd_arg.reunitarize = REUNITARIZE_YES;

}
