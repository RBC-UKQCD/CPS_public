/* Quick staggered R2 code, which measures the plaquette on each trajectory. */

#include <util/qcdio.h>
#include <stdlib.h>	// exit()
#include <config.h>
#include <math.h>

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
#include <util/vector.h>
#include <alg/alg_remez.h>
#include <util/random.h>

const int nx = 4;
const int ny = 4;
const int nz = 4;
const int nt = 4;

USING_NAMESPACE_CPS

// some function prototypes
void setup_do_arg(DoArg& do_arg) ; 
void setup_asqtad_arg(DoArg& do_arg, Float plaq) ; 
void setup_hmd_arg(HmdArg& hmd_arg) ;

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
  VRB.Level(0);
  VRB.ActivateLevel(VERBOSE_RESULT_LEVEL);

  char *cname = "stag_r2";
  char *fname = "main";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  HmdArg hmd_arg;
  setup_hmd_arg(hmd_arg) ;

  // parameters for the simulation
  const int no_measure_sweep = 1 ; 
  int sweep_counter = 0 ;
  const int total_measure = 10;
  
  // average plaquette
  Float plaquette=0.0;

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------

  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  
  char *total_sites_st = "total sites = ";
  VRB.Flow(cname,fname,"%s%f\n",total_sites_st,IFloat(total_sites));

  const char *plaqfile = "plaquette_r2.dat";

  FILE *plaq = Fopen(plaqfile,"w");

  //----------------------------------------------------------------
  // Run Rational Hybrid Monte Carlo
  //----------------------------------------------------------------
  {
    // Wilson gauge action with staggered fermions
    GwilsonFstag lat;

    {

      for (int i = 0; i < total_measure ; i += no_measure_sweep ) {

	if (i>100) hmd_arg.metropolis = METROPOLIS_YES;
	else hmd_arg.metropolis = METROPOLIS_NO;

        VRB.Flow(cname,fname,"iteration # = %d\n", i);

	//------------------------------------------------------------
	// Run staggered R2
	//------------------------------------------------------------

	VRB.Flow(cname,fname,"AlgHmdR2 starts....\n");
	{
 	  AlgHmdR2 r2(lat,&common_arg,&hmd_arg);
          for (int n = 0 ; n < no_measure_sweep; n++)
	    {
	      VRB.Flow(cname,fname,"HMD sweep n= %d/%d\n",n,no_measure_sweep) ;
	      r2.run();
	      sweep_counter++;
	    }

	} // end of R2 scope

	//----------------------------------------------------------------
	// calculate action and write it to file
	//----------------------------------------------------------------

	Float sum_plaq0 = lat.SumReTrPlaq();
	Float aver_plaq0 = sum_plaq0/(18.0*total_sites);
	VRB.Flow(cname,fname,"%d plaquette = %e\n",i,(Float)aver_plaq0); 
	Fprintf(plaq,"%d %0.16e\n",i,(Float)aver_plaq0);
	plaquette += aver_plaq0;

	// Dynamic setting of u0
	//do_arg.u0 = pow(3.0*aver_plaq0/(Float)(i+1),-0.25);
	VRB.Flow(cname,fname,"u0 = %e\n", do_arg.u0);

      } // end of loop over i

    }
  }

  Fclose(plaq); 

  End();
  
  return(0);
}

  
void setup_do_arg(DoArg& do_arg)
{
  
  do_arg.x_node_sites = 4;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 4;
  do_arg.t_node_sites = 4;
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = 1;

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_INPUT;
  do_arg.start_seed_value = 1213;
  do_arg.beta = 4.8;
  do_arg.dwf_height = 1.8;
  do_arg.clover_coeff = 2.3;

  do_arg.xi_bare = 1;
  do_arg.xi_v = 1;
  do_arg.xi_v_xi = 1;

  setup_asqtad_arg(do_arg, 1.0);

}

void setup_asqtad_arg(DoArg& do_arg, Float plaq) {
  do_arg.u0 = pow(3.0*plaq,-0.25);
  do_arg.asqtad_KS      = (1.0/8.0)+(3.0/8.0)+(1.0/8.0);
  do_arg.asqtad_naik    = (-1.0/24.0)*pow(do_arg.u0,-2);
  do_arg.asqtad_lepage  = (-1.0/16.0)*pow(do_arg.u0,-4);
  do_arg.asqtad_3staple = (-1.0/8.0)*(1.0/2.0)*pow(do_arg.u0,-2);
  do_arg.asqtad_5staple = ( 1.0/8.0)*(1.0/8.0)*pow(do_arg.u0,-4);
  do_arg.asqtad_7staple = (-1.0/8.0)*(1.0/48.0)*pow(do_arg.u0,-6);
}

void setup_hmd_arg(HmdArg& hmd_arg)
{
  Float tau = 1.0;
  hmd_arg.n_frm_masses = 2;
  hmd_arg.frm_mass[0] = 0.25;
  hmd_arg.frm_mass[1] = 0.05;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.steps_per_traj = 10;
  hmd_arg.step_size = tau/hmd_arg.steps_per_traj;
  hmd_arg.metropolis = METROPOLIS_YES;
  hmd_arg.reunitarize = REUNITARIZE_YES;
    
  // set any other common variables
  hmd_arg.max_num_iter[0] = 5000;
  hmd_arg.max_num_iter[1] = 5000;
  hmd_arg.stop_rsd[0] = 1.0E-8;
  hmd_arg.stop_rsd[1] = 1.0E-8;
  
}
