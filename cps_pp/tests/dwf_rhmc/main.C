/* 
Quick DWF RHMC code, which measures the plaquette on each trajectory.
Performs a one flavour calculation, and includes the one flavour 
Pauli-Villars field to the action.
*/

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

USING_NAMESPACE_CPS

// some function prototypes
void setup_do_arg(DoArg& do_arg) ; 
void setup_hmd_arg(HmdArg& hmd_arg) ;
void setup_eig_arg(EigArg& eig_arg) ;

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

  char *cname = "dwf_rhmc";
  char *fname = "main";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  HmdArg hmd_arg;
  setup_hmd_arg(hmd_arg) ;
  EigArg eig_arg;
  setup_eig_arg(eig_arg);

  // parameters for the simulation
  const int no_measure_sweep = 1 ; 
  int sweep_counter = 0 ;
  const int total_measure = 10;
  
  // acceptance probability
  Float acceptance=0.0;

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------

  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  
  char *total_sites_st = "total sites = ";
  VRB.Flow(cname,fname,"%s%f\n",total_sites_st,IFloat(total_sites));

  const char *plaqfile = "plaquette.dat";
  const char *accfile = "acceptance.dat";

  FILE *plaq = Fopen(plaqfile,"w");
  FILE *acc = Fopen(accfile,"w");

  //----------------------------------------------------------------
  // Run Rational Hybrid Monte Carlo
  //----------------------------------------------------------------
  {
    // Wilson gauge action with dwf fermions
    GwilsonFdwf lat;

    {

      for (int i = 0; i < total_measure ; i += no_measure_sweep ) {

	if (i>100) hmd_arg.metropolis = METROPOLIS_YES;
	else hmd_arg.metropolis = METROPOLIS_NO;

        VRB.Flow(cname,fname,"iteration # = %d\n", i);

	//------------------------------------------------------------
	// Run asqtad RHMC
	//------------------------------------------------------------

	VRB.Flow(cname,fname,"AlgHmcRHMC starts....\n");
	{
 	  AlgHmcRHMC rhmc(lat,&common_arg,&hmd_arg,&eig_arg);
 	  //AlgHmcPhi rhmc(lat,&common_arg,&hmd_arg);
          for (int n = 0 ; n < no_measure_sweep; n++)
	    {
	      VRB.Flow(cname,fname,"HMD sweep n= %d/%d\n",n,no_measure_sweep) ;
	      Float accept = rhmc.run();
	      acceptance += accept;
	      sweep_counter++;
	      VRB.Flow(cname,fname,"%d Acceptance = %e\n",i, accept);
	      Fprintf(acc,"%d %e\n",i+n,accept);
	    }

	} // end of RHMC scope

	//----------------------------------------------------------------
	// calculate action and write it to file
	//----------------------------------------------------------------

	Float sum_plaq0 = lat.SumReTrPlaq();
	Float aver_plaq0 = sum_plaq0/(18.0*total_sites);
	VRB.Flow(cname,fname,"%d plaquette = %e\n",i,(Float)aver_plaq0); 
	Fprintf(plaq,"%d %0.16e\n",i,(Float)aver_plaq0);

      } // end of loop over i

    }
  }

  Fclose(plaq); 
  Fclose(acc);

  End();
  
  return(0);
}

  
void setup_do_arg(DoArg& do_arg)
{
  
  do_arg.x_node_sites = 4;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 4;
  do_arg.t_node_sites = 2;
  do_arg.s_node_sites = 8;
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
  do_arg.start_seed_kind = START_SEED_INPUT;
  do_arg.start_seed_value = 1213;
  do_arg.beta = 6.0;
  do_arg.dwf_height = 1.8;

  do_arg.xi_bare = 1;
  do_arg.xi_v = 1;
  do_arg.xi_v_xi = 1;

}

void setup_hmd_arg(HmdArg& hmd_arg)
{
  // Fermion fields first, then boson fields
  Float tau = 0.5;
  int frm = 1; // Number of fermion fields
  hmd_arg.n_frm_masses = 2;

  hmd_arg.frm_mass[0] = 0.04;
  hmd_arg.lambda_low[0] = 0.0004;
  hmd_arg.lambda_high[0] = 2.4;
  hmd_arg.frm_power_num[0] = 1;
  hmd_arg.frm_power_den[0] = 2;
  hmd_arg.FRatDeg[0] = 9;
  hmd_arg.SRatDeg[0] = 14;
  hmd_arg.field_type[0] = FERMION;

  hmd_arg.frm_mass[1] = 1.0;
  hmd_arg.lambda_low[1] = 0.02;
  hmd_arg.lambda_high[1] = 2.4;
  hmd_arg.frm_power_num[1] = 1;
  hmd_arg.frm_power_den[1] = 2;
  hmd_arg.FRatDeg[1] = 4;
  hmd_arg.SRatDeg[1] = 7;
  hmd_arg.field_type[1] = BOSON;

  hmd_arg.n_bsn_masses = 0;
  hmd_arg.steps_per_traj = 10;
  hmd_arg.step_size = tau/hmd_arg.steps_per_traj;
  hmd_arg.metropolis = METROPOLIS_YES;
  hmd_arg.reunitarize = REUNITARIZE_YES;
  hmd_arg.isz = 0; // Location of smallest polar shift
  hmd_arg.sw = 2; // Sexton-Weingarten term (gauge contribution per fermion)

  // Set the required degree of approximation

  hmd_arg.precision = 50;
  hmd_arg.approx_type = CONSTANT;
  hmd_arg.spread = 0.10;

  // set any other common variables
  for (int i=0; i<hmd_arg.n_frm_masses; i++) {
    hmd_arg.max_num_iter[i] = 5000;
    hmd_arg.stop_rsd[i] = 1.0e-6;
    hmd_arg.stop_rsd_md[i] = 1e-6;
    hmd_arg.stop_rsd_mc[i] = 1e-10;
    hmd_arg.valid_approx[i] = 0;
  }
  
}

void setup_eig_arg(EigArg& eig_arg)
{
  eig_arg.N_eig = 1;
  eig_arg.RsdR_a = 1e-3;
  eig_arg.RsdR_r = 1e-3;
  eig_arg.Rsdlam = 1e-3;
  eig_arg.Kalk_Sim = 0;
  eig_arg.N_min = 0;
  eig_arg.N_max = 0;
  eig_arg.N_KS_max = 0;
  eig_arg.n_renorm = 15;
  eig_arg.Cv_fact = 0;
  eig_arg.MaxCG = 5000;
  eig_arg.ProjApsiP = 0;
  eig_arg.hsum_dir = 0;
  eig_arg.print_hsum = 0;
}

