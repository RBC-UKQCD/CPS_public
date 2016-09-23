/*
  $Id: main.C,v 1.20 2008/02/08 18:35:08 chulwoo Exp $
*/

/* Quick Asqtad Monte Carlo code, which measures the plaquette on each trajectory. */

#include <util/qcdio.h>
#include <config.h>

#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <alg/alg_hmd.h>
#include <alg/do_arg.h>
#include <util/random.h>
#include <comms/sysfunc_cps.h>

const int nx = 4;
const int ny = 4;
const int nz = 4;
const int nt = 4;

const int SAVE_DOARG = 1;
const int LOAD_DOARG = 0;

USING_NAMESPACE_CPS


// some function prototypes
void setup_do_arg(DoArg& do_arg) ; 
void setup_hmd_arg(HmdArg& hmd_arg) ;


int main(int argc,char *argv[])
{
  Start();

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;
  if (LOAD_DOARG) do_arg.Decode("asqtad_hmc.in","do_arg");
  else setup_do_arg(do_arg); 
  if (SAVE_DOARG) do_arg.Encode("asqtad_hmc.out","do_arg");
  GJP.Initialize(do_arg);

  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------

  VRB.Level(0);
  VRB.ActivateLevel(VERBOSE_RESULT_LEVEL);
  char *cname = "asqtad_hmd_r";
  char *fname = "main";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  HmdArg hmd_arg;
  setup_hmd_arg(hmd_arg) ;


  // parameters for the simulation
  const int no_warmup_sweep = 0 ; 
  const int no_measure_sweep = 1 ; 
  int sweep_counter = 0 ;
  const int total_measure = 10;
  
  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------


  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  
  char *total_sites_st = "total sites = ";
  VRB.Flow(cname,fname,"%s%f\n",total_sites_st,IFloat(total_sites));
  FILE *plaq;
  const char *filename = "plaquette.dat";
  plaq = Fopen(filename,"w");
  
  //----------------------------------------------------------------
  // Run Wilson Monte Carlo
  //----------------------------------------------------------------
  {
    GwilsonFasqtad lat;

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
	  VRB.Flow(cname,fname,"%d plaquette = %0.16e\n",n,(Float)aver_plaq0); 
	  Fprintf(plaq,"%d %0.16e\n",n,(Float)aver_plaq0);

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
	Fprintf(plaq,"%d %0.16e\n",i,(Float)aver_plaq0);


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

  Fclose(plaq); 
  }
  End();
  return(0);
}


void setup_do_arg(DoArg& do_arg)
{

  int size[4]; 
  int sites[4];
  size[0] = SizeX(); sites[0] =nx;
  size[1] = SizeY(); sites[1] =ny;
  size[2] = SizeZ(); sites[2] =nz;
  size[3] = SizeT(); sites[3] =nt;
  for (int i = 0;i<4;i++) 
    if( (size[i]*2) > sites[i] ) size[i] = sites[i]/2;
  for (int i = 0;i<4;i++) 
    if( sites[i] % size[i] !=0) size[i] = 1;
  do_arg.x_node_sites = nx/size[0];
  do_arg.y_node_sites = ny/size[1];
  do_arg.z_node_sites = nz/size[2];
  do_arg.t_node_sites = nt/size[3];
  do_arg.x_nodes = size[0];
  do_arg.y_nodes = size[1];
  do_arg.z_nodes = size[2];
  do_arg.t_nodes = size[3];
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 5.6;
  do_arg.xi_bare = 1;
  do_arg.xi_v = 1;
  do_arg.xi_v_xi = 1;
  do_arg.asqtad_KS = (1.0/8.0)+(6.0/16.0)+(1.0/8.0);
  do_arg.asqtad_naik = -1.0/24.0;
  do_arg.asqtad_lepage = -1.0/16;
  do_arg.asqtad_3staple = (1.0/8.0)*0.5;
  do_arg.asqtad_5staple = ( 1.0/8.0)*0.25*0.5;
  do_arg.asqtad_7staple = (1.0/8.0)*0.125*(1.0/6.0);
}



void setup_hmd_arg(HmdArg& hmd_arg)
{

  hmd_arg.n_frm_masses = 1;
  hmd_arg.frm_mass[0] = 0.1;
  hmd_arg.frm_flavors[0] = 2;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.max_num_iter[0] = 5000;
  hmd_arg.stop_rsd[0] = 1.0E-06;
  hmd_arg.step_size = 0.002;
  hmd_arg.steps_per_traj = 25;
//  hmd_arg.step_size = 0.005;
//  hmd_arg.steps_per_traj = 100;
  hmd_arg.metropolis = METROPOLIS_YES;
  hmd_arg.reunitarize = REUNITARIZE_YES;

}
