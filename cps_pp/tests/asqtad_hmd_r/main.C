#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-09-21 20:16:50 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/asqtad_hmd_r/main.C,v 1.10 2004-09-21 20:16:50 chulwoo Exp $
//  $Id: main.C,v 1.10 2004-09-21 20:16:50 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/asqtad_hmd_r/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/* R algorithm with Asqtad fermions and Wilson gauge action
   which measures the plaquette on each trajectory. */




#include <util/lattice.h>
#include <util/gjp.h>

#include <util/qcdio.h>
#include <alg/alg_hmd.h>
#include <alg/do_arg.h>


const int nx = 4;
const int ny = 4;
const int nz = 4;
const int nt = 4;


USING_NAMESPACE_CPS



// some function prototypes
void setup_do_arg(DoArg& do_arg) ; 
void setup_hmd_arg(HmdArg& hmd_arg) ;

int main(int argc,char *argv[])
{
#if TARGET==QCDOC
  DefaultSetup(); 
#endif


  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;
  setup_do_arg(do_arg); 
#if TARGET==cpsMPI
  MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes);
  using MPISCU::fprintf;
  using MPISCU::printf;
#endif
  GJP.Initialize(do_arg);

  VRB.DeactivateAll();
  VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);
  VRB.ActivateLevel(VERBOSE_CLOCK_LEVEL);
  

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
    const int total_measure = 2 ;
  
  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------


  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  
  char *total_sites_st = "total sites = ";
  VRB.Flow(cname,fname,"%s%f\n",total_sites_st,IFloat(total_sites));
  FILE *plaq;
#if TARGET == QCDOC
  char filename[200];
  const char *plaq_filename = CWDPREFIX("plaquette");
  sprintf(filename, "%s%d%d%d%d%d%d.test2",
	plaq_filename,CoorX(), CoorY(), CoorZ(), CoorT(), CoorS(), CoorW());
#else
  const char *filename = "plaquette.dat";
#endif
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
	AlgHmdR hmc_r(lat,&common_arg,&hmd_arg);
        for (int n = 0 ; n < no_warmup_sweep ; n++) {
          VRB.Flow(cname,fname,"warmup iteration # = %d/%d\n",n,no_warmup_sweep);
	  Float sum_plaq0 = lat.SumReTrPlaq();
	  Float aver_plaq0 = sum_plaq0/(18.0*total_sites);
	  VRB.Flow(cname,fname,"%d plaquette = %0.16e\n",n,(Float)aver_plaq0); 
	  Fprintf(plaq,"%d %0.16e\n",n,(Float)aver_plaq0);fflush(plaq);

       	  hmc_r.run();
	  sweep_counter++; 
        }

      }
      for (int i = 0; i < total_measure ; i ++ ) {

        VRB.Flow(cname,fname,"iteration # = %d\n", i);

	//----------------------------------------------------------------
	// calculate plaquette and write it to file
	//----------------------------------------------------------------

	Float sum_plaq0 = lat.SumReTrPlaq();
	Float aver_plaq0 = sum_plaq0/(18.0*total_sites);
	VRB.Flow(cname,fname,"%d plaquette = %e\n",i,(Float)aver_plaq0); 
	Fprintf(plaq,"%d %0.16e\n",i,(Float)aver_plaq0);
	fflush(plaq);


	//------------------------------------------------------------
	// Run HMC
	//------------------------------------------------------------

	VRB.Flow(cname,fname,"AlgHmcR starts....\n");
	{
	  hmd_arg.metropolis = METROPOLIS_NO;
	  AlgHmdR hmc_r(lat,&common_arg,&hmd_arg);
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

  return(0);
}


void setup_do_arg(DoArg& do_arg)
{


#if TARGET==QCDOC
  do_arg.x_node_sites = nx/SizeX();
  do_arg.y_node_sites = ny/SizeY();
  do_arg.z_node_sites = nz/SizeZ();
  do_arg.t_node_sites = nt/SizeT();
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
#else
  do_arg.x_node_sites = nx ;
  do_arg.y_node_sites = ny ;
  do_arg.z_node_sites = nz ;
  do_arg.t_node_sites = nt ;
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
#endif
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 5.6;

  do_arg.asqtad_KS = (1.0/8.0)+(6.0/16.0)+(1.0/8.0);
  do_arg.asqtad_naik = -1.0/24.0;
  do_arg.asqtad_lepage = -1.0/16;
  do_arg.asqtad_3staple = (-1.0/8.0)*0.5;
  do_arg.asqtad_5staple = ( 1.0/8.0)*0.25*0.5;
  do_arg.asqtad_7staple = (-1.0/8.0)*0.125*(1.0/6.0);
#if 0
  do_arg.asqtad_KS = 1e-20;
  do_arg.asqtad_naik = 1.0;
  do_arg.asqtad_lepage = 1e-20;
  do_arg.asqtad_3staple = 1e-20;
  do_arg.asqtad_5staple = 1e-20;
  do_arg.asqtad_7staple = 1e-20;
#endif
}



void setup_hmd_arg(HmdArg& hmd_arg)
{

  hmd_arg.n_frm_masses = 1;
  hmd_arg.frm_mass[0] = 0.1;
  hmd_arg.frm_flavors[0] = 2;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.max_num_iter[0] = 500;
  hmd_arg.stop_rsd[0] = 1.0E-7;
  hmd_arg.step_size = 0.02;
  hmd_arg.steps_per_traj = 50;
  hmd_arg.metropolis = METROPOLIS_NO;
  hmd_arg.reunitarize = REUNITARIZE_YES;

}
