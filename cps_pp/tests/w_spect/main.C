#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2007-06-25 21:39:14 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/w_spect/main.C,v 1.11 2007-06-25 21:39:14 chulwoo Exp $
//  $Id: main.C,v 1.11 2007-06-25 21:39:14 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.11 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/w_spect/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <util/qcdio.h>
#include<util/lattice.h>
#include<alg/alg_hmd.h>
#include<alg/alg_pbp.h>
#include<alg/alg_w_spect.h>
#include<alg/alg_plaq.h>
#include<alg/do_arg.h>
#include<alg/no_arg.h>


USING_NAMESPACE_CPS


int main(int argc,char *argv[])
{

  Start();
  FILE *fp;
  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 4;
  do_arg.s_node_sites = 4;
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
  do_arg.s_nodes = 1;

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_conf_load_addr = 0x12f13;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.start_seed_value = 1111;
  do_arg.beta = 5.8;
  do_arg.dwf_height = 0.9;

#if TARGET==cpsMPI
    MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes);
    using MPISCU::printf;
    using MPISCU::fprintf;
#endif

  GJP.Initialize(do_arg);

  //----------------------------------------------------------------
  // Initialize arguments for the new Wilson spectrum code. -- Ping
  //----------------------------------------------------------------
  AlgWspect::SetCounter(1 /* 1st measurement */, 1 /* stride */);
  // The above is also the default behavior.

  WspectArg w_spect_arg;
  CgArg w_spect_cg_arg;
  w_spect_cg_arg.mass = 0.1;
  w_spect_cg_arg.stop_rsd = 1.0E-12;
  w_spect_cg_arg.max_num_iter = 5000;
  w_spect_arg.prop_dir = 3;
  w_spect_arg.num_mom = 1;              // zero spatial momentum only
  w_spect_arg.source_kind = POINT_W; 
  w_spect_arg.aots_start = 0;
  w_spect_arg.aots_step = 1;
  w_spect_arg.aots_num = 1;
  w_spect_arg.baryons_on = 1;
  w_spect_arg.normal_mesons_on = 1;

  // this scales the quark source by 1, the default
  // mechanism seems to be broken in w_quark
  w_spect_arg.rescale_factor = 1.0 ;
  

  WspectOutput wout;
  wout.fold = BARYON_PAST;      // BARYON_FOLD, BARYON_RAW, BARYON_PAST
  wout.cg = "w_spect.dat";
  wout.pbp = "w_pbp.dat";

  wout.a0 = "a0.dat";
  wout.a0_prime = "a0_prime.dat";
  wout.a1_x = "a1_x.dat";
  wout.a1_y = "a1_y.dat";
  wout.a1_z = "a1_z.dat";
  wout.b1_x = "b1_x.dat"; 
  wout.b1_y = "b1_y.dat"; 
  wout.b1_z = "b1_z.dat";  
  wout.pion = "pion.dat";  
  wout.pion_prime = "pion_prime.dat"; 
  wout.rho_x = "rho_x.dat";  	
  wout.rho_y = "rho_y.dat";  	
  wout.rho_z = "rho_z.dat";  	 	    
  wout.rho_x_prime = "rho_x_prime.dat";  	 	      	    
  wout.rho_y_prime = "rho_y_prime.dat";  	 	     
  wout.rho_z_prime = "rho_z_prime.dat";  	 	      	     
  wout.nucleon = "nucleon.dat";  	 	      	     
  wout.nucleon_prime = "nucleon_prime.dat";  	 	      	           
  wout.delta_x = "delta_x.dat";  	 	      	      	     
  wout.delta_y = "delta_y.dat";  	 	      	      	      	     
  wout.delta_z = "delta_z.dat";  	 	      	      	      	     
  wout.delta_t = "delta_t.dat";  	 	      	      	      	     


  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;

  PbpArg pbp_arg;
  pbp_arg.pattern_kind = ARRAY;
  pbp_arg.n_masses = 1;
  pbp_arg.mass[0] = 0.1;
  pbp_arg.stop_rsd = 1.0E-12;
  pbp_arg.max_num_iter = 5000;
  pbp_arg.src_u_s = 0;
  pbp_arg.src_l_s = GJP.SnodeSites()-1;
  pbp_arg.snk_u_s = GJP.SnodeSites()-1;
  pbp_arg.snk_l_s = 0;
  pbp_arg.snk_loop = 0;


  HmdArg hmd_arg;
  hmd_arg.n_frm_masses = 1;
  hmd_arg.frm_mass[0] = 0.1;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.max_num_iter[0] = 5;
  hmd_arg.stop_rsd[0] = 1.0E-12;
  hmd_arg.step_size = 0.01;
  hmd_arg.steps_per_traj = 50;
  hmd_arg.metropolis = METROPOLIS_YES;
  hmd_arg.reunitarize = REUNITARIZE_YES;

  NoArg plaq_arg;


  //----------------------------------------------------------------
  // Initialize run parameters
  //----------------------------------------------------------------
  int num_therm_no_met = 1;
  int num_therm_met = 1;
  int pbp_num_hits = 1;
  int spect_yes = 1;
  int num_meas = 1;
  int num_sweeps = 1;
  int min_ls = 3;
  int num_ls = 1;
  int ls_step = 1;
  Float mass_start = 0.0;
  Float mass_step = -0.0;
  int num_mass = 1;


  //----------------------------------------------------------------
  // Initialize counters
  //----------------------------------------------------------------
  int iter = 0;
  int init_g_upd_cnt = 0;
  int cur_g_upd_cnt = 0;

  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();

  if( (fp = Fopen("info.dat", "a")) == NULL ) {
    ERR.FileA(" ","main", "info.dat");
  }
  Fprintf(fp, "Total sites :    %d\n", total_sites); 
  Fclose(fp);

  //----------------------------------------------------------------
  // Run AlgPlaq 
  //----------------------------------------------------------------
  {
    GwilsonFnone lat;
    common_arg.set_filename("plaq.dat");
    AlgPlaq plaq(lat,&common_arg,&plaq_arg);
    
    plaq.run();
  }


  //----------------------------------------------------------------
  // Thermalize gauge field using HMC Phi with no metropolis
  //----------------------------------------------------------------
  hmd_arg.metropolis = METROPOLIS_NO;
  if(num_therm_no_met != 0)
    {
      GwilsonFnone lat;
      common_arg.results = 0;
      AlgHmcPhi hmc(lat,&common_arg,&hmd_arg);
      
      init_g_upd_cnt = lat.GupdCnt();
      int therm;
      for(therm=0; therm<num_therm_no_met; therm++){
	hmc.run();
	cur_g_upd_cnt = lat.GupdCnt() - init_g_upd_cnt;
      }
      if( (fp = Fopen("info.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "info.dat");
      }
      Fprintf(fp,"%d %d %f\n",
	      therm, cur_g_upd_cnt, 
	      Float(cur_g_upd_cnt) / Float(therm) );
      Fclose(fp);

      init_g_upd_cnt = lat.GupdCnt();
   }


  //----------------------------------------------------------------
  // Thermalize gauge field using HMC Phi with metropolis
  //----------------------------------------------------------------
  hmd_arg.metropolis = METROPOLIS_YES;
  if(num_therm_met != 0)
    {
      GwilsonFnone lat;
      common_arg.results = 0;
      AlgHmcPhi hmc(lat,&common_arg,&hmd_arg);
      
      init_g_upd_cnt = lat.GupdCnt();
      int therm;
      for(therm=0; therm<num_therm_met; therm++){
	hmc.run();
	cur_g_upd_cnt = lat.GupdCnt() - init_g_upd_cnt;
      }
      
      if( (fp = Fopen("info.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "info.dat");
      }
      Fprintf(fp,"%d %d %f\n",
	      therm, cur_g_upd_cnt, 
	      Float(cur_g_upd_cnt) / Float(therm) );
      Fclose(fp);

      init_g_upd_cnt = lat.GupdCnt();
    }


  //----------------------------------------------------------------
  // Loop over the number of measurements
  //----------------------------------------------------------------
  for(int meas = 0; meas< num_meas; meas++){

    //--------------------------------------------------------------
    // Evolve gauge field using HMC Phi
    //--------------------------------------------------------------
    hmd_arg.metropolis = METROPOLIS_YES;
    if(num_sweeps != 0)
      {
	GwilsonFnone lat;
	common_arg.results = 0;
	AlgHmcPhi hmc(lat,&common_arg,&hmd_arg);
	
	for(int sweeps=0; sweeps<num_sweeps; sweeps++){
	  hmc.run();
	  iter = iter + 1;
	}
	
	cur_g_upd_cnt = lat.GupdCnt() - init_g_upd_cnt;

	if( (fp = Fopen("info.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "info.dat");
	}
	Fprintf(fp,"%d %d %f\n",
		iter, cur_g_upd_cnt, 
		Float(cur_g_upd_cnt) / Float(iter) );	
	Fclose(fp);
      }


    //----------------------------------------------------------------
    // Loop over mass and ls 
    //----------------------------------------------------------------
    Float mass = mass_start;
    for(int m=0; m<num_mass; m++){
      int ls = min_ls;
      for(int ils=0; ils<num_ls; ils++){
	GJP.SnodeSites(ls);


	GwilsonFdwf lat;
	{
	  //------------------------------------------------------------
	  // Run AlgPbp 
	  //------------------------------------------------------------
	  if(pbp_num_hits != 0)
	    {
	      common_arg.set_filename("pbp.dat");
	      pbp_arg.mass[0] = mass;
	      pbp_arg.src_u_s = 0;
	      pbp_arg.src_l_s = ls-1;
	      pbp_arg.snk_u_s = ls-1;
	      pbp_arg.snk_l_s = 0;
	      AlgPbp pbp(lat,&common_arg,&pbp_arg);
	      for(int i=0; i<pbp_num_hits; i++){
		pbp.run();
	      }
	    }

	  //------------------------------------------------------------
	  // Run spectroscopy
	  //------------------------------------------------------------
	  if(spect_yes == 1){
	    {
	      common_arg.results = &wout;
	      w_spect_cg_arg.mass = mass;
	      AlgWspect ws(lat,&common_arg,&w_spect_arg,&w_spect_cg_arg);
	      
	      ws.run();
	    }
	  }


	}

	ls = ls + ls_step;
      }
      mass = mass + mass_step;
    }



  }
  End();
  return(0);
}





