#include<config.h>

//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-06-02 09:36:41 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_dwf_pbp/main.C,v 1.3 2004-06-02 09:36:41 zs Exp $
//  $Id: main.C,v 1.3 2004-06-02 09:36:41 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_dwf_pbp/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/random.h>
#include<alg/alg_hmd.h>
#include<alg/alg_pbp.h>
#include<alg/alg_plaq.h>
#include<alg/do_arg.h>
#include<alg/no_arg.h>

CPS_START_NAMESPACE
GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;
CPS_END_NAMESPACE

USING_NAMESPACE_CPS

int main(int argc,char *argv[])
{

  FILE *fp;
  
  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

#ifdef PARALLEL
  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 4;
  do_arg.t_node_sites = 4;
  do_arg.s_node_sites = 4;
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
  do_arg.s_nodes = 1;
#else
  do_arg.x_node_sites = 4;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 4;
  do_arg.t_node_sites = 4;
  do_arg.s_node_sites = 4;
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
  do_arg.s_nodes = 1;
#endif

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_APRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_conf_load_addr = (Matrix *) 0x12f13;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.start_seed_value = 1111;
  do_arg.beta = 5.8;
  do_arg.dwf_height = 0.9;

#if TARGET==cpsMPI
    MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes);    
    using MPISCU::fprintf;
    using MPISCU::printf;
#endif

  GJP.Initialize(do_arg);




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
  pbp_arg.src_l_s = GJP.SnodeSites() * GJP.Snodes() - 1;
  pbp_arg.snk_u_s = GJP.SnodeSites() * GJP.Snodes() - 1;
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
  int num_hits = 1;
  int num_meas = 1;
  int num_sweeps = 1;
  int min_ls = 3;
  int num_ls = 2;
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

  if( (fp = fopen("info.dat", "a")) == NULL ) {
    ERR.FileA(" ","main", "info.dat");
  }
  fprintf(fp, "Total sites :    %d\n", total_sites); 
  fclose(fp);


  //----------------------------------------------------------------
  // Run AlgPlaq 
  //----------------------------------------------------------------
  {
    GwilsonFnone lat;
    common_arg.results = CAST_AWAY_CONST("plaq.dat");
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

      if( (fp = fopen("info.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "info.dat");
      }
      fprintf(fp,"%d %d %f\n",
	      therm, cur_g_upd_cnt, 
	      float(cur_g_upd_cnt) / float(therm) );
      fclose(fp);

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
      
      if( (fp = fopen("info.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "info.dat");
      }
      fprintf(fp,"%d %d %f\n",
	      therm, cur_g_upd_cnt, 
	      float(cur_g_upd_cnt) / float(therm) );
      fclose(fp);

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

	if( (fp = fopen("info.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "info.dat");
	}
	fprintf(fp,"%d %d %f\n",
		iter, cur_g_upd_cnt, 
		float(cur_g_upd_cnt) / float(iter) );	
	fclose(fp);
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
	  if(num_hits != 0)
	    {
	      common_arg.results = CAST_AWAY_CONST("pbp.dat");
	      pbp_arg.mass[0] = mass;
	      pbp_arg.src_u_s = 0;
	      pbp_arg.src_l_s = ls-1;
	      pbp_arg.snk_u_s = ls-1;
	      pbp_arg.snk_l_s = 0;
	      AlgPbp pbp(lat,&common_arg,&pbp_arg);
	      for(int i=0; i<num_hits; i++){
		pbp.run();
	      }
	    }
	}

	ls = ls + ls_step;
      }
      mass = mass + mass_step;
    }



  }

  return(0);
}






