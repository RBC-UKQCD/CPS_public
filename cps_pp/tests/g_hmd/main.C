#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:10 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/g_hmd/main.C,v 1.7 2004-08-18 11:58:10 zs Exp $
//  $Id: main.C,v 1.7 2004-08-18 11:58:10 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.7 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/g_hmd/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#include<util/lattice.h>
#include<alg/alg_hmd.h>
#include<alg/alg_plaq.h>
#include<alg/alg_rect.h>
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
  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  do_arg.x_node_sites = 4;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 4;
  do_arg.t_node_sites = 4;
  do_arg.s_node_sites = 2;

#ifdef PARALLEL
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 2;
  do_arg.z_nodes = 2;
  do_arg.t_nodes = 2;
  do_arg.s_nodes = 1;
#else
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
  do_arg.s_nodes = 1;
#endif 

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 4.8;
  do_arg.c_1 = -0.05;
  do_arg.dwf_height = 0.9;
  do_arg.power_plaq_cutoff = 1.0;
  do_arg.power_plaq_exponent = 2;
  do_arg.power_rect_cutoff = 0.9;
  do_arg.power_rect_exponent = 4;

  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg_hmc;
  common_arg_hmc.set_filename("hmc.dat");

  CommonArg common_arg_plaq;
  common_arg_plaq.set_filename("plaq.dat");

  CommonArg common_arg_p_plaq;
  common_arg_p_plaq.set_filename("p_plaq.dat");

  CommonArg common_arg_rn_plaq;
  common_arg_rn_plaq.set_filename("rn_plaq.dat");

  CommonArg common_arg_rn_rect;
  common_arg_rn_rect.set_filename("rn_rect.dat");

  CommonArg common_arg_pr_plaq;
  common_arg_pr_plaq.set_filename("pr_plaq.dat");

  CommonArg common_arg_pr_rect;
  common_arg_pr_rect.set_filename("pr_rect.dat");

  NoArg plaq_arg;

  HmdArg hmd_arg;
  hmd_arg.n_frm_masses = 1;
  hmd_arg.frm_mass[0] = 0.1;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.max_num_iter[0] = 5;
  hmd_arg.stop_rsd[0] = 1.0E-12;
  hmd_arg.step_size = 0.02;
  hmd_arg.steps_per_traj = 20;
  hmd_arg.metropolis = METROPOLIS_YES;
  hmd_arg.reunitarize = REUNITARIZE_YES;

  //----------------------------------------------------------------
  // Run HMC Phi for Gwilson
  //----------------------------------------------------------------
  {
    GwilsonFnone lat;
    AlgHmcPhi hmc(lat,&common_arg_hmc,&hmd_arg);
    AlgPlaq plaq(lat,&common_arg_plaq,&plaq_arg);
    
    plaq.run();
    for (int i = 0; i < 1; ++i) {
      hmc.run();
      plaq.run();
    }
  }

  //----------------------------------------------------------------
  // Run HMC Phi for GpowerPlaq
  //----------------------------------------------------------------
  {
    GpowerPlaqFnone lat;
    AlgHmcPhi hmc(lat,&common_arg_hmc,&hmd_arg);
    AlgPlaq plaq(lat,&common_arg_p_plaq,&plaq_arg);
    
    plaq.run();
    for (int i = 0; i < 1; ++i) {
      hmc.run();
      plaq.run();
    }
  }

  //----------------------------------------------------------------
  // Run HMC Phi for GimprRect
  //----------------------------------------------------------------
  {
    GimprRectFnone lat;
    AlgHmcPhi hmc(lat,&common_arg_hmc,&hmd_arg);
    AlgPlaq plaq(lat,&common_arg_rn_plaq,&plaq_arg);
    AlgRect rect(lat,&common_arg_rn_rect,&plaq_arg);
    
    plaq.run();
    rect.run();
    for (int i = 0; i < 1; ++i) {
      hmc.run();
      plaq.run();
      rect.run();
    }
  }

  //----------------------------------------------------------------
  // Run HMC Phi for GpowerRect
  //----------------------------------------------------------------
  {
    GpowerRectFnone lat;
    AlgHmcPhi hmc(lat,&common_arg_hmc,&hmd_arg);
    AlgPlaq plaq(lat,&common_arg_pr_plaq,&plaq_arg);
    AlgRect rect(lat,&common_arg_pr_rect,&plaq_arg);
    
    plaq.run();
    rect.run();
    for (int i = 0; i < 1; ++i) {
      hmc.run();
      plaq.run();
      rect.run();
    }
  }

  return(0);
}






