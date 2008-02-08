#include<config.h>
#include<stdlib.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-02-08 18:35:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/g_hmd/main.C,v 1.14 2008-02-08 18:35:08 chulwoo Exp $
//  $Id: main.C,v 1.14 2008-02-08 18:35:08 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.14 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/g_hmd/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#include<util/gjp.h>
#include<util/lattice.h>
#include<alg/alg_hmd.h>
#include<alg/alg_plaq.h>
#include<alg/alg_rect.h>
#include<alg/do_arg.h>
#include<alg/no_arg.h>
#include<comms/sysfunc_cps.h>




USING_NAMESPACE_CPS





int main(int argc,char *argv[])
{
  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  Start();
  DoArg do_arg;

  int nx,ny,nz,nt;
  //----------------------------------------------------------------
  // tests are generated with (nx,ny,nz,nt) = (8 8 8 8)
  //----------------------------------------------------------------
  if (argc<5) {printf("usage: %s nx ny nz nt\n",argv[0]);exit(-2);}
  sscanf(argv[1],"%d",&nx);
  sscanf(argv[2],"%d",&ny);
  sscanf(argv[3],"%d",&nz);
  sscanf(argv[4],"%d",&nt);
  printf("sizes = %d %d %d %d\n",nx,ny,nz,nt);
  do_arg.x_node_sites = nx/SizeX();
  do_arg.y_node_sites = ny/SizeY();
  do_arg.z_node_sites = nz/SizeZ();
  do_arg.t_node_sites = nt/SizeT();
  do_arg.s_node_sites = 2;

  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = 1;

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
  do_arg.u0 = 0.9; // for GimprOLSym

  GJP.Initialize(do_arg);
  VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);
//  VRB.ActivateLevel(VERBOSE_FUNC_LEVEL);

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
  hmd_arg.metropolis = METROPOLIS_NO;
  hmd_arg.reunitarize = REUNITARIZE_YES;

#if 1
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
#endif

  //----------------------------------------------------------------
  // Run HMC Phi for GimprOLSym
  //----------------------------------------------------------------
  {
    GimprOLSymFnone lat;
    AlgHmcPhi hmc(lat,&common_arg_hmc,&hmd_arg);
    AlgPlaq plaq(lat,&common_arg_pr_plaq,&plaq_arg);
    AlgRect rect(lat,&common_arg_pr_rect,&plaq_arg);
    
    plaq.run();
    rect.run();
    for (int i = 0; i < 3; ++i) {
      hmc.run();
      plaq.run();
      rect.run();
    }
  }

  End();
  return(0);
}






