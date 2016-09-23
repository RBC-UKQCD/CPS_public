#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/tests/various/task_func/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// test.C
//
// Test program that does the appropriate initializations
// and then calls tasks.
//
//------------------------------------------------------------------
CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/do_arg.h>
#include<alg/pbp_arg.h>
#include<alg/hmd_arg.h>
#include<alg/common_arg.h>
#include<task_func.h>
CPS_START_NAMESPACE


GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;


main(int argc,char *argv[])
{
  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 2;
  do_arg.s_node_sites = 2;
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 2;
  do_arg.z_nodes = 2;
  do_arg.t_nodes = 2;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.colors = 3;
  do_arg.beta = 6.0;
  do_arg.verbose_level = -10050402; // = 100;

  GJP.Initialize(do_arg);

  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------
  VRB.Level(GJP.VerboseLevel());

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  PbpArg pbp_arg;
  HmdArg hmd_arg;

  pbp_arg.cg.mass = 0.1;
  pbp_arg.cg.stop_rsd = 1.0E-12;
  pbp_arg.cg.max_num_iter = 500;

  hmd_arg.n_frm_masses = 2;
  hmd_arg.frm_mass[0] = 0.1;
  hmd_arg.frm_mass[1] = 0.2;
  hmd_arg.n_bsn_masses = 2;
  hmd_arg.bsn_mass[0] = 0.7;
  hmd_arg.bsn_mass[1] = 0.8;
  hmd_arg.max_num_iter[0] = 5;
  hmd_arg.max_num_iter[1] = 5;
  hmd_arg.stop_rsd[0] = 1.0E-12;
  hmd_arg.stop_rsd[1] = 1.0E-12;
  hmd_arg.step_size = 1.0E-03;
  hmd_arg.steps_per_traj = 1;

  //----------------------------------------------------------------
  // Run pbp
  //----------------------------------------------------------------
  task_pbp_g_wilson_f_stag(&common_arg, &pbp_arg);

  //----------------------------------------------------------------
  // Run HMC Phi
  //----------------------------------------------------------------
  task_hmc_phi_g_wilson_f_stag(&common_arg, &hmd_arg);

}


CPS_END_NAMESPACE
