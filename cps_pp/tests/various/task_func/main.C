#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/task_func/main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2002/12/04 17:16:27  zs
//  Merged the new 2^4 RNG into the code.
//  This new RNG is implemented in the LatRanGen class.
//  The following algorithm and utility classes are affected:
//
//  AlgEig                  Fdwf
//  AlgGheatBath            Fstag
//  AlgHmd                  GlobalJobParameter
//  AlgNoise                Lattice
//  AlgPbp                  Matrix
//  AlgThreept              RandomGenerator
//                          Vector
//
//  Revision 1.4  2001/08/17 20:03:39  anj
//  Multiple (extra) changes to make the test suite smaller (16CPUs
//  required, not 64) and faster.  Anj
//
//  Revision 1.3  2001/08/16 12:54:20  anj
//  Some fixes follosin the float-> float change, mostly of the (variable
//  anme) float_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.2  2001/06/19 18:12:33  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:04  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: main.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/task_func/main.C,v $
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
