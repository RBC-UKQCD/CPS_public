#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:12 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/glb_sum/main.C,v 1.5 2004-08-18 11:58:12 zs Exp $
//  $Id: main.C,v 1.5 2004-08-18 11:58:12 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/glb_sum/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <stdlib.h>	// exit()
#include<config.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<comms/glb.h>
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
  do_arg.s_node_sites = 0;

#ifdef PARALLEL
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 2;
  do_arg.z_nodes = 2;
  do_arg.t_nodes = 2;
#else
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
#endif

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_seed_kind = START_SEED_FIXED;

  do_arg.beta = 6.0;
  do_arg.dwf_height = 0.9;


  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------

  VRB.DeactivateLevel(VERBOSE_RNGSEED_LEVEL);
  VRB.ActivateLevel(VERBOSE_FUNC_LEVEL);
  

  //----------------------------------------------------------------
  // Test global sum
  //----------------------------------------------------------------
  Float sum = 1.1;

  VRB.Debug("sum = %g\n", sum);
  glb_sum(&sum);
  VRB.Debug("after global sum sum = %g\n", sum);

  glb_sum(&sum);
  VRB.Debug("after global sum sum = %g\n", sum);

  glb_sum(&sum);
  VRB.Debug("after global sum sum = %g\n", sum);

  glb_sum(&sum);
  VRB.Debug("after global sum sum = %g\n", sum);



}






CPS_END_NAMESPACE
