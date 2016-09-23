#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/02/08 18:35:09 $
//  $Header: /space/cvs/cps/cps++/tests/various/ghb/main.C,v 1.7 2008/02/08 18:35:09 chulwoo Exp $
//  $Id: main.C,v 1.7 2008/02/08 18:35:09 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.7 $
//  $Source: /space/cvs/cps/cps++/tests/various/ghb/main.C,v $
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
#include<alg/alg_plaq.h>
#include<alg/alg_ghb.h>
#include<alg/do_arg.h>
#include<alg/hmd_arg.h>
#include<alg/common_arg.h>
#include<alg/no_arg.h>
CPS_START_NAMESPACE

#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif


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
#else
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
#endif 

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_seed_kind = START_SEED_FIXED;

  do_arg.beta = 4.8;
  do_arg.dwf_height = 0.9;

  GJP.Initialize(do_arg);

  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------
  VRB.DeactivateLevel(VERBOSE_RNGSEED_LEVEL);
  VRB.ActivateLevel(VERBOSE_FUNC_LEVEL);
  

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  NoArg plaq_arg;

  GhbArg ghb_arg;
    ghb_arg.num_iter = 10;

  
  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  VRB.Result(" ", "main", "Total sites :    %d\n", total_sites); 

  //----------------------------------------------------------------
  // Run gauge heat bath
  //----------------------------------------------------------------
  {
    GwilsonFnone lat;
    AlgGheatBath ghb(lat,&common_arg,&ghb_arg);
    AlgPlaq plaq(lat,&common_arg,&plaq_arg);
    
    plaq.run();
    
    ghb.run();
    plaq.run();
  }


}






CPS_END_NAMESPACE
