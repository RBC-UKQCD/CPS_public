#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-02-08 18:35:09 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/utils/getseeds.C,v 1.7 2008-02-08 18:35:09 chulwoo Exp $
//  $Id: getseeds.C,v 1.7 2008-02-08 18:35:09 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: getseeds.C,v $
//  $Revision: 1.7 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/utils/getseeds.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/random.h>
#include<alg/alg_pbp.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/pbp_arg.h>

#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif

USING_NAMESPACE_CPS

int main(int argc,char *argv[]){

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  do_arg.x_node_sites = 4;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 4;
  do_arg.s_node_sites = 0;

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

  do_arg.beta = 6.0;
  do_arg.dwf_height = 0.9;

  GJP.Initialize(do_arg);



  VRB.Deactivate_level(VERBOSE_RNGSEED_LEVEL);

  //----------------------------------------------------------------
  // Get all of the seeds:
  //----------------------------------------------------------------

#ifdef PARALLEL

  printf("[%i,%i,%i,%i] %i %i %i %i\n",
	 CoorT(),CoorX(),CoorY(),CoorZ(),
	 Seed(),SeedS(),SeedT(),SeedST()
	 );

#else

  printf("[0,0,0,0] %i\n", SERIAL_SEED );

#endif

  return(0);
}




