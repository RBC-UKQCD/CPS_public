#include<config.h>

//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004/09/21 20:16:54 $
//  $Header: /space/cvs/cps/cps++/tests/f_wilson_pbp/main.C,v 1.7 2004/09/21 20:16:54 chulwoo Exp $
//  $Id: main.C,v 1.7 2004/09/21 20:16:54 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.7 $
//  $Source: /space/cvs/cps/cps++/tests/f_wilson_pbp/main.C,v $
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



USING_NAMESPACE_CPS

int main(int argc,char *argv[])
{

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



  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  common_arg.results = CAST_AWAY_CONST("pbp.dat");

  PbpArg pbp_arg;
  pbp_arg.pattern_kind = ARRAY;
  pbp_arg.n_masses = 1;
  pbp_arg.mass[0] = 0.1;
  pbp_arg.stop_rsd = 1.0E-12;
  pbp_arg.max_num_iter = 500;

  //----------------------------------------------------------------
  // Run AlgPbp
  //----------------------------------------------------------------
  {
    GwilsonFwilson lat;
    common_arg.results = CAST_AWAY_CONST("pbp.dat");
    AlgPbp pbp(lat,&common_arg,&pbp_arg);
    int num_hits;
    int i;

    num_hits = 1;
    for(i=0; i<num_hits; i++){
      pbp.run();
    }

  }

  return(0);
}







