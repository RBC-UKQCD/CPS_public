#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-09-21 20:16:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_stag_pbp/main.C,v 1.9 2004-09-21 20:16:52 chulwoo Exp $
//  $Id: main.C,v 1.9 2004-09-21 20:16:52 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.9 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_stag_pbp/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#include<util/lattice.h>
#include<alg/alg_pbp.h>
#include<alg/do_arg.h>




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
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
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
  do_arg.beta = 6.0;

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
  pbp_arg.max_num_iter = 500;

  //----------------------------------------------------------------
  // Run AlgPbp
  //----------------------------------------------------------------
  {
    GwilsonFstag lat;
    common_arg.set_filename("pbp.dat");
    AlgPbp pbp(lat,&common_arg,&pbp_arg);
    int num_hits;
    int i;

    num_hits = 1;
    for(i=0; i<num_hits; i++){
      pbp.run();
    }

  }

  return 1;
}







