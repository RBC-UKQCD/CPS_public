#include<config.h>

//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:10 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_dwfso_pbp/main.C,v 1.5 2004-08-18 11:58:10 zs Exp $
//  $Id: main.C,v 1.5 2004-08-18 11:58:10 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_dwfso_pbp/main.C,v $
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

#ifdef PARALLEL
  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 4;
  do_arg.s_node_sites = 1;
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
  do_arg.s_nodes = 2;
  do_arg.s_axis = SCU_T;
#else
  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 4;
  do_arg.s_node_sites = 2;
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
  do_arg.start_conf_load_addr = (Matrix *) 0x12f13;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.start_seed_value = 1111;
  do_arg.beta = 5.8;
  do_arg.dwf_height = 0.9;

#if TARGET==cpsMPI
    MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.s_nodes);    
    using MPISCU::fprintf;
    using MPISCU::printf;
#endif

  GJP.Initialize(do_arg);



  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  common_arg.results = CAST_AWAY_CONST("pbp.dat");

  PbpArg pbp_arg;
  pbp_arg.pattern_kind = ARRAY;
  pbp_arg.n_masses = 1;
  pbp_arg.mass[0] = 0.01;
  pbp_arg.stop_rsd = 1.0E-10;
  pbp_arg.max_num_iter = 5000;
  pbp_arg.src_u_s = 0;
  pbp_arg.src_l_s = GJP.SnodeSites() * GJP.Snodes() - 1;
  pbp_arg.snk_u_s = GJP.SnodeSites() * GJP.Snodes() - 1;
  pbp_arg.snk_l_s = 0;
  pbp_arg.snk_loop = 0;

  int num_hits = 1;
  
  //------------------------------------------------------------
  // Run AlgPbp, a5 = 1.0
  //------------------------------------------------------------  
  GJP.DwfA5Inv(1.0);
  {
    GwilsonFdwf lat;
    AlgPbp pbp(lat,&common_arg,&pbp_arg);
    for(int i=0; i<num_hits; i++){
      pbp.run();
    }
  }

  //------------------------------------------------------------
  // Run AlgPbp, a5 = 0.5
  //------------------------------------------------------------  
  GJP.DwfA5Inv(2.0);
  {
    GwilsonFdwf lat;
    AlgPbp pbp(lat,&common_arg,&pbp_arg);
    for(int i=0; i<num_hits; i++){
      pbp.run();
    }
  }


  return(1);
}






