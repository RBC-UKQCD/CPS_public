#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:17 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/g_hb/main.C,v 1.5 2004-06-04 21:14:17 chulwoo Exp $
//  $Id: main.C,v 1.5 2004-06-04 21:14:17 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/g_hb/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#include<util/lattice.h>
#include<alg/do_arg.h>
#include<alg/alg_ghb.h>
#include<alg/no_arg.h>
#include<alg/alg_plaq.h>

CPS_START_NAMESPACE
GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;
CPS_END_NAMESPACE

USING_NAMESPACE_CPS


int main(int argc,char *argv[])
{
  int ITERATIONS = 1;

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  //-----------------------------------------------------
  // Timing results:
  // 20 it, 16^3x4 on 1mb 75sec (w/asm)
  // 20 it, 16^3x4 on 1mb 45sec (w/cram w/asm)
  // 20 it, 16^3x4 on 1mb 37sec (w/o cmhb kernal)
  // 20 it, 16^3x4 on 1mb  8sec (w/o update - action only)
  do_arg.t_node_sites = 2;
  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.s_node_sites = 2;

#ifdef PARALLEL
  do_arg.t_nodes = 2;
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 2;
  do_arg.z_nodes = 2;
  do_arg.s_nodes = 1;
#else
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
  do_arg.s_nodes = 1;
#endif 
  do_arg.t_bc = BND_CND_PRD;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_APRD;

  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 5.71;
// do_arg.beta = 1000;
  do_arg.dwf_height = 0.9;


  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg_ghb;

  CommonArg common_arg_plaq;

  GhbArg ghb_arg;
  ghb_arg.num_iter = 1;

  NoArg plaq_arg;


  //----------------------------------------------------------------
  // Run GHB
  //----------------------------------------------------------------
  {
    GwilsonFstag lat;
    
    {
      common_arg_ghb.set_filename("ghb.dat");
      AlgGheatBath ghb(lat,&common_arg_ghb,&ghb_arg);
      
      common_arg_plaq.set_filename("plaq.dat");
      AlgPlaq plaq(lat,&common_arg_plaq,&plaq_arg);

      for (int i = 0; i < ITERATIONS; ++i) {
	ghb.run();
	plaq.run();
      }
    }
  }

  return(0);
}






