#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:12 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/inst/main.C,v 1.5 2004-08-18 11:58:12 zs Exp $
//  $Id: main.C,v 1.5 2004-08-18 11:58:12 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/inst/main.C,v $
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
#include<alg/alg_inst.h>
#include<alg/alg_plaq.h>
#include<alg/alg_pbp.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/inst_arg.h>
#include<alg/pbp_arg.h>
#include<alg/no_arg.h>
CPS_START_NAMESPACE


GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;

main(int argc,char *argv[])
{

  //----------------------------------------------------------------
  // Initialize run parameters
  //----------------------------------------------------------------
  int lx = 4;
  int ly = 4;
  int lz = 4;
  int lt = 4;

//  int verbose = 0;
  int verbose = -7; //Debug

  Float mass = 0.1;
  int num_hits = 0;

  int ls = 4;
  int dwf_height = 0.9;

//  InstType inst_kind = NONSINGULARSQUASHEDTRANSFORMED;
  InstType inst_kind = SINGULAR;
  Float charge = 1.0;
  Float rho = 4.0;
  Float rho_cutoff = 3.0;
  Float noise = 0.0;



  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

#ifdef PARALLEL
  do_arg.x_node_sites = lx;
  do_arg.y_node_sites = ly;
  do_arg.z_node_sites = lz;
  do_arg.t_node_sites = lt;
  do_arg.s_node_sites = ls;
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 2;
  do_arg.z_nodes = 2;
  do_arg.t_nodes = 2;
#else
  do_arg.x_node_sites = lx;
  do_arg.y_node_sites = ly;
  do_arg.z_node_sites = lz;
  do_arg.t_node_sites = lt;
  do_arg.s_node_sites = ls;
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
  do_arg.dwf_height = dwf_height;
  

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

  PbpArg pbp_arg;
  pbp_arg.cg.mass = mass;
  pbp_arg.cg.stop_rsd = 1.0E-12;
  pbp_arg.cg.max_num_iter = 5000;
  pbp_arg.src_u_s = 0;
  pbp_arg.src_l_s = GJP.SnodeSites()-1;
  pbp_arg.snk_u_s = GJP.SnodeSites()-1;
  pbp_arg.snk_l_s = 0;

  NoArg plaq_arg;

  InstArg inst_arg;
  inst_arg.inst_kind = inst_kind;
  inst_arg.charge = charge;
  inst_arg.rho = rho;
  inst_arg.rho_cutoff = rho_cutoff;
  inst_arg.noise = noise;



  //----------------------------------------------------------------
  // Run AlgInst
  //----------------------------------------------------------------
  {
    GwilsonFnone lat;
    AlgInst inst(lat,&common_arg,&inst_arg);
    common_arg.results = "plaq.dat";
    AlgPlaq plaq(lat,&common_arg,&plaq_arg);

    Float *u = (Float *) lat.GaugeField();

    inst.run();
    

    //--------------------------------------------------------------
    // Print out gauge field
    //--------------------------------------------------------------
    for(int i=0; i< 4*GJP.VolSites(); i=i+18){
      VRB.Debug("(%e, %e)  (%e, %e)  (%e, %e)\n", 
		float(u[i]), float(u[i+1]),
		float(u[i+2]), float(u[i+3]),
		float(u[i+4]), float(u[i+5]));
      VRB.Debug("(%e, %e)  (%e, %e)  (%e, %e)\n", 
		float(u[i+6]), float(u[i+7]),
		float(u[i+8]), float(u[i+9]),
		float(u[i+10]), float(u[i+11]));
      VRB.Debug("(%e, %e)  (%e, %e)  (%e, %e)\n\n\n", 
		float(u[i+12]), float(u[i+13]),
		float(u[i+14]), float(u[i+15]),
		float(u[i+16]), float(u[i+17]));
    }
    


    //--------------------------------------------------------------
    // Print out 1 - plaquette/3 field
    //--------------------------------------------------------------
    int coor[4];
    float plq;
    for(coor[0]=0; coor[0]<GJP.XnodeSites(); coor[0]++)
    for(coor[1]=0; coor[1]<GJP.YnodeSites(); coor[1]++)
    for(coor[2]=0; coor[2]<GJP.ZnodeSites(); coor[2]++)
    for(coor[3]=0; coor[3]<GJP.TnodeSites(); coor[3]++)
      for(int mu=0; mu<4; mu++)
      for(int nu=0; nu<4; nu++){
	plq = lat.ReTrPlaq(coor, mu, nu);
	VRB.Debug("plaq(%d,%d,%d,%d; %d,%d) = %e\n",
		  coor[0],
		  coor[1],
		  coor[2],
		  coor[3],
		  mu,
		  nu,
		  1.0 - (plq/3.0));
      }

    plaq.run();
  }


  //----------------------------------------------------------------
  // Run AlgPlaq
  //----------------------------------------------------------------
  {
    GwilsonFnone lat;
    common_arg.results = "plaq.dat";
    AlgPlaq plaq(lat,&common_arg,&plaq_arg);

    plaq.run();
  }

  //----------------------------------------------------------------
  // Run AlgPbp
  //----------------------------------------------------------------
  {
    GwilsonFstag lat;
    common_arg.results = "pbp.dat";
    AlgPbp pbp(lat,&common_arg,&pbp_arg);
    int i;

    for(i=0; i<num_hits; i++){
      pbp.run();
    }

  }

}






CPS_END_NAMESPACE
