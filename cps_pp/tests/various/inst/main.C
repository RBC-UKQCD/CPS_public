#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/inst/main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.8  2002/12/04 17:16:27  zs
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
//  Revision 1.7  2001/08/17 20:03:38  anj
//  Multiple (extra) changes to make the test suite smaller (16CPUs
//  required, not 64) and faster.  Anj
//
//  Revision 1.6  2001/08/16 12:54:20  anj
//  Some fixes follosin the float-> float change, mostly of the (variable
//  anme) float_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.5  2001/08/16 10:50:07  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "float".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and float).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.4  2001/07/03 17:00:58  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:15  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:12:30  anj
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
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/inst/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <stdio.h>
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
  do_arg.colors = 3;
  do_arg.beta = 6.0;
  do_arg.dwf_height = dwf_height;
  do_arg.verbose_level = -10050402; // = verbose;

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
