#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-07-29 10:07:04 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_wilson_eig/main.C,v 1.2 2003-07-29 10:07:04 mcneile Exp $
//  $Id: main.C,v 1.2 2003-07-29 10:07:04 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.1.1  2003/06/22 13:34:49  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.11  2002/12/04 17:16:27  zs
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
//  Revision 1.10  2002/03/11 22:26:52  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.7.2.1  2002/03/08 16:36:20  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.7  2001/09/06 11:50:57  anj
//  Minor modifications to the test suite, e.g. standardizing the
//  verbosity and such.  Collected the output from the original and the
//  latest versions using the new test suite, and checked them. Anj
//
//  Revision 1.6  2001/08/17 20:03:34  anj
//  Multiple (extra) changes to make the test suite smaller (16CPUs
//  required, not 64) and faster.  Anj
//
//  Revision 1.5  2001/08/16 12:54:16  anj
//  Some fixes follosin the float-> float change, mostly of the (variable
//  anme) float_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.4  2001/07/03 17:00:56  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:13  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:12:22  anj
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
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_wilson_eig/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>	// exit()
#include<config.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_ghb.h>
#include<alg/alg_plaq.h>
#include<alg/alg_eig.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/ghb_arg.h>
#include<alg/eig_arg.h>
#include<alg/no_arg.h>

namespace cps
{
GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;
}

using namespace cps ;


int main(int argc,char *argv[])
{

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 2;
  do_arg.s_node_sites = 1;

#ifdef PARALLEL
//  do_arg.x_nodes = 2;
//  do_arg.y_nodes = 2;
//  do_arg.z_nodes = 2;
//  do_arg.t_nodes = 2;
  do_arg.x_nodes = 2; //SizeX();
  do_arg.y_nodes = 2; //SizeY();
  do_arg.z_nodes = 2; //SizeZ();
  do_arg.t_nodes = 2; //SizeT();
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
  do_arg.t_bc = BND_CND_PRD;
//  do_arg.t_bc = BND_CND_APRD;
//  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_conf_kind = START_CONF_ORD;
//  do_arg.start_conf_kind = START_CONF_LOAD;
//  do_arg.start_conf_load_addr = (Matrix *)0x5f700;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.colors = 3;
  do_arg.beta = 5.3;
  do_arg.dwf_height = 0.9;
  //do_arg.verbose_level = DEFAULT_VERBOSE_LEVEL;
  do_arg.verbose_level = 10;


  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------
  VRB.Level(GJP.VerboseLevel());


  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg_ghb;
  CommonArg common_arg_eig;
  CommonArg common_arg_plaq;

  GhbArg ghb_arg;
  ghb_arg.num_iter = 3;

  NoArg plaq_arg;

  EigArg eig_arg;
  eig_arg.N_eig = 2;
//  eig_arg.Mass_init = eig_arg.Mass_final = -1.5;
//  eig_arg.Mass_init = eig_arg.Mass_final = 0.0;
//  eig_arg.Mass_init = 0.0;
//  eig_arg.Mass_final = -2.01;
  eig_arg.Mass_init = -1.0;
//  eig_arg.Mass_final = -2.01;
  eig_arg.Mass_final = -1.051;
  eig_arg.Mass_step = -0.05;
  eig_arg.RsdR_a = 1.0E-5;
  eig_arg.RsdR_r = 1.0E-4;
  eig_arg.Rsdlam = 1.0E-5;
  eig_arg.Kalk_Sim = 1;
  eig_arg.Cv_fact = 0.1;
  eig_arg.N_min = 10;
  eig_arg.N_max = 200;
  eig_arg.N_KS_max = 50;
  eig_arg.n_renorm = 15;
  eig_arg.MaxCG = 500;
  eig_arg.ProjApsiP = 1;
  eig_arg.RitzMatOper = MAT_HERM;
  eig_arg.print_hsum = 1;
  eig_arg.hsum_dir = 3;

  //----------------------------------------------------------------
  // Run AlgEig
  //----------------------------------------------------------------
  {
    GwilsonFwilson lat;

    int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
      GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
    VRB.Result(" ", "main", "Total sites :    %d\n", total_sites); 
    VRB.Result(" ", "main", " Xn=%d Yn=%d Zn=%d Tn=%d\n", 
	       GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes());
    VRB.Result(" ", "main", " Lx=%d Ly=%d Lz=%d Lt=%d\n", 
	       GJP.Xnodes()*GJP.XnodeSites(), GJP.Ynodes()*GJP.YnodeSites(), 
	       GJP.Znodes()*GJP.ZnodeSites(), GJP.Tnodes()*GJP.TnodeSites());
    VRB.Result(" ", "main", " Bx=%d By=%d Bz=%d Bt=%d\n", 
	       GJP.Xbc(), GJP.Ybc(), GJP.Zbc(), GJP.Tbc());
    VRB.Result(" ", "main", " Bxn=%d Byn=%d Bzn=%d Btn=%d\n", 
	       GJP.XnodeBc(), GJP.YnodeBc(), GJP.ZnodeBc(), GJP.TnodeBc());

    {
      Matrix *gauge = lat.GaugeField();
      VRB.Result(" ","main","gauge = %x\n",gauge);
    }

#if 1
    {

      int ITERATIONS = 10;

      common_arg_ghb.results = CAST_AWAY_CONST("ghb.dat");
      AlgGheatBath ghb(lat,&common_arg_ghb,&ghb_arg);
      common_arg_plaq.results = CAST_AWAY_CONST("plaq.dat");
      AlgPlaq plaq(lat,&common_arg_plaq,&plaq_arg);

      for (int i = 0; i < ITERATIONS; ++i) 
      {
	plaq.run();
	ghb.run();
      }

      plaq.run();
    }
#else
    {
      common_arg_plaq.results = CAST_AWAY_CONST("plaq.dat");
      AlgPlaq plaq(lat,&common_arg_plaq,&plaq_arg);
	plaq.run();

    }
#endif
    {
      common_arg_eig.results = CAST_AWAY_CONST("eig.dat");
      AlgEig eig(lat,&common_arg_eig,&eig_arg);
      eig.run();
    }
  }

  return(0);
}







