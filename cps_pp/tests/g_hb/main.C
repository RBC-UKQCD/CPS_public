#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-07-29 11:03:06 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/g_hb/main.C,v 1.2 2003-07-29 11:03:06 mcneile Exp $
//  $Id: main.C,v 1.2 2003-07-29 11:03:06 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.1.1  2003/06/22 13:34:49  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
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
//  Revision 1.7  2001/09/06 11:50:58  anj
//  Minor modifications to the test suite, e.g. standardizing the
//  verbosity and such.  Collected the output from the original and the
//  latest versions using the new test suite, and checked them. Anj
//
//  Revision 1.6  2001/08/17 20:03:35  anj
//  Multiple (extra) changes to make the test suite smaller (16CPUs
//  required, not 64) and faster.  Anj
//
//  Revision 1.5  2001/08/16 12:54:17  anj
//  Some fixes follosin the float-> float change, mostly of the (variable
//  anme) float_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.4  2001/07/03 17:00:57  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:14  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:12:23  anj
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
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/g_hb/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  test for alg
 */

#include <stdio.h>
#include <stdlib.h>	// exit()
#include<config.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/cg_arg.h>
#include<alg/ghb_arg.h>
#include<alg/alg_ghb.h>
#include<alg/no_arg.h>
#include<alg/alg_plaq.h>

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
  do_arg.colors = 3;
  do_arg.beta = 5.71;
// do_arg.beta = 1000;
  do_arg.dwf_height = 0.9;
  do_arg.verbose_level = DEFAULT_VERBOSE_LEVEL;

  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------
  VRB.Level(GJP.VerboseLevel());


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
      common_arg_ghb.results = CAST_AWAY_CONST("ghb.dat");
      AlgGheatBath ghb(lat,&common_arg_ghb,&ghb_arg);
      
      common_arg_plaq.results = CAST_AWAY_CONST("plaq.dat");
      AlgPlaq plaq(lat,&common_arg_plaq,&plaq_arg);

      for (int i = 0; i < ITERATIONS; ++i) {
	ghb.run();
	plaq.run();
      }
    }
  }

  return(0);
}






