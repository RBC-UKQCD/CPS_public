#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/threept/main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
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
//  Revision 1.10  2002/03/11 22:26:57  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.7.2.1  2002/03/08 16:36:28  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.7  2001/09/06 11:51:38  anj
//  Minor modifications to the test suite, e.g. standardizing the
//  verbosity and such.  Collected the output from the original and the
//  latest versions using the new test suite, and checked them. Anj
//
//  Revision 1.6  2001/08/17 20:03:37  anj
//  Multiple (extra) changes to make the test suite smaller (16CPUs
//  required, not 64) and faster.  Anj
//
//  Revision 1.5  2001/08/16 12:54:19  anj
//  Some fixes follosin the float-> float change, mostly of the (variable
//  anme) float_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.4  2001/07/03 17:00:58  anj
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
//  Revision 1.2  2001/06/19 18:12:28  anj
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
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/threept/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  test for alg
 */

CPS_END_NAMESPACE
#include <stdio.h>
#include <stdlib.h>	// exit()
#include<config.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif
CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_ghb.h>
#include<alg/alg_plaq.h>
#include<alg/alg_threept.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/ghb_arg.h>
#include<alg/threept_arg.h>
#include<alg/no_arg.h>
#include<alg/fix_gauge_arg.h>
#include<alg/alg_fix_gauge.h>
CPS_START_NAMESPACE

GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;
int main(int argc,char *argv[])
{
 
  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

#ifdef PARALLEL
// For the 100 Gflop:
  do_arg.x_node_sites = 2;
       do_arg.x_nodes = 2;
  do_arg.y_node_sites = 2;
       do_arg.y_nodes = 2;
  do_arg.z_node_sites = 2;
       do_arg.z_nodes = 2;
  do_arg.t_node_sites = 2;
       do_arg.t_nodes = 2;
  do_arg.s_node_sites = 4;
       do_arg.s_nodes = 1;


#else
  do_arg.x_node_sites = 4;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 4;
  do_arg.t_node_sites = 4;
  do_arg.s_node_sites = 4;
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
  do_arg.s_nodes = 1;
#endif

  //----------------------------------------------------------------
  //  Get seed and lattice location
  //----------------------------------------------------------------

  do_arg.start_seed_kind = START_SEED_FIXED;

  do_arg.start_conf_kind = START_CONF_DISORD;

  do_arg.colors = 3;
  do_arg.beta = 6.0;
  do_arg.dwf_height = 1.8;
  do_arg.verbose_level = DEFAULT_VERBOSE_LEVEL;
  //do_arg.verbose_level = 10;


  //----------------------------------------------------------------
  // Initialize the GJP class
  //----------------------------------------------------------------
  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------
  VRB.Level(GJP.VerboseLevel());


  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  CommonArg common_arg_ghb;
  CommonArg common_arg_plaq;
  ThreePtArg threept_arg;

  NoArg plaq_arg;

  //----------------------------------------------------------------
  // Run Wilson type spectrum                               
  //----------------------------------------------------------------

  threept_arg.cg.stop_rsd = 1.0E-6;
  threept_arg.cg.max_num_iter = 4;
  threept_arg.seed = 16;
	// obsolete
  threept_arg.t_src = 2;
  threept_arg.t_Op = 3;
  threept_arg.t_Op_2 = 5;
  threept_arg.t_sink = 6;
  threept_arg.num_masses = 1;
  threept_arg.mass[0] = .50;
	// strange quark
  threept_arg.mass[1] = .40;
  threept_arg.mass[2] = .030;
  threept_arg.mass[3] = .040;
  threept_arg.mass[4] = .050;
  threept_arg.mass[5] = .075;
  threept_arg.mass[6] = .100;
  threept_arg.mass[7] = .125;

  FixGaugeArg fix_arg;
  fix_arg.fix_gauge_kind=FIX_GAUGE_COULOMB_T;
  fix_arg.hyperplane_start=0;
  fix_arg.hyperplane_step=1;
  fix_arg.hyperplane_num=GJP.Tnodes()*GJP.TnodeSites();
  fix_arg.stop_cond=1e-3;
  fix_arg.max_iter_num=10;
  
  {

    GwilsonFdwf lat;

    common_arg.results = CAST_AWAY_CONST("threept.dat");
    FILE *fp;
    if( (fp = fopen((char *)(common_arg.results), "a")) == NULL ) {
      ERR.FileA("main", "main", (char *)(common_arg.results) );
    }
    //    fprintf(fp, "seed= %d\n",do_arg.start_seed_value);
    fprintf(fp, "do_arg.dwf_height= %e\n", do_arg.dwf_height);
    fprintf(fp, "prop stop_rsd= %e\n",threept_arg.cg.stop_rsd);
    fprintf(fp, "gauge fix type= %d\n",fix_arg.fix_gauge_kind);
    fprintf(fp, "gauge fix stop_cond= %e\n",fix_arg.stop_cond);
    fprintf(fp, "t_src= %d\n",threept_arg.t_src);
    fprintf(fp, "t_sink= %d\n",threept_arg.t_sink);
    fprintf(fp, "t_Op= %d\n",threept_arg.t_Op);
    fclose(fp);

    AlgPlaq plaq(lat,&common_arg,&plaq_arg);
    plaq.run();

    if(fix_arg.fix_gauge_kind==FIX_GAUGE_COULOMB_T ||
       fix_arg.fix_gauge_kind==FIX_GAUGE_LANDAU){
         AlgFixGauge fix_gauge(lat,&common_arg,&fix_arg);
         fix_gauge.run();
    }


    AlgThreePt threept(lat,&common_arg,&threept_arg);
    threept.run();


  }

  return(0);
}






CPS_END_NAMESPACE
