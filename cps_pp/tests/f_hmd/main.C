#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-10-23 13:38:59 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_hmd/main.C,v 1.3 2003-10-23 13:38:59 zs Exp $
//  $Id: main.C,v 1.3 2003-10-23 13:38:59 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2.6.2  2003/09/22 04:18:44  zs
//  Back to how it was originally.
//
//  Revision 1.2.6.1  2003/09/04 14:32:35  zs
//  Temporary changes for asqtad debugging.
//
//  Revision 1.2  2003/07/29 09:57:15  mcneile
//  I have added namespace support to the codes.
//
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
//  Revision 1.10  2001/09/06 11:50:56  anj
//  Minor modifications to the test suite, e.g. standardizing the
//  verbosity and such.  Collected the output from the original and the
//  latest versions using the new test suite, and checked them. Anj
//
//  Revision 1.9  2001/08/17 23:23:21  anj
//
//  Cut down on the verbosity, but also no RNG seed output. Anj
//
//  Revision 1.8  2001/08/17 20:03:33  anj
//  Multiple (extra) changes to make the test suite smaller (16CPUs
//  required, not 64) and faster.  Anj
//
//  Revision 1.7  2001/08/17 19:38:51  anj
//  Minor alterations to the test suite to make it more practical.
//
//  Revision 1.6  2001/08/16 12:54:15  anj
//  Some fixes follosin the float-> float change, mostly of the (variable
//  anme) float_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.5  2001/08/16 10:50:05  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "float".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and float).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.4  2001/07/03 17:00:55  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:11  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:12:21  anj
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
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_hmd/main.C,v $
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
#include<alg/alg_hmd.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/cg_arg.h>
#include<alg/hmd_arg.h>
#include<alg/ghb_arg.h>

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
  FILE *fp;

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 2;
  do_arg.s_node_sites = 6;
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
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.colors = 3;
  do_arg.beta = 5.5;
  do_arg.dwf_height = 0.9;
  do_arg.clover_coeff = 2.0171;
  do_arg.verbose_level = DEFAULT_VERBOSE_LEVEL;

  // asqtad stuff
  do_arg.asqtad_KS = 1.0;
  do_arg.asqtad_naik = 0.0;
  do_arg.asqtad_lepage = 0.0;
  do_arg.asqtad_3staple = 0.0;
  do_arg.asqtad_5staple = 0.0;
  do_arg.asqtad_7staple = 0.0;
  
  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------
  VRB.Level(GJP.VerboseLevel());


  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  HmdArg hmd_arg;

  hmd_arg.n_frm_masses = 1;
  {
    Float kappa = 0.14;
    hmd_arg.frm_mass[0] = (0.5/kappa) - 4.0;
  }
  hmd_arg.frm_flavors[0] = 2;  
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.max_num_iter[0] = 500;
  hmd_arg.stop_rsd[0] = 1.0E-6;
  hmd_arg.step_size = 0.02;
  hmd_arg.steps_per_traj = 25;
  hmd_arg.metropolis = METROPOLIS_NO;
  hmd_arg.reunitarize = REUNITARIZE_YES;


  //----------------------------------------------------------------
  // Run HMC Phi Wilson
  //----------------------------------------------------------------
  {
    GwilsonFwilson lat;

    {
      common_arg.results = CAST_AWAY_CONST("hmc.dat");
      AlgHmcPhi hmc_phi(lat,&common_arg,&hmd_arg);
      
      int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
	GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
      VRB.Result(" ", "main", "Total sites :    %d\n", total_sites); 
      
      Float sum_plaq0 = lat.SumReTrPlaq();
      Float aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);

      if( (fp = fopen("plaq.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "plaq.dat");
      }
      fprintf(fp, "%e\n", float(aver_plaq0));
      fclose(fp);

      for (int i = 0; i < 1; ++i) {
	hmc_phi.run();
	sum_plaq0 = lat.SumReTrPlaq();
	aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);
	if( (fp = fopen("plaq.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "plaq.dat");
	}
	fprintf(fp,"%e\n", float(aver_plaq0));
	fclose(fp);
      }
    }
  }



  //----------------------------------------------------------------
  // Run HMD R Staggered
  //----------------------------------------------------------------
  {
    GwilsonFstag lat;
    
    {
      common_arg.results = CAST_AWAY_CONST("hmd.dat");
      AlgHmdR hmd_r(lat,&common_arg,&hmd_arg);
      
      int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
	GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
      VRB.Result(" ", "main", "Total sites :    %d\n", total_sites); 
      
      Float sum_plaq0 = lat.SumReTrPlaq();
      Float aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);

      if( (fp = fopen("plaq.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "plaq.dat");
      }
      fprintf(fp, "%e\n", float(aver_plaq0));
      fclose(fp);

      for (int i = 0; i < 1; ++i) {
	hmd_r.run();
	sum_plaq0 = lat.SumReTrPlaq();
	aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);
	if( (fp = fopen("plaq.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "plaq.dat");
	}
	fprintf(fp,"%e\n", float(aver_plaq0));
	fclose(fp);
      }
    }
  }


  //----------------------------------------------------------------
  // Run HMC Phi DWF
  //----------------------------------------------------------------
  {
    GwilsonFdwf lat;

    hmd_arg.n_bsn_masses = 1;
    hmd_arg.bsn_mass[0] = 1.0;
    hmd_arg.frm_mass[0] = 0.05;

    {
      common_arg.results = CAST_AWAY_CONST("hmc.dat");
      AlgHmcPhi hmc_phi(lat,&common_arg,&hmd_arg);
      
      int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
	GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
      VRB.Result(" ", "main", "Total sites :    %d\n", total_sites); 
      
      Float sum_plaq0 = lat.SumReTrPlaq();
      Float aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);

      if( (fp = fopen("plaq.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "plaq.dat");
      }
      fprintf(fp, "%e\n", float(aver_plaq0));
      fclose(fp);

      for (int i = 0; i < 1; ++i) {
	hmc_phi.run();
	sum_plaq0 = lat.SumReTrPlaq();
	aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);
	if( (fp = fopen("plaq.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "plaq.dat");
	}
	fprintf(fp,"%e\n", float(aver_plaq0));
	fclose(fp);
      }
    }
  }


  //----------------------------------------------------------------
  // Run HMC Phi Clover
  //----------------------------------------------------------------
  {
    GwilsonFclover lat;

    hmd_arg.n_bsn_masses = 0;
    {
      Float kappa = 0.14;
      hmd_arg.frm_mass[0] = (0.5/kappa) - 4.0;
    }

    {
      common_arg.results = CAST_AWAY_CONST("hmc.dat");
      AlgHmcPhi hmc_phi(lat,&common_arg,&hmd_arg);
      
      int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
	GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
      VRB.Result(" ", "main", "Total sites :    %d\n", total_sites); 
      
      Float sum_plaq0 = lat.SumReTrPlaq();
      Float aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);

      if( (fp = fopen("plaq.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "plaq.dat");
      }
      fprintf(fp, "%e\n", float(aver_plaq0));
      fclose(fp);

      for (int i = 0; i < 1; ++i) {
	hmc_phi.run();
	sum_plaq0 = lat.SumReTrPlaq();
	aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);
	if( (fp = fopen("plaq.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "plaq.dat");
	}
	fprintf(fp,"%e\n", float(aver_plaq0));
	fclose(fp);
      }
    }
  }

  
  return(1);  
}





  


