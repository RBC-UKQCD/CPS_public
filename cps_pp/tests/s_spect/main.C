#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-07-28 13:56:38 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/s_spect/main.C,v 1.2 2003-07-28 13:56:38 mcneile Exp $
//  $Id: main.C,v 1.2 2003-07-28 13:56:38 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.1.1  2003/06/22 13:34:47  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.14  2003/01/30 14:30:53  mcneile
//  I have made two minor changes so that the code runs on a scalar
//  machine.
//
//  Revision 1.13  2002/12/04 17:16:27  zs
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
//  Revision 1.12  2002/03/11 22:26:55  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.9.2.1  2002/03/08 16:36:25  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.9  2001/09/06 11:51:38  anj
//  Minor modifications to the test suite, e.g. standardizing the
//  verbosity and such.  Collected the output from the original and the
//  latest versions using the new test suite, and checked them. Anj
//
//  Revision 1.8  2001/08/17 20:03:37  anj
//  Multiple (extra) changes to make the test suite smaller (16CPUs
//  required, not 64) and faster.  Anj
//
//  Revision 1.7  2001/08/16 12:54:19  anj
//  Some fixes follosin the float-> float change, mostly of the (variable
//  anme) float_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.6  2001/08/16 10:50:06  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "float".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and float).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.5  2001/07/03 17:00:57  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.4  2001/06/28 14:34:14  anj
//
//  The core ANSIfication should now be complete.  There are a few
//  remaining issues, but this version should compile anywhere and be
//  backward compatable with QCDSP (although this requires the top source
//  directory (.../phys/ to be added to the include path).
//
//  The serial GCC version has also been tested, and all test programs
//  appear to behave as they should (not to imply that they all work, but
//  I believe those that should work are ok).  There are minor differences
//  in the results due to rounding, (see example pbp_gccsun.dat files),
//  but that is all.
//
//  Anj.
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
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/s_spect/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  test for alg
 */

#include<util/random.h>


#include <stdio.h>
#include <stdlib.h>	// exit()
#include<config.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_hmd.h>
#include<alg/alg_s_spect.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/cg_arg.h>
#include<alg/hmd_arg.h>
#include<alg/s_spect_arg.h>
#include<alg/fix_gauge_arg.h>
#include<alg/alg_fix_gauge.h>
#include<alg/aots_s.h>
#include<util/vector.h>
#include<util/random.h>



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

#ifdef PARALLEL
  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 4;
  do_arg.t_node_sites = 6;
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 2;
  do_arg.z_nodes = 2;
  do_arg.t_nodes = 2;
  do_arg.s_nodes = 1;
#else
  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 4;
  do_arg.t_node_sites = 6;
  do_arg.s_node_sites = 0;
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
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.colors = 3;
  do_arg.beta = 5.375;
  do_arg.verbose_level = DEFAULT_VERBOSE_LEVEL;

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
  hmd_arg.frm_mass[0] = 0.1;
  hmd_arg.frm_mass[1] = 0.2;
  hmd_arg.frm_flavors[0] = 2;
  hmd_arg.frm_flavors[1] = 4;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.bsn_mass[0] = 0.8;
  hmd_arg.bsn_mass[1] = 0.9;
  hmd_arg.max_num_iter[0] = 3000;
  hmd_arg.max_num_iter[1] = 3000;
  hmd_arg.stop_rsd[0] = 1.0E-8;
  hmd_arg.stop_rsd[1] = 1.0E-12;
  hmd_arg.step_size = 0.02;
  hmd_arg.steps_per_traj = 1;
  hmd_arg.metropolis = METROPOLIS_YES;

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------

  // Fix gauge arguments
  //----------------------------------------------------------------
  FixGaugeArg farg;
  farg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
  farg.hyperplane_start = 0;
  farg.hyperplane_step = 1;
  //  farg.hyperplane_num = 24;
  // --mcneile the scalar code did not like 24
  farg.hyperplane_num = do_arg.t_node_sites ;
  farg.stop_cond = 1.E-8;
  farg.max_iter_num = 10000;

  // Staggered quark propagator arguments (point source)
  //----------------------------------------------------------------
  StagQuarkArg p_qarg1;
  p_qarg1.qid = 1;
  p_qarg1.cg.mass = 0.1;
  p_qarg1.cg.max_num_iter = 5000;
  p_qarg1.cg.stop_rsd = 1e-8;
  p_qarg1.src.type = POINT;
  p_qarg1.src.origin[0] = 0;
  p_qarg1.src.origin[1] = 0;
  p_qarg1.src.origin[2] = 0;
  p_qarg1.src.origin[3] = 0;
  p_qarg1.src.end[0] = 0;
  p_qarg1.src.end[1] = 0;
  p_qarg1.src.end[2] = 0;
  p_qarg1.src.end[3] = 0;
  p_qarg1.src.dir = HDM_T;


  // Staggered quark propagator arguments (wall2z source)
  //----------------------------------------------------------------
  StagQuarkArg wz_qarg1;
  wz_qarg1.qid = 3;
  wz_qarg1.cg.mass = 0.25;
  wz_qarg1.cg.max_num_iter = 5000;
  wz_qarg1.cg.stop_rsd = 1e-8;
  wz_qarg1.src.type = WALL2Z;
  wz_qarg1.src.origin[0] = 0;
  wz_qarg1.src.origin[1] = 0;
  wz_qarg1.src.origin[2] = 0;
  wz_qarg1.src.origin[3] = 0;
#ifdef PARALLEL
  wz_qarg1.src.end[0] = 7;
  wz_qarg1.src.end[1] = 7;
  wz_qarg1.src.end[2] = 7;
  wz_qarg1.src.end[3] = 0;
#else
  wz_qarg1.src.end[0] = 1;
  wz_qarg1.src.end[1] = 1;
  wz_qarg1.src.end[2] = 1;
  wz_qarg1.src.end[3] = 0;
#endif
  wz_qarg1.src.dir = HDM_T;
  wz_qarg1.sln = LOCAL;
  
  // Staggered quark propagator arguments (WALLZ source)
  //----------------------------------------------------------------
  StagQuarkArg wz_qarg2;
  wz_qarg2.qid = 4;
  wz_qarg2.cg.mass = 0.25;
  wz_qarg2.cg.max_num_iter = 5000;
  wz_qarg2.cg.stop_rsd = 1e-8;
  wz_qarg2.src.type = WALLZ;
  wz_qarg2.src.origin[0] = 0;
  wz_qarg2.src.origin[1] = 0;
  wz_qarg2.src.origin[2] = 0;
  wz_qarg2.src.origin[3] = 0;
#ifdef PARALLEL
  wz_qarg2.src.end[0] = 7;
  wz_qarg2.src.end[1] = 7;
  wz_qarg2.src.end[2] = 7;
  wz_qarg2.src.end[3] = 0;
#else
  wz_qarg2.src.end[0] = 1;
  wz_qarg2.src.end[1] = 1;
  wz_qarg2.src.end[2] = 1;
  wz_qarg2.src.end[3] = 0;
#endif
  wz_qarg2.src.dir = HDM_T;
  wz_qarg2.sln = NONLOCAL;


  //----------------------------------------------------------------
  // Staggered meson propagator arguments
  //----------------------------------------------------------------
  StagMesonArg marg1;
  marg1.qid0 = 1;
  marg1.qid1 = 1;
  marg1.dir = HDM_T;

  StagMesonArg marg3;
  marg3.qid0 = 3;
  marg3.qid1 = 3;
  marg3.dir = HDM_T;

  StagMesonArg marg4;
  marg4.qid0 = 4;
  marg4.qid1 = 4;
  marg4.dir = HDM_T;

  // Staggered nucleon propagator arguments
  //----------------------------------------------------------------
  StagNucleonArg narg1;
  narg1.qid0 = 1;
  narg1.qid1 = 1;
  narg1.qid2 = 1;
  narg1.dir = HDM_T;

  StagNucleonArg narg3;
  narg3.qid0 = 3;
  narg3.qid1 = 3;
  narg3.qid2 = 3;
  narg3.dir = HDM_T;

  StagNucleonArg narg4;
  narg4.qid0 = 4;
  narg4.qid1 = 4;
  narg4.qid2 = 4;
  narg4.dir = HDM_T;
    
  StagNonLocalArg darg4;
  darg4.qid0 = 4;
  darg4.qid1 = 4;
  darg4.qid2 = 4;
  darg4.dir = HDM_T;
    
  //----------------------------------------------------------------
  // Run Staggered type spectroscopy                               
  //----------------------------------------------------------------
  {
    GwilsonFstag lat;
    {
      int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
			GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();

      printf("main: Total sites :    %d\n", total_sites); 


      //-----------------------------------------------------------------
      // warming up 
      //-----------------------------------------------------------------
      {
	AlgHmdR hmd_r(lat,&common_arg,&hmd_arg);
        for (int n = 0 ; n < 0; n++) {
          printf("warmup iteration # = %d\n", n);
       	  hmd_r.run();
        }
      }
      int stride = 2;
      for (int i = 0; i < 4; i += stride ) {

        printf("iteration # = %d\n", i);

	//----------------------------------------------------------------
	// calculate action and write it to file
	//----------------------------------------------------------------

	Float sum_plaq0 = lat.SumReTrPlaq();
	Float aver_plaq0 = sum_plaq0/(18.0*total_sites);
        FILE *action;
        if( 0 == (action = fopen("action.dat", "a")) ) {
          printf("cannot open file action.dat\n");
          exit (1);
        }
	fprintf(action, "%e\n", (float)aver_plaq0);
        fclose(action);

        //------------------------------------------------------------
	// Run stag spectrocopy with point source
	//------------------------------------------------------------

        printf("main: spectroscopy staggered with Point Source\n");


        printf("Aots starts....\n");
        for(Aots it(0, 1, 2); it; ++it) {
          {
	    common_arg.results = CAST_AWAY_CONST("quark1.dat");
            AlgStagQuark sq(lat, &common_arg, &p_qarg1, it);
            sq.run();
          }
          {
	    common_arg.results = CAST_AWAY_CONST("mes11.dat");
            AlgStagMeson meson(lat, &common_arg, &marg1, it);
            meson.run();
          }
          {
	    common_arg.results = CAST_AWAY_CONST("nuc111.dat");
            AlgStagNucleon nucleon(lat, &common_arg, &narg1, it);
            nucleon.run();
          }
          {
	    common_arg.results = CAST_AWAY_CONST("quark1.dat");
            AlgStagQuark sq(lat, &common_arg, &p_qarg1, it);
            sq.free();
          }
        }


        //------------------------------------------------------------
	// fix lattice to Coulomb gauge
	//------------------------------------------------------------
	printf("AlgFixGauge.run() starts....\n");
	{ AlgFixGauge fg(lat, &common_arg, &farg);
	  fg.run();
	}

	//------------------------------------------------------------
	// Run staggered 2z wall source spectroscopy
	//------------------------------------------------------------
	printf("main: spectroscopy staggered with 2ZWall Source\n");

	printf("Aots starts....\n");
	// --mcneile failed on a scalar node
	//	for(Aots it2(0, 23, 6); it2; ++it2) {


	for(Aots it2(0, 1, 2); it2; ++it2) {
	  {
	    common_arg.results = CAST_AWAY_CONST("quark3.dat");
	    AlgStagQuark sq(lat, &common_arg, &wz_qarg1, it2);
	    sq.run();
	  }

	  {
	    common_arg.results = CAST_AWAY_CONST("mes33.dat");
	    AlgStagMeson meson(lat, &common_arg, &marg3, it2);
	    meson.run();
	  }
	  {
	    common_arg.results = CAST_AWAY_CONST("nuc333.dat");
	    AlgStagNucleon nucleon(lat, &common_arg, &narg3, it2);
	    nucleon.run();
	  }
	  {
	    common_arg.results = CAST_AWAY_CONST("quark3.dat");
	    AlgStagQuark sq(lat, &common_arg, &wz_qarg1, it2);
	    sq.free();
	  }

	}	

	//------------------------------------------------------------
	// Run staggered zwall source spectroscopy
	//------------------------------------------------------------
	printf("main: spectroscopy staggered with ZWall Source\n");

	printf("Aots starts....\n");

	for(Aots it3(0, 1 , 2); it3; ++it3) {
	  {
	    common_arg.results = CAST_AWAY_CONST("quark4.dat");
	    AlgStagQuark sq(lat, &common_arg, &wz_qarg2, it3);
	    sq.run();
	  }
	  {
	    common_arg.results = CAST_AWAY_CONST("mes44.dat");
	    AlgStagMeson meson(lat, &common_arg, &marg4, it3);
	    meson.run();
	  }
	  {
	    common_arg.results = CAST_AWAY_CONST("nuc444.dat");
	    AlgStagNucleon nucleon(lat, &common_arg, &narg4, it3);
	    nucleon.run();
	  }
	  {
	    common_arg.results = CAST_AWAY_CONST("nonlocal444.dat");
	    AlgStagNonLocal nonlocal(lat, &common_arg, &darg4, it3);
	    nonlocal.run();
	  }
	  {
	    common_arg.results = CAST_AWAY_CONST("quark4.dat");
	    AlgStagQuark sq(lat, &common_arg, &wz_qarg2, it3);
	    sq.free();
	  }
	}

	//------------------------------------------------------------
	// free Fix gauge matrices
	//------------------------------------------------------------
	printf("AlgFixGauge.free() starts...\n");
	{ AlgFixGauge fg(lat, &common_arg, &farg);
	  fg.free();
	}

	//------------------------------------------------------------
	// Run staggered HMC R_alg
	//------------------------------------------------------------

	printf("AlgHmdR starts....\n");
	{
	  AlgHmdR hmd_r(lat,&common_arg,&hmd_arg);
          for (int n = 0 ; n < 0; n++)
       	    hmd_r.run();
	}
      }

    }
  }

  return(0);
}

