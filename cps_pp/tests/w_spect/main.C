#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-07-29 09:48:40 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/w_spect/main.C,v 1.2 2003-07-29 09:48:40 mcneile Exp $
//  $Id: main.C,v 1.2 2003-07-29 09:48:40 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.1.1  2003/06/22 13:34:47  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.12  2003/02/17 22:29:40  mcneile
//  I have turned the normal mesons on.
//
//  Revision 1.11  2003/01/31 16:23:58  mcneile
//  I have set the scaling factor to 1.
//  In the routine that created the source for the quark propagator
//  it had some weird default mechanism.
//  This was a problem on lomond (sun SMP) and florence (raid box).
//
//  Revision 1.10  2002/12/04 17:16:27  zs
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
//  Revision 1.9  2001/09/06 11:51:41  anj
//  Minor modifications to the test suite, e.g. standardizing the
//  verbosity and such.  Collected the output from the original and the
//  latest versions using the new test suite, and checked them. Anj
//
//  Revision 1.8  2001/08/17 20:03:41  anj
//  Multiple (extra) changes to make the test suite smaller (16CPUs
//  required, not 64) and faster.  Anj
//
//  Revision 1.7  2001/08/16 12:54:22  anj
//  Some fixes follosin the float-> float change, mostly of the (variable
//  anme) float_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.6  2001/08/16 10:50:09  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "float".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and float).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.5  2001/07/03 17:01:00  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.4  2001/06/28 14:34:16  anj
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
//  Revision 1.3  2001/06/21 15:40:16  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:12:35  anj
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
//  Revision 1.2  2001/05/25 06:16:05  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: main.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/w_spect/main.C,v $
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
#include<alg/alg_hmd.h>
#include<alg/alg_pbp.h>
#include<alg/alg_w_spect.h>
#include<alg/alg_plaq.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/hmd_arg.h>
#include<alg/pbp_arg.h>
#include<alg/w_spect_arg.h>
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

  FILE *fp;

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

#ifdef PARALLEL
  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 2;
  do_arg.s_node_sites = 4;
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 2;
  do_arg.z_nodes = 2;
  do_arg.t_nodes = 2;
  do_arg.s_nodes = 1;
#else
  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 2;
  do_arg.s_node_sites = 4;
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
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_conf_load_addr = (Matrix *) 0x12f13;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.start_seed_value = 1111;
  do_arg.colors = 3;
  do_arg.beta = 5.8;
  do_arg.dwf_height = 0.9;
  do_arg.verbose_level = DEFAULT_VERBOSE_LEVEL;

  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------
  VRB.Level(GJP.VerboseLevel());


  //----------------------------------------------------------------
  // Initialize arguments for the new Wilson spectrum code. -- Ping
  //----------------------------------------------------------------
  AlgWspect::SetCounter(1 /* 1st measurement */, 1 /* stride */);
  // The above is also the default behavior.

  WspectArg w_spect_arg;
  w_spect_arg.cg.mass = 0.1;
  w_spect_arg.cg.stop_rsd = 1.0E-12;
  w_spect_arg.cg.max_num_iter = 5000;
  w_spect_arg.prop_dir = 3;
  w_spect_arg.num_mom = 1;              // zero spatial momentum only
  w_spect_arg.source_kind = POINT_W; 
  w_spect_arg.aots_start = 0;
  w_spect_arg.aots_step = 1;
  w_spect_arg.aots_num = 1;
  w_spect_arg.baryons_on = 1;
  w_spect_arg.normal_mesons_on = 1;

  // this scales the quark source by 1, the default
  // mechanism seems to be broken in w_quark
  w_spect_arg.rescale_factor = 1.0 ;
  

  WspectOutput wout;
  wout.fold = BARYON_PAST;      // BARYON_FOLD, BARYON_RAW, BARYON_PAST
  wout.cg = "w_spect.dat";
  wout.pbp = "w_pbp.dat";

  wout.a0 = "a0.dat";
  wout.a0_prime = "a0_prime.dat";
  wout.a1_x = "a1_x.dat";
  wout.a1_y = "a1_y.dat";
  wout.a1_z = "a1_z.dat";
  wout.b1_x = "b1_x.dat"; 
  wout.b1_y = "b1_y.dat"; 
  wout.b1_z = "b1_z.dat";  
  wout.pion = "pion.dat";  
  wout.pion_prime = "pion_prime.dat"; 
  wout.rho_x = "rho_x.dat";  	
  wout.rho_y = "rho_y.dat";  	
  wout.rho_z = "rho_z.dat";  	 	    
  wout.rho_x_prime = "rho_x_prime.dat";  	 	      	    
  wout.rho_y_prime = "rho_y_prime.dat";  	 	     
  wout.rho_z_prime = "rho_z_prime.dat";  	 	      	     
  wout.nucleon = "nucleon.dat";  	 	      	     
  wout.nucleon_prime = "nucleon_prime.dat";  	 	      	           
  wout.delta_x = "delta_x.dat";  	 	      	      	     
  wout.delta_y = "delta_y.dat";  	 	      	      	      	     
  wout.delta_z = "delta_z.dat";  	 	      	      	      	     
  wout.delta_t = "delta_t.dat";  	 	      	      	      	     


  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;

  PbpArg pbp_arg;
  pbp_arg.pattern_kind = ARRAY;
  pbp_arg.n_masses = 1;
  pbp_arg.mass[0] = 0.1;
  pbp_arg.stop_rsd = 1.0E-12;
  pbp_arg.max_num_iter = 5000;
  pbp_arg.src_u_s = 0;
  pbp_arg.src_l_s = GJP.SnodeSites()-1;
  pbp_arg.snk_u_s = GJP.SnodeSites()-1;
  pbp_arg.snk_l_s = 0;
  pbp_arg.snk_loop = 0;


  HmdArg hmd_arg;
  hmd_arg.n_frm_masses = 1;
  hmd_arg.frm_mass[0] = 0.1;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.max_num_iter[0] = 5;
  hmd_arg.stop_rsd[0] = 1.0E-12;
  hmd_arg.step_size = 0.01;
  hmd_arg.steps_per_traj = 50;
  hmd_arg.metropolis = METROPOLIS_YES;
  hmd_arg.reunitarize = REUNITARIZE_YES;

  NoArg plaq_arg;


  //----------------------------------------------------------------
  // Initialize run parameters
  //----------------------------------------------------------------
  int num_therm_no_met = 1;
  int num_therm_met = 1;
  int pbp_num_hits = 1;
  int spect_yes = 1;
  int num_meas = 1;
  int num_sweeps = 1;
  int min_ls = 3;
  int num_ls = 1;
  int ls_step = 1;
  Float mass_start = 0.0;
  Float mass_step = -0.0;
  int num_mass = 1;


  //----------------------------------------------------------------
  // Initialize counters
  //----------------------------------------------------------------
  int iter = 0;
  int init_g_upd_cnt = 0;
  int cur_g_upd_cnt = 0;

  int total_sites = GJP.VolNodeSites() * GJP.Xnodes() *
    GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();

  if( (fp = fopen("info.dat", "a")) == NULL ) {
    ERR.FileA(" ","main", "info.dat");
  }
  fprintf(fp, "Total sites :    %d\n", total_sites); 
  fclose(fp);

  //----------------------------------------------------------------
  // Run AlgPlaq 
  //----------------------------------------------------------------
  {
    GwilsonFnone lat;
    common_arg.results = CAST_AWAY_CONST("plaq.dat");
    AlgPlaq plaq(lat,&common_arg,&plaq_arg);
    
    plaq.run();
  }


  //----------------------------------------------------------------
  // Thermalize gauge field using HMC Phi with no metropolis
  //----------------------------------------------------------------
  hmd_arg.metropolis = METROPOLIS_NO;
  if(num_therm_no_met != 0)
    {
      GwilsonFnone lat;
      common_arg.results = 0;
      AlgHmcPhi hmc(lat,&common_arg,&hmd_arg);
      
      init_g_upd_cnt = lat.GupdCnt();
      int therm;
      for(therm=0; therm<num_therm_no_met; therm++){
	hmc.run();
	cur_g_upd_cnt = lat.GupdCnt() - init_g_upd_cnt;
      }
      if( (fp = fopen("info.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "info.dat");
      }
      fprintf(fp,"%d %d %f\n",
	      therm, cur_g_upd_cnt, 
	      float(cur_g_upd_cnt) / float(therm) );
      fclose(fp);

      init_g_upd_cnt = lat.GupdCnt();
   }


  //----------------------------------------------------------------
  // Thermalize gauge field using HMC Phi with metropolis
  //----------------------------------------------------------------
  hmd_arg.metropolis = METROPOLIS_YES;
  if(num_therm_met != 0)
    {
      GwilsonFnone lat;
      common_arg.results = 0;
      AlgHmcPhi hmc(lat,&common_arg,&hmd_arg);
      
      init_g_upd_cnt = lat.GupdCnt();
      int therm;
      for(therm=0; therm<num_therm_met; therm++){
	hmc.run();
	cur_g_upd_cnt = lat.GupdCnt() - init_g_upd_cnt;
      }
      
      if( (fp = fopen("info.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "info.dat");
      }
      fprintf(fp,"%d %d %f\n",
	      therm, cur_g_upd_cnt, 
	      float(cur_g_upd_cnt) / float(therm) );
      fclose(fp);

      init_g_upd_cnt = lat.GupdCnt();
    }


  //----------------------------------------------------------------
  // Loop over the number of measurements
  //----------------------------------------------------------------
  for(int meas = 0; meas< num_meas; meas++){

    //--------------------------------------------------------------
    // Evolve gauge field using HMC Phi
    //--------------------------------------------------------------
    hmd_arg.metropolis = METROPOLIS_YES;
    if(num_sweeps != 0)
      {
	GwilsonFnone lat;
	common_arg.results = 0;
	AlgHmcPhi hmc(lat,&common_arg,&hmd_arg);
	
	for(int sweeps=0; sweeps<num_sweeps; sweeps++){
	  hmc.run();
	  iter = iter + 1;
	}
	
	cur_g_upd_cnt = lat.GupdCnt() - init_g_upd_cnt;

	if( (fp = fopen("info.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "info.dat");
	}
	fprintf(fp,"%d %d %f\n",
		iter, cur_g_upd_cnt, 
		float(cur_g_upd_cnt) / float(iter) );	
	fclose(fp);
      }


    //----------------------------------------------------------------
    // Loop over mass and ls 
    //----------------------------------------------------------------
    Float mass = mass_start;
    for(int m=0; m<num_mass; m++){
      int ls = min_ls;
      for(int ils=0; ils<num_ls; ils++){
	GJP.SnodeSites(ls);


	GwilsonFdwf lat;
	{
	  //------------------------------------------------------------
	  // Run AlgPbp 
	  //------------------------------------------------------------
	  if(pbp_num_hits != 0)
	    {
	      common_arg.results = CAST_AWAY_CONST("pbp.dat");
	      pbp_arg.mass[0] = mass;
	      pbp_arg.src_u_s = 0;
	      pbp_arg.src_l_s = ls-1;
	      pbp_arg.snk_u_s = ls-1;
	      pbp_arg.snk_l_s = 0;
	      AlgPbp pbp(lat,&common_arg,&pbp_arg);
	      for(int i=0; i<pbp_num_hits; i++){
		pbp.run();
	      }
	    }

	  //------------------------------------------------------------
	  // Run spectroscopy
	  //------------------------------------------------------------
	  if(spect_yes == 1){
	    {
	      common_arg.results = &wout;
	      w_spect_arg.cg.mass = mass;
	      AlgWspect ws(lat,&common_arg,&w_spect_arg);
	      
	      ws.run();
	    }
	  }


	}

	ls = ls + ls_step;
      }
      mass = mass + mass_step;
    }



  }
  return(0);
}





