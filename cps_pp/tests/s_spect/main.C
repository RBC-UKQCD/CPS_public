#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-09-21 20:16:55 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/s_spect/main.C,v 1.10 2004-09-21 20:16:55 chulwoo Exp $
//  $Id: main.C,v 1.10 2004-09-21 20:16:55 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/s_spect/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#include <util/qcdio.h>
#include<util/lattice.h>
#include<alg/alg_hmd.h>
#include<alg/alg_s_spect.h>
#include<alg/do_arg.h>
#include<alg/alg_fix_gauge.h>
#include<alg/aots_s.h>




USING_NAMESPACE_CPS




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
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 2;
  do_arg.t_nodes = 1;
#else
  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 4;
  do_arg.t_node_sites = 6;
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
#endif
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 5.375;

#if TARGET==cpsMPI
    MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes);
    using MPISCU::fprintf;
    using MPISCU::printf;
#endif

  GJP.Initialize(do_arg);

  
  

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
        if( 0 == (action = Fopen("action.dat", "a")) ) 
	    ERR.FileA("", "main", "action.dat");
		
	Fprintf(action, "%e\n", (Float)aver_plaq0);
        Fclose(action);

        //------------------------------------------------------------
	// Run stag spectrocopy with point source
	//------------------------------------------------------------

        printf("main: spectroscopy staggered with Point Source\n");


        printf("Aots starts....\n");
        for(Aots it(0, 1, 2); it; ++it) {
          {
	    common_arg.set_filename("quark1.dat");
            AlgStagQuark sq(lat, &common_arg, &p_qarg1, it);
            sq.run();
          }
          {
	    common_arg.set_filename("mes11.dat");
            AlgStagMeson meson(lat, &common_arg, &marg1, it);
            meson.run();
          }
          {
	    common_arg.set_filename("nuc111.dat");
            AlgStagNucleon nucleon(lat, &common_arg, &narg1, it);
            nucleon.run();
          }
          {
	    common_arg.set_filename("quark1.dat");
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
	    common_arg.set_filename("quark3.dat");
	    AlgStagQuark sq(lat, &common_arg, &wz_qarg1, it2);
	    sq.run();
	  }

	  {
	    common_arg.set_filename("mes33.dat");
	    AlgStagMeson meson(lat, &common_arg, &marg3, it2);
	    meson.run();
	  }
	  {
	    common_arg.set_filename("nuc333.dat");
	    AlgStagNucleon nucleon(lat, &common_arg, &narg3, it2);
	    nucleon.run();
	  }
	  {
	    common_arg.set_filename("quark3.dat");
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
	    common_arg.set_filename("quark4.dat");
	    AlgStagQuark sq(lat, &common_arg, &wz_qarg2, it3);
	    sq.run();
	  }
	  {
	    common_arg.set_filename("mes44.dat");
	    AlgStagMeson meson(lat, &common_arg, &marg4, it3);
	    meson.run();
	  }
	  {
	    common_arg.set_filename("nuc444.dat");
	    AlgStagNucleon nucleon(lat, &common_arg, &narg4, it3);
	    nucleon.run();
	  }
	  {
	    common_arg.set_filename("nonlocal444.dat");
	    AlgStagNonLocal nonlocal(lat, &common_arg, &darg4, it3);
	    nonlocal.run();
	  }
	  {
	    common_arg.set_filename("quark4.dat");
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

