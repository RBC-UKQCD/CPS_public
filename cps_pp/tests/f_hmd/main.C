#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-04-30 12:18:01 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_hmd/main.C,v 1.5 2004-04-30 12:18:01 zs Exp $
//  $Id: main.C,v 1.5 2004-04-30 12:18:01 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_hmd/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#include <stdio.h>
#include<util/lattice.h>
#include<alg/alg_hmd.h>
#include<alg/do_arg.h>
#include<alg/ghb_arg.h>

CPS_START_NAMESPACE
GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;
CPS_END_NAMESPACE

USING_NAMESPACE_CPS 


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

  do_arg.beta = 5.5;
  do_arg.dwf_height = 0.9;
  do_arg.clover_coeff = 2.0171;


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

  VRB.ActivateLevel(VERBOSE_RESULT_LEVEL);
  VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);
  VRB.ActivateLevel(VERBOSE_CLOCK_LEVEL);
  VRB.ActivateLevel(VERBOSE_RNGSEED_LEVEL);


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
	common_arg.set_filename("hmc.dat");
	
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
      common_arg.set_filename("hmd.dat");
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
      common_arg.set_filename("hmc.dat");
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
      common_arg.set_filename("hmc.dat");
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





  


