#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:10 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_hmd/main.C,v 1.10 2004-08-18 11:58:10 zs Exp $
//  $Id: main.C,v 1.10 2004-08-18 11:58:10 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_hmd/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#include <util/qcdio.h>
#include <stdlib.h>
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

static int nx,ny,nz,nt,ns;
int main(int argc,char *argv[])
{
#if TARGET==QCDOC
DefaultSetup();
#endif
  FILE *fp;
  if (argc<6) {printf("usage: %s nx ny nz nt ns\n",argv[0]);exit(-2);}
  sscanf(argv[1],"%d",&nx);
  sscanf(argv[2],"%d",&ny);
  sscanf(argv[3],"%d",&nz);
  sscanf(argv[4],"%d",&nt);
  sscanf(argv[5],"%d",&ns);
  printf("sizes = %d %d %d %d %d\n",nx,ny,nz,nt,ns);

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

#ifdef PARALLEL
  do_arg.x_node_sites = nx/SizeX();
  do_arg.y_node_sites = ny/SizeY();
  do_arg.z_node_sites = nz/SizeZ();
  do_arg.t_node_sites = nt/SizeT();
  do_arg.s_node_sites = ns/SizeS();
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = SizeS();
#else
  do_arg.x_node_sites = nx;
  do_arg.y_node_sites = ny;
  do_arg.z_node_sites = nz;
  do_arg.t_node_sites = nt;
  do_arg.s_node_sites = ns;
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

  do_arg.beta = 5.5;
  do_arg.dwf_height = 0.9;
  do_arg.clover_coeff = 2.0171;

#if TARGET==cpsMPI
    MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes);    
    using MPISCU::fprintf;
    using MPISCU::printf;
#endif
  
  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------

  VRB.ActivateLevel(VERBOSE_RESULT_LEVEL);
  VRB.ActivateLevel(VERBOSE_FUNC_LEVEL);
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


#if 1
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

      if( (fp = Fopen("plaq.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "plaq.dat");
      }
      Fprintf(fp, "%e\n", float(aver_plaq0));
      Fclose(fp);

      for (int i = 0; i < 1; ++i) {
	hmc_phi.run();
	sum_plaq0 = lat.SumReTrPlaq();
	aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);
	if( (fp = Fopen("plaq.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "plaq.dat");
	}
	Fprintf(fp,"%e\n", float(aver_plaq0));
	Fclose(fp);
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

      if( (fp = Fopen("plaq.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "plaq.dat");
      }
      Fprintf(fp, "%e\n", float(aver_plaq0));
      Fclose(fp);

      for (int i = 0; i < 1; ++i) {
	hmd_r.run();
	sum_plaq0 = lat.SumReTrPlaq();
	aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);
	if( (fp = Fopen("plaq.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "plaq.dat");
	}
	Fprintf(fp,"%e\n", float(aver_plaq0));
	Fclose(fp);
      }
    }
  }
#endif


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

      if( (fp = Fopen("plaq.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "plaq.dat");
      }
      Fprintf(fp, "%e\n", float(aver_plaq0));
      Fclose(fp);

      for (int i = 0; i < 1; ++i) {
	hmc_phi.run();
	sum_plaq0 = lat.SumReTrPlaq();
	aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);
	if( (fp = Fopen("plaq.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "plaq.dat");
	}
	Fprintf(fp,"%e\n", float(aver_plaq0));
	Fclose(fp);
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

      if( (fp = Fopen("plaq.dat", "a")) == NULL ) {
	ERR.FileA(" ","main", "plaq.dat");
      }
      Fprintf(fp, "%e\n", float(aver_plaq0));
      Fclose(fp);

      for (int i = 0; i < 1; ++i) {
	hmc_phi.run();
	sum_plaq0 = lat.SumReTrPlaq();
	aver_plaq0 = 1.-sum_plaq0/(18.0*total_sites);
	if( (fp = Fopen("plaq.dat", "a")) == NULL ) {
	  ERR.FileA(" ","main", "plaq.dat");
	}
	Fprintf(fp,"%e\n", float(aver_plaq0));
	Fclose(fp);
      }
    }
  }

  
  return(1);  
}





  


