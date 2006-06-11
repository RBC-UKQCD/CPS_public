#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006-06-11 05:35:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_wilson_eig/main.C,v 1.10 2006-06-11 05:35:08 chulwoo Exp $
//  $Id: main.C,v 1.10 2006-06-11 05:35:08 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_wilson_eig/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#include<util/lattice.h>
#include<alg/alg_ghb.h>
#include<alg/alg_plaq.h>
#include<alg/alg_eig.h>
#include<alg/do_arg.h>


USING_NAMESPACE_CPS


int main(int argc,char *argv[])
{
  
  Start();
  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  do_arg.x_node_sites = 2;
  do_arg.y_node_sites = 2;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 2;

#ifdef PARALLEL
  do_arg.x_nodes = 2; //SizeX();
  do_arg.y_nodes = 1; //SizeY();
  do_arg.z_nodes = 1; //SizeZ();
  do_arg.t_nodes = 1; //SizeT();
#else
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
#endif

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 5.3;

#if TARGET==cpsMPI
    MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes);    
    using MPISCU::fprintf;
    using MPISCU::printf;
#endif

  GJP.Initialize(do_arg);

  VRB.DeactivateAll();
  VRB.Level(VERBOSE_CLOCK_LEVEL);
  VRB.Level(VERBOSE_FUNC_LEVEL);


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

      common_arg_ghb.set_filename("ghb.dat");
      AlgGheatBath ghb(lat,&common_arg_ghb,&ghb_arg);
      common_arg_plaq.set_filename("plaq.dat");
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
      common_arg_plaq.set_filename("plaq.dat");
      AlgPlaq plaq(lat,&common_arg_plaq,&plaq_arg);
	plaq.run();

    }
#endif
    {
      common_arg_eig.set_filename("eig.dat");
      AlgEig eig(lat,&common_arg_eig,&eig_arg);
      eig.run();
    }
  }
  End();
  return(0);
}







