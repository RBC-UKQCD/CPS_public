#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:31 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/g_hb/main.C,v 1.10 2013-04-05 17:46:31 chulwoo Exp $
//  $Id: main.C,v 1.10 2013-04-05 17:46:31 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/g_hb/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#include<util/lattice.h>
#include<alg/do_arg.h>
#include<alg/alg_ghb.h>
#include<alg/no_arg.h>
#include<alg/alg_plaq.h>
#include<comms/sysfunc_cps.h>


USING_NAMESPACE_CPS


int main(int argc,char *argv[])
{
  int ITERATIONS = 10;

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  Start(&argc,&argv);
  DoArg do_arg;

  int nx,ny,nz,nt;
  if (argc<5) {printf("usage: %s nx ny nz nt\n",argv[0]);exit(-2);}
  sscanf(argv[1],"%d",&nx);
  sscanf(argv[2],"%d",&ny);
  sscanf(argv[3],"%d",&nz);
  sscanf(argv[4],"%d",&nt);

  printf("sizes = %d %d %d %d\n",nx,ny,nz,nt);
  do_arg.x_node_sites = nx/SizeX();
  do_arg.y_node_sites = ny/SizeY();
  do_arg.z_node_sites = nz/SizeZ();
  do_arg.t_node_sites = nt/SizeT();
  do_arg.s_node_sites = 2;

  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = 1;

  do_arg.t_bc = BND_CND_PRD;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_APRD;

  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 5.71;
// do_arg.beta = 1000;
  do_arg.dwf_height = 0.9;


  GJP.Initialize(do_arg);


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
      common_arg_ghb.set_filename("ghb.dat");
      AlgGheatBath ghb(lat,&common_arg_ghb,&ghb_arg);
      
      common_arg_plaq.set_filename("plaq.dat");
      AlgPlaq plaq(lat,&common_arg_plaq,&plaq_arg);

      for (int i = 0; i < ITERATIONS; ++i) {
	ghb.run();
	plaq.run();
      }
    }
  }

  return(0);
}






