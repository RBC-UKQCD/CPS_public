#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-03-09 19:18:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_p4_pbp/main.C,v 1.2 2005-03-09 19:18:47 chulwoo Exp $
//  $Id: main.C,v 1.2 2005-03-09 19:18:47 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_p4_pbp/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#include<util/lattice.h>
#include<alg/alg_pbp.h>
#include<alg/do_arg.h>
#include<util/qcdio.h>

const int MAX_FILENAME = 200;
char fprefix[MAX_FILENAME] = "";


USING_NAMESPACE_CPS

int main(int argc,char *argv[])
{

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;
  int nx,ny,nz,nt;

    if (argc < 5){
        ERR.General("f_stag_test","main()","usage: %s nx ny nz nt\n",argv[0]);
    }
    sscanf(argv[1],"%d",&nx);
    sscanf(argv[2],"%d",&ny);
    sscanf(argv[3],"%d",&nz);
    sscanf(argv[4],"%d",&nt);


#ifdef PARALLEL
  do_arg.x_node_sites = nx/SizeX();
  do_arg.y_node_sites = ny/SizeY();
  do_arg.z_node_sites = nz/SizeZ();
  do_arg.t_node_sites = nz/SizeT();
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
#else
  do_arg.x_node_sites = nx;
  do_arg.y_node_sites = ny;
  do_arg.z_node_sites = nz;
  do_arg.t_node_sites = nt;
  do_arg.s_node_sites = 0;
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
  do_arg.beta = 6.0;

  do_arg.p4_KS = 3.0/8.0;
  do_arg.p4_knight = 1.0/48.0;
  do_arg.p4_3staple = 1e-12;
  do_arg.p4_5staple = 1e-12;
  do_arg.p4_7staple = 1e-12;
  do_arg.p4_lepage = 1e-12;


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
  PbpArg pbp_arg;

  pbp_arg.pattern_kind = ARRAY;
  pbp_arg.n_masses = 2;
  pbp_arg.mass[0] = 0.1;
  pbp_arg.mass[1] = .01;
  pbp_arg.stop_rsd = 1.0E-12;
  pbp_arg.max_num_iter = 500;

  //----------------------------------------------------------------
  // Run AlgPbp
  //----------------------------------------------------------------
  {
    GnoneFp4 lat;
    common_arg.set_filename("pbp.dat");
    AlgPbp pbp(lat,&common_arg,&pbp_arg);
    int num_hits;
    int i;

    num_hits = 1;
    for(i=0; i<num_hits; i++){
      pbp.run();
    }

  }

  return 1;
}







