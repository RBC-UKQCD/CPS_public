#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-04-30 12:18:01 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/threept/main.C,v 1.3 2004-04-30 12:18:01 zs Exp $
//  $Id: main.C,v 1.3 2004-04-30 12:18:01 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/threept/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  test for alg
 */

CPS_END_NAMESPACE
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc.h>
CPS_START_NAMESPACE
#endif
CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/random.h>
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

  do_arg.beta = 6.0;
  do_arg.dwf_height = 1.8;


  //----------------------------------------------------------------
  // Initialize the GJP class
  //----------------------------------------------------------------
  GJP.Initialize(do_arg);


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
