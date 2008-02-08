#include<config.h>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/random.h>
#include<alg/alg_ghb.h>
#include<alg/alg_plaq.h>
#include<alg/qpropw.h>
#include<alg/qpropw_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_fix_gauge.h>


USING_NAMESPACE_CPS

int main(int argc,char *argv[])
{
 
  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;
  FixGaugeArg fix_arg;
  QPropWArg wq_arg;
  QPropWArg wq_arg_back;

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
  do_arg.s_node_sites = 2;
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

  //  do_arg.beta = 0.8;
  do_arg.beta = 5.8;
  //  do_arg.colors = 3  ;
  //  do_arg.c_1    = -1.4069 ;
  // do_arg.Iwasaki_style = 1 ;
  // do_arg.u0 = 1.0 ;
  do_arg.dwf_height = 1.80;


  //----------------------------------------------------------------
  // Initialize the GJP class
  //----------------------------------------------------------------
  GJP.Initialize(do_arg);

  //----------------------------------------------------------------
  // Set arguments and evolution information to run light quarks
  //----------------------------------------------------------------
  wq_arg.t = 0;
  wq_arg.x = 0;
  wq_arg.y = 0;
  wq_arg.z = 0;
  wq_arg.cg.mass = 0.04;
  wq_arg.cg.max_num_iter = 1200;
  wq_arg.cg.stop_rsd = 1.0e-8;
  wq_arg.GaugeFixSrc = 0;
  wq_arg.GaugeFixSnk = 0;
  wq_arg.StoreMidpointProp = 0;
  wq_arg.SaveProp = 0;
  wq_arg.DoHalfFermion = 0;

  wq_arg_back.t = 3;
  wq_arg_back.x = 0;
  wq_arg_back.y = 0;
  wq_arg_back.z = 0;
  wq_arg_back.cg.mass = wq_arg.cg.mass;
  wq_arg_back.cg.max_num_iter = 1200;
  wq_arg_back.cg.stop_rsd = 1.0e-8;
  wq_arg_back.GaugeFixSrc = 0;
  wq_arg_back.GaugeFixSnk = 0;
  wq_arg_back.StoreMidpointProp = 0;
  wq_arg_back.SaveProp = 0;
  wq_arg_back.DoHalfFermion = 0;
  
  //----------------------------------------------------------------
  // Set fix gauge argument and default values
  //----------------------------------------------------------------
  fix_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
  fix_arg.hyperplane_start = wq_arg.t;
  fix_arg.hyperplane_step = 1;
  fix_arg.hyperplane_num = GJP.Tnodes()*GJP.TnodeSites();
  fix_arg.stop_cond = 1e-3;
  fix_arg.max_iter_num = 10;


  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------
  CommonArg common_arg;
  CommonArg common_arg_ghb;
  CommonArg common_arg_plaq;

  NoArg plaq_arg;

//-------------- Check plaq --------------------
{ 
  GimprRectFnone lattice;
  lattice.EnableLinkBuffer(1000);

  AlgPlaq plaq(lattice,&common_arg,&plaq_arg);
  plaq.run();
}

  
  {

    //    GimprRectFdwf lat;
    GwilsonFdwf lat;
    //   GwilsonFwilson lat;
    WilsonMatrix temp_mtrx;

    common_arg.results = CAST_AWAY_CONST("qpropw.dat");
    FILE *fp;
    if( (fp = Fopen((char *)(common_arg.results), "a")) == NULL ) {
      ERR.FileA("main", "main", (char *)(common_arg.results) );
    }

    Fprintf(fp, "do_arg.dwf_height= %e\n", do_arg.dwf_height);
    Fprintf(fp, "prop stop_rsd= %e\n",wq_arg.cg.stop_rsd);
    Fprintf(fp, "gauge fix type= %d\n",fix_arg.fix_gauge_kind);
    Fprintf(fp, "gauge fix stop_cond= %e\n",fix_arg.stop_cond);
    Fprintf(fp, "t_src= %d\n",wq_arg.t);
    Fprintf(fp, "t_sink= %d\n",wq_arg_back.t);
    Fclose(fp);

    AlgFixGauge fix_gauge(lat,&common_arg,&fix_arg);
    fix_gauge.run();
    

    
    //------------------------------------------------------------
    // Calculate the propagators
    //------------------------------------------------------------
    char *src_type = (char*) "Point";
    QPropWWallSrc prop_wq(lat, &wq_arg, &common_arg);
    QPropWWallSrc prop_wq_back(lat, &wq_arg_back, &common_arg);
 
    
    //------------------------------------------------------------
    // Calculate two-point function
    //------------------------------------------------------------
    int t, i;
    int time_size = GJP.Tnodes()*GJP.TnodeSites();
    char *cname = "TwoPointFunc";
    char *fname = "calcTwoPointFunc()";


    Rcomplex* func_tp=(Rcomplex*)smalloc(time_size*sizeof(Rcomplex));
    if(func_tp == 0)
    ERR.Pointer(cname, fname, "func_tp");
    VRB.Smalloc(cname, fname,
              "func_tp",func_tp, time_size*sizeof(Rcomplex));

    for(t=0; t<time_size; t++) func_tp[t] = 0.0;

    for(i=0; i< GJP.VolNodeSites(); i++)
    {
      t = i/(GJP.VolNodeSites()/GJP.TnodeSites()); // t coordinate in node 
      t += GJP.TnodeCoor()*GJP.TnodeSites(); // plus t node offset physical t

      temp_mtrx = prop_wq[i];
      temp_mtrx.hconj() ;
   
      func_tp[t] += Trace(temp_mtrx, prop_wq_back[i]) ;

    } // volume sum

    
    for(t=0; t<time_size; t++){
      slice_sum((Float*)&func_tp[t], 2, 99);
    }

    //------------------ Write to file ---------------------------
    char cf_name[1024];
    sprintf(cf_name, "twopt.dat");
    fp = fopen(cf_name, "a");
    fprintf(fp,"Two-point function \n MASSES:  %e   %e SOURCE: %s\n",
	     (float)wq_arg.cg.mass ,(float)wq_arg_back.cg.mass, src_type);
    for(t=0; t<time_size; t++){
      fprintf(fp,"%d %e %e\n", t, func_tp[t].real(),func_tp[t].imag() );	         
    }
    fclose(fp);
    
    //------------- Free memory ----------------------------------
    printf("AlgFixGauge.free() starts... \n");
    fix_gauge.free();
  }
  
  return(0);
}
