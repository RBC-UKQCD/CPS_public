//In this code we demonstrate the relationship 
// \prop^(1-j,1-k) (y,z) = (-1)^|j-k| \gamma^5 C [ \prop^(j,k) (y,z) ]* C^\dagger \gamma^5

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include <util/lattice.h>
#include <util/time_cps.h>
#include <alg/do_arg.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>
#include <alg/alg_plaq.h>
#include <alg/alg_rnd_gauge.h>
#include <alg/threept_arg.h>
#include <alg/threept_prop_arg.h>
#include <alg/alg_threept.h>
#include <util/smalloc.h>
#if(0==1)
 #include <ReadLattice.h>
 #include <WriteLattice.h>
#endif

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <util/command_line.h>
#include <sstream>

#include<unistd.h>
#include<config.h>

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include <alg/alg_fix_gauge.h>
#include <alg/fix_gauge_arg.h>

#include <util/data_shift.h>

#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>
#include <alg/prop_dft.h>

using namespace std;
USING_NAMESPACE_CPS

void print(const WilsonMatrix &w){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      Complex c = w(i,0,j,0);
      printf("(%.4f %.4f) ",c.real(),c.imag());
    }
    printf("\n");
  }
  printf("\n");
}

bool test_equals(const WilsonMatrix &a, const WilsonMatrix &b, const double &eps){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int aa=0;aa<3;aa++){
	for(int bb=0;bb<3;bb++){
	  Complex ca = a(i,aa,j,bb);
	  Complex cb = b(i,aa,j,bb);
	  if( fabs(ca.real()-cb.real()) > eps || fabs(ca.imag()-cb.imag()) > eps ) return false;
	}
      }
    }
  }
  return true;
}



void global_coord(const int &site, int *into_vec){
  int shift_x = GJP.XnodeCoor()*GJP.XnodeSites();
  int shift_y = GJP.YnodeCoor()*GJP.YnodeSites();
  int shift_z = GJP.ZnodeCoor()*GJP.ZnodeSites();
  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  //Local lattice dimensions:
  int size_x = GJP.XnodeSites();
  int size_y = GJP.YnodeSites();
  int size_z = GJP.ZnodeSites();
  int size_t = GJP.TnodeSites();
  int size_xy = size_x*size_y;
  int spatial_vol = (GJP.VolNodeSites()/GJP.TnodeSites()); // =size_x*size_y_size_z

  into_vec[3] = site/spatial_vol + shift_t;
  into_vec[2] = (site%spatial_vol)/size_xy + shift_z;
  into_vec[1] = (site%size_xy)/size_x + shift_y;
  into_vec[0] = site%size_x + shift_x;
}

int get_site(int pos[4]){
  int locpos[4];
  for(int i=0;i<4;i++){
    int shift = GJP.NodeCoor(i)*GJP.NodeSites(i);
    locpos[i] = pos[i] - shift;
    if(locpos[i] < 0 || locpos[i] >= GJP.NodeSites(i)) return -1; //not on node
  }
  return locpos[0] + GJP.XnodeSites()*(locpos[1] + GJP.YnodeSites()*(locpos[2] + GJP.ZnodeSites()*locpos[3]));
}
int get_site(const int &x, const int &y, const int &z, const int &t){
  int pos[4] = {x,y,z,t};
  return get_site(pos);
}


const Float Pi_const(3.141592654);

Rcomplex sink_phasefac(int *momphase, int *pos){
  //momphase is the sum of the phase factors from the propagators forming the contraction
  //NOTE: In G-parity directions, momentum is discretised in odd units of \pi/2L rather than even/odd units of 2\pi/L (periodic/antiperiodic).
  Float pdotx = 0.0;
  
  for(int d=0;d<3;d++){
    Float mom_unit;
    if(GJP.Bc(d) == BND_CND_GPARITY) mom_unit = Pi_const/( (Float) 2*GJP.Nodes(d)*GJP.NodeSites(d));
    else if(GJP.Bc(d) == BND_CND_PRD) mom_unit = 2.0*Pi_const/( (Float) GJP.Nodes(d)*GJP.NodeSites(d));
    else if(GJP.Bc(d) == BND_CND_APRD) mom_unit = Pi_const/( (Float) GJP.Nodes(d)*GJP.NodeSites(d));
    else ERR.General("","sink_phasefac(int *,const int &)","Unknown boundary condition\n");
    
    pdotx += momphase[d]*pos[d]*mom_unit;
  }
  return Rcomplex(cos(pdotx),sin(pdotx));
}
Rcomplex sink_phasefac(int *momphase,const int &site){
  int pos[4]; global_coord(site,pos);
  return sink_phasefac(momphase,pos);
}
void sum_momphase(int *into, PropagatorContainer &prop, const bool &is_cconj){
  int propmom[3]; prop.momentum(propmom);
  if(is_cconj){
    for(int i=0;i<3;i++) into[i]-=propmom[i];
  }else{
    for(int i=0;i<3;i++) into[i]+=propmom[i];
  }
}

void test_comb(const JobPropagatorArgs &prop_args, const int &Pidx, const int &Aidx, const int &Combidx, const PropCombination &comb, Lattice &latt){
  PropagatorContainer &q_comb_pc = PropManager::getProp(prop_args.props.props_val[Combidx].generics.tag);
  QPropW &q_comb_qpw = q_comb_pc.getProp(latt);

  PropagatorContainer &q_P_pc = PropManager::getProp(prop_args.props.props_val[Pidx].generics.tag);
  QPropW &q_P_qpw = q_P_pc.getProp(latt);

  PropagatorContainer &q_A_pc = PropManager::getProp(prop_args.props.props_val[Aidx].generics.tag);
  QPropW &q_A_qpw = q_A_pc.getProp(latt);

  QPropW cmb(q_P_qpw);
  if(comb == A_PLUS_B) cmb.Average(q_A_qpw);
  else cmb.LinComb(q_A_qpw,0.5,-0.5);

  bool fail(false);
  for(int i=0;i<GJP.VolNodeSites();i++){
    for(int f=0;f<2;f++){
      WilsonMatrix P(q_P_qpw.SiteMatrix(i,f));
      WilsonMatrix A(q_A_qpw.SiteMatrix(i,f));
      
      WilsonMatrix cpa;
      if(comb == A_PLUS_B) cpa = 0.5*(P+A);
      else cpa =0.5*(P-A);

      WilsonMatrix C(q_comb_qpw.SiteMatrix(i,f));

      WilsonMatrix cmbval(cmb.SiteMatrix(i,f));

      if(!test_equals(cpa,C,1e-08)){
	printf("Combination fail at x=%d, f=%d: \n",i,f);
	printf("Combining WilsonMatrices got:\n");
	print(cpa);
	printf("Combined prop got:\n");
	print(C);
	printf("Combining QPropW got:\n");
	print(cmbval);
	
	fail=true;
      }
    }
  }
  if(fail){ exit(-1); }
}


void test_props(const JobPropagatorArgs &prop_args, const int &f0idx, const int &f1idx, Lattice &latt){
  PropagatorContainer &q_f0_pc = PropManager::getProp(prop_args.props.props_val[f0idx].generics.tag);
  PropagatorContainer &q_f1_pc = PropManager::getProp(prop_args.props.props_val[f1idx].generics.tag);
  
  QPropW &q_f0_qpw = q_f0_pc.getProp(latt);
  QPropW &q_f1_qpw = q_f1_pc.getProp(latt);
  
  int shift_x = GJP.XnodeCoor()*GJP.XnodeSites();
  int shift_y = GJP.YnodeCoor()*GJP.YnodeSites();
  int shift_z = GJP.ZnodeCoor()*GJP.ZnodeSites();
  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  //Local lattice dimensions:
  int size_x = GJP.XnodeSites();
  int size_y = GJP.YnodeSites();
  int size_z = GJP.ZnodeSites();
  int size_t = GJP.TnodeSites();
  int size_xy = size_x*size_y;
  int spatial_vol = (GJP.VolNodeSites()/GJP.TnodeSites()); // =size_x*size_y_size_z

  for(int i=0;i<GJP.VolNodeSites();i++){
    int pos[4];
    pos[3] = i/spatial_vol + shift_t;
    pos[2] = (i%spatial_vol)/size_xy + shift_z;
    pos[1] = (i%size_xy)/size_x + shift_y;
    pos[0] = i%size_x + shift_x;

    printf("(%d,%d,%d,%d):\n",pos[0],pos[1],pos[2],pos[3]);

    WilsonMatrix g_0_0 = q_f0_qpw.SiteMatrix(i,0);
    WilsonMatrix g_1_0 = q_f0_qpw.SiteMatrix(i,1);
    WilsonMatrix g_0_1 = q_f1_qpw.SiteMatrix(i,0);
    WilsonMatrix g_1_1 = q_f1_qpw.SiteMatrix(i,1);

    //note, due to CPS oddity, A.ccl(1) = C^-1 A,  A.ccl(-1) = C A
    //                whereas, A.ccr(1) = A C,          A.ccr(-1) = A C^-1

    // \prop^(1-j,1-k) (y,z) = (-1)^|j-k| \gamma^5 C [ \prop^(j,k) (y,z) ]* C^\dagger \gamma^5
    {
      printf("Testing 1,1 = f(0,0):");
      WilsonMatrix f_0_0 = g_0_0;
      f_0_0.cconj();
      f_0_0.ccl(-1).gl(-5).ccr(-1).gr(-5);
      bool result = test_equals(g_1_1,f_0_0,1e-08);
      if(!result){
	printf(" false\n");
	print(g_1_1);
	print(f_0_0);
	exit(-1);
      }else printf(" true\n");
    }
    {
      printf("Testing 1,0 = f(0,1):");
      WilsonMatrix f_0_1 = g_0_1;
      f_0_1.cconj();
      f_0_1.ccl(-1).gl(-5).ccr(-1).gr(-5)*=Complex(-1.0,0.0);
      bool result = test_equals(g_1_0,f_0_1,1e-08);
      if(!result){
	printf(" false\n");
	print(g_1_0);
	print(f_0_1);
	exit(-1);
      }else printf(" true\n");
    }
    {
      printf("Testing 0,1 = f(1,0):");
      WilsonMatrix f_1_0 = g_1_0;
      f_1_0.cconj();
      f_1_0.ccl(-1).gl(-5).ccr(-1).gr(-5)*=Complex(-1.0,0.0);
      bool result = test_equals(g_0_1,f_1_0,1e-08);
      if(!result){
	printf(" false\n");
	print(g_0_1);
	print(f_1_0);
	exit(-1);
      }else printf(" true\n");
    }
    {
      printf("Testing 0,0 = f(1,1):");
      WilsonMatrix f_1_1 = g_1_1;
      f_1_1.cconj();
      f_1_1.ccl(-1).gl(-5).ccr(-1).gr(-5);
      bool result = test_equals(g_0_0,f_1_1,1e-08);
      if(!result){
	printf(" false\n");
	print(g_0_0);
	print(f_1_1);
	exit(-1);
      }else printf(" true\n");
    }


    
    {
      printf("Testing new 1,1 = f'(0,0):");
      WilsonMatrix f_0_0 = g_0_0;
      f_0_0.cconj();
      f_0_0.ccl(-1).gl(-5).ccr(1).gr(-5);
      bool result = test_equals(g_1_1,f_0_0,1e-08);
      if(!result){
	printf(" false\n");
	print(g_1_1);
	print(f_0_0);
	exit(-1);
      }else printf(" true\n");
    }
    {
      printf("Testing new 1,0 = f(0,1):");
      WilsonMatrix f_0_1 = g_0_1;
      f_0_1.cconj();
      f_0_1.ccl(-1).gl(-5).ccr(1).gr(-5);
      bool result = test_equals(g_1_0,f_0_1,1e-08);
      if(!result){
	printf(" false\n");
	print(g_1_0);
	print(f_0_1);
	exit(-1);
      }else printf(" true\n");
    }
    {
      printf("Testing new 0,1 = f(1,0):");
      WilsonMatrix f_1_0 = g_1_0;
      f_1_0.cconj();
      f_1_0.ccl(-1).gl(-5).ccr(1).gr(-5);
      bool result = test_equals(g_0_1,f_1_0,1e-08);
      if(!result){
	printf(" false\n");
	print(g_0_1);
	print(f_1_0);
	exit(-1);
      }else printf(" true\n");
    }
    {
      printf("Testing new 0,0 = f(1,1):");
      WilsonMatrix f_1_1 = g_1_1;
      f_1_1.cconj();
      f_1_1.ccl(-1).gl(-5).ccr(1).gr(-5);
      bool result = test_equals(g_0_0,f_1_1,1e-08);
      if(!result){
	printf(" false\n");
	print(g_0_0);
	print(f_1_1);
	exit(-1);
      }else printf(" true\n");
    }

  }
    
}



int main(int argc,char *argv[])
{
  Start(&argc,&argv);
  CommandLine::is(argc,argv);

  bool gparity_X(false);
  bool gparity_Y(false);

  int arg0 = CommandLine::arg_as_int(0);
  printf("Arg0 is %d\n",arg0);
  if(arg0==0){
    gparity_X=true;
    printf("Doing G-parity HMC test in X direction\n");
  }else if(arg0==1){
    printf("Doing G-parity HMC test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
  }

  bool dbl_latt_storemode(false);
  bool save_config(false);
  bool load_config(false);
  bool load_lrg(false);
  bool save_lrg(false);
  char *load_config_file;
  char *save_config_file;
  char *save_lrg_file;
  char *load_lrg_file;
  bool gauge_fix(false);
  bool verbose(false);

  int size[] = {2,2,2,2,2};

  int i=2;
  while(i<argc){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-save_config",15) == 0){
      if(i==argc-1){ printf("-save_config requires an argument\n"); exit(-1); }
      save_config=true;
      save_config_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-load_config",15) == 0){
      if(i==argc-1){ printf("-save_config requires an argument\n"); exit(-1); }
      load_config=true;
      load_config_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-latt",10) == 0){
      if(i>argc-6){
	printf("Did not specify enough arguments for 'latt' (require 5 dimensions)\n"); exit(-1);
      }
      size[0] = CommandLine::arg_as_int(i); //CommandLine ignores zeroth input arg (i.e. executable name)
      size[1] = CommandLine::arg_as_int(i+1);
      size[2] = CommandLine::arg_as_int(i+2);
      size[3] = CommandLine::arg_as_int(i+3);
      size[4] = CommandLine::arg_as_int(i+4);
      i+=6;
    }else if( strncmp(cmd,"-save_double_latt",20) == 0){
      dbl_latt_storemode = true;
      i++;
    }else if( strncmp(cmd,"-load_lrg",15) == 0){
      if(i==argc-1){ printf("-load_lrg requires an argument\n"); exit(-1); }
      load_lrg=true;
      load_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-save_lrg",15) == 0){
      if(i==argc-1){ printf("-save_lrg requires an argument\n"); exit(-1); }
      save_lrg=true;
      save_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-gauge_fix",15) == 0){
      gauge_fix=true;
      i++;   
    }else if( strncmp(cmd,"-verbose",15) == 0){
      verbose=true;
      i++;
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }
  
  printf("Lattice size is %d %d %d %d\n",size[0],size[1],size[2],size[3],size[4]);

  DoArg do_arg;
  do_arg.x_sites = size[0];
  do_arg.y_sites = size[1];
  do_arg.z_sites = size[2];
  do_arg.t_sites = size[3];
  do_arg.s_sites = size[4];
  do_arg.x_node_sites = 0;
  do_arg.y_node_sites = 0;
  do_arg.z_node_sites = 0;
  do_arg.t_node_sites = 0;
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = 0;
  do_arg.y_nodes = 0;
  do_arg.z_nodes = 0;
  do_arg.t_nodes = 0;
  do_arg.s_nodes = 0;
  do_arg.updates = 0;
  do_arg.measurements = 0;
  do_arg.measurefreq = 0;
  do_arg.cg_reprod_freq = 10;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_conf_load_addr = 0x0;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.start_seed_filename = "../rngs/ckpoint_rng.0";
  do_arg.start_conf_filename = "../configurations/ckpoint_lat.0";
  do_arg.start_conf_alloc_flag = 6;
  do_arg.wfm_alloc_flag = 2;
  do_arg.wfm_send_alloc_flag = 2;
  do_arg.start_seed_value = 83209;
  do_arg.beta =   2.25;
  do_arg.c_1 =   -3.3100000000000002e-01;
  do_arg.u0 =   1.0000000000000000e+00;
  do_arg.dwf_height =   1.8000000000000000e+00;
  do_arg.dwf_a5_inv =   1.0000000000000000e+00;
  do_arg.power_plaq_cutoff =   0.0000000000000000e+00;
  do_arg.power_plaq_exponent = 0;
  do_arg.power_rect_cutoff =   0.0000000000000000e+00;
  do_arg.power_rect_exponent = 0;
  do_arg.verbose_level = -1202;//VERBOSE_DEBUG_LEVEL; //-1202;
  do_arg.checksum_level = 0;
  do_arg.exec_task_list = 0;
  do_arg.xi_bare =   1.0000000000000000e+00;
  do_arg.xi_dir = 3;
  do_arg.xi_v =   1.0000000000000000e+00;
  do_arg.xi_v_xi =   1.0000000000000000e+00;
  do_arg.clover_coeff =   0.0000000000000000e+00;
  do_arg.clover_coeff_xi =   0.0000000000000000e+00;
  do_arg.xi_gfix =   1.0000000000000000e+00;
  do_arg.gfix_chkb = 1;
  do_arg.asqtad_KS =   0.0000000000000000e+00;
  do_arg.asqtad_naik =   0.0000000000000000e+00;
  do_arg.asqtad_3staple =   0.0000000000000000e+00;
  do_arg.asqtad_5staple =   0.0000000000000000e+00;
  do_arg.asqtad_7staple =   0.0000000000000000e+00;
  do_arg.asqtad_lepage =   0.0000000000000000e+00;
  do_arg.p4_KS =   0.0000000000000000e+00;
  do_arg.p4_knight =   0.0000000000000000e+00;
  do_arg.p4_3staple =   0.0000000000000000e+00;
  do_arg.p4_5staple =   0.0000000000000000e+00;
  do_arg.p4_7staple =   0.0000000000000000e+00;
  do_arg.p4_lepage =   0.0000000000000000e+00;

  if(verbose) do_arg.verbose_level = VERBOSE_DEBUG_LEVEL;

  if(gparity_X) do_arg.x_bc = BND_CND_GPARITY;
  if(gparity_Y) do_arg.y_bc = BND_CND_GPARITY;

  GJP.Initialize(do_arg);

  SerialIO::dbl_latt_storemode = dbl_latt_storemode;
  
  LRG.Initialize();

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }


  GwilsonFdwf lattice;
					       
  if(!load_config){
    printf("Creating gauge field\n");
    lattice.SetGfieldDisOrd(); //unit gauge
  }else{
    ReadLatticeParallel readLat;
    if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
    readLat.read(lattice,load_config_file);
    if(UniqueID()==0) printf("Config read.\n");
  }

  if(save_config){
    if(UniqueID()==0) printf("Saving config to %s\n",save_config_file);

    QioArg wt_arg(save_config_file,0.001);
    
    wt_arg.ConcurIONumber=32;
    WriteLatticeParallel wl;
    wl.setHeader("disord_id","disord_label",0);
    wl.write(lattice,wt_arg);
    
    if(!wl.good()) ERR.General("main","()","Failed write lattice %s",save_config_file);

    if(UniqueID()==0) printf("Config written.\n");
  }

  if(gauge_fix){
    lattice.FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    lattice.FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }

#define SETUP_ARRAY(OBJ,ARRAYNAME,TYPE,SIZE)	\
  OBJ . ARRAYNAME . ARRAYNAME##_len = SIZE; \
  OBJ . ARRAYNAME . ARRAYNAME##_val = new TYPE [SIZE]

#define ELEM(OBJ,ARRAYNAME,IDX) OBJ . ARRAYNAME . ARRAYNAME##_val[IDX]

  JobPropagatorArgs prop_args;
  SETUP_ARRAY(prop_args,props,PropagatorArg,6);
  
  char* names[6] = {"prop_f0_0000","prop_f0_1111", "prop_f0_1010",
		    "prop_f1_0000", "prop_f1_1111","prop_f1_1010"};
  int pos[6][4] = { {0,0,0,0}, {1,1,1,1}, {1,0,1,0},
		    {0,0,0,0}, {1,1,1,1}, {1,0,1,0} };

  int flav[6] = {0,0,0,1,1,1};

  for(int i=0;i<6;i++){
    PropagatorArg &parg = prop_args.props.props_val[i];
    
    parg.generics.tag = names[i];
    parg.generics.mass = 0.1;
    parg.generics.bc[0] = GJP.Xbc();
    parg.generics.bc[1] = GJP.Ybc();
    parg.generics.bc[2] = GJP.Zbc();
    parg.generics.bc[3] = GJP.Tbc();

    SETUP_ARRAY(parg,attributes,AttributeContainer,2);
    
    ELEM(parg,attributes,0).type = POINT_SOURCE_ATTR;
    PointSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.point_source_attr;
    for(int j=0;j<4;j++) srcarg.pos[j] = pos[i][j];

    ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
    GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
    gparg.flavor = flav[i];
  }

  if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args.props.props_len);

  PropManager::setup(prop_args);
  PropManager::calcProps(lattice);

  {
    //test  g5 G(y,x)^dag g5 = G(x,y)

    { //x = 0000 y=1111
      WilsonMatrix xval_f0_p_src_f0_y = 0.0;
      WilsonMatrix xval_f0_p_src_f1_y = 0.0;
      WilsonMatrix xval_f1_p_src_f0_y = 0.0;
      WilsonMatrix xval_f1_p_src_f1_y = 0.0;

      WilsonMatrix yval_f0_p_src_f0_x = 0.0;
      WilsonMatrix yval_f0_p_src_f1_x = 0.0;
      WilsonMatrix yval_f1_p_src_f0_x = 0.0;
      WilsonMatrix yval_f1_p_src_f1_x = 0.0;

      int ysite = get_site(1,1,1,1);
      int xsite = get_site(0,0,0,0);
      char* src_f0_x = names[0];
      char* src_f1_x = names[3];
      char* src_f0_y = names[1];
      char* src_f1_y = names[4];

      if(ysite!=-1){ //on node
	PropagatorContainer &p_src_f0_x = PropManager::getProp(src_f0_x);
	yval_f0_p_src_f0_x = p_src_f0_x.getProp(lattice).SiteMatrix(ysite,0);
	yval_f1_p_src_f0_x = p_src_f0_x.getProp(lattice).SiteMatrix(ysite,1);

	PropagatorContainer &p_src_f1_x = PropManager::getProp(src_f1_x);
	yval_f0_p_src_f1_x = p_src_f1_x.getProp(lattice).SiteMatrix(ysite,0);
	yval_f1_p_src_f1_x = p_src_f1_x.getProp(lattice).SiteMatrix(ysite,1);
      }
      _FourierProp_helper<WilsonMatrix>::lattice_sum(yval_f0_p_src_f0_x);
      _FourierProp_helper<WilsonMatrix>::lattice_sum(yval_f0_p_src_f1_x);
      _FourierProp_helper<WilsonMatrix>::lattice_sum(yval_f1_p_src_f0_x);
      _FourierProp_helper<WilsonMatrix>::lattice_sum(yval_f1_p_src_f1_x);

      if(xsite!=-1){ //on node
	PropagatorContainer &p_src_f0_y = PropManager::getProp(src_f0_y);
	xval_f0_p_src_f0_y = p_src_f0_y.getProp(lattice).SiteMatrix(xsite,0);
	xval_f1_p_src_f0_y = p_src_f0_y.getProp(lattice).SiteMatrix(xsite,1);

	PropagatorContainer &p_src_f1_y = PropManager::getProp(src_f1_y);
	xval_f0_p_src_f1_y = p_src_f1_y.getProp(lattice).SiteMatrix(xsite,0);
	xval_f1_p_src_f1_y = p_src_f1_y.getProp(lattice).SiteMatrix(xsite,1);
      }
      _FourierProp_helper<WilsonMatrix>::lattice_sum(xval_f0_p_src_f0_y);
      _FourierProp_helper<WilsonMatrix>::lattice_sum(xval_f0_p_src_f1_y);
      _FourierProp_helper<WilsonMatrix>::lattice_sum(xval_f1_p_src_f0_y);
      _FourierProp_helper<WilsonMatrix>::lattice_sum(xval_f1_p_src_f1_y);

      if(!UniqueID()){
	//do comparison only on head node now that everything is synced
	bool fail(false);	

	//test  g5 G(y,x)^dag g5 = G(x,y)
	WilsonMatrix g5_Gyx_g5_f0_f0 = yval_f0_p_src_f0_x;
	g5_Gyx_g5_f0_f0.hconj();
	g5_Gyx_g5_f0_f0.gl(-5).gr(-5);
	
	if(!test_equals(g5_Gyx_g5_f0_f0, xval_f0_p_src_f0_y, 1e-06)){
	  printf("f0 f0 test fail\n");
	  print(g5_Gyx_g5_f0_f0);
	  print(xval_f0_p_src_f0_y);
	  fail = true;
	}

	WilsonMatrix g5_Gyx_g5_f0_f1 = yval_f0_p_src_f1_x;
	g5_Gyx_g5_f0_f1.hconj();
	g5_Gyx_g5_f0_f1.gl(-5).gr(-5);
	
	if(!test_equals(g5_Gyx_g5_f0_f1, xval_f1_p_src_f0_y, 1e-06)){
	  printf("f0 f1 test fail\n");
	  print(g5_Gyx_g5_f0_f1);
	  print(xval_f1_p_src_f0_y);
	  fail = true;
	}

	WilsonMatrix g5_Gyx_g5_f1_f0 = yval_f1_p_src_f0_x;
	g5_Gyx_g5_f1_f0.hconj();
	g5_Gyx_g5_f1_f0.gl(-5).gr(-5);
	
	if(!test_equals(g5_Gyx_g5_f1_f0, xval_f0_p_src_f1_y, 1e-06)){
	  printf("f1 f0 test fail\n");
	  print(g5_Gyx_g5_f1_f0);
	  print(xval_f0_p_src_f1_y);
	  fail = true;
	}

	WilsonMatrix g5_Gyx_g5_f1_f0 = yval_f1_p_src_f0_x;
	g5_Gyx_g5_f1_f0.hconj();
	g5_Gyx_g5_f1_f0.gl(-5).gr(-5);
	
	if(!test_equals(g5_Gyx_g5_f1_f0, xval_f0_p_src_f1_y, 1e-06)){
	  printf("f1 f0 test fail\n");
	  print(g5_Gyx_g5_f1_f0);
	  print(xval_f0_p_src_f1_y);
	  fail = true;
	}

	if(fail){
	  printf("g5 herm test failed\n");
	  exit(-1);
	}else printf("g5 herm test passed\n");
	  
      }
	


    }

      
      
  }










  if(gauge_fix) lattice.FixGaugeFree();

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}



