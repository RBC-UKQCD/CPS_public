//CK: In this test we check that the DWF fermion force term is evaluated correctly

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

#include <util/gparity_singletodouble.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

void setup_double_latt(Lattice &double_latt, Matrix* orig_gfield, bool gparity_X, bool gparity_Y);
void setup_double_rng(bool gparity_X, bool gparity_Y);
void setup_double_matrixfield(Matrix* double_mat, Matrix* orig_mat, int nmat_per_site, bool gparity_X, bool gparity_Y);

void EvolveMomFforceUstarGparity(Lattice &lat, Matrix *mom, Vector *chi, Float mass, Float dt);
void GaugeTransformU(Matrix *gtrans, Lattice &lat);
void GaugeTransformP(Matrix *gtrans, Matrix* mom, Lattice &lat);
void GaugeTransformF(Matrix *gtrans,Vector* ferm,Lattice &lat);

IFloat vect_norm(Vector *v, int n){
  IFloat out(0.0);

  for (int i=0; i< n; i++)
   out += v[i].NormSqGlbSum(6);

  return out;
}

int main(int argc,char *argv[])
{
  Start(&argc,&argv); //initialises QMP

#ifdef HAVE_BFM
  Chroma::initialize(&argc,&argv);
#endif

  CommandLine::is(argc,argv);

  bool gparity_X(false);
  bool gparity_Y(false);

  int arg0 = CommandLine::arg_as_int(0);
  printf("Arg0 is %d\n",arg0);
  if(arg0==0){
    gparity_X=true;
    printf("Doing G-parity HMC test in X direction\n");
  }else{
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

  {
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
  do_arg.verbose_level = -1202;// -1202;//VERBOSE_DEBUG_LEVEL; //-1202;
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

  if(verbose) do_arg.verbose_level = VERBOSE_CLOCK_LEVEL;

  if(gparity_X) do_arg.x_bc = BND_CND_GPARITY;
  if(gparity_Y) do_arg.y_bc = BND_CND_GPARITY;

  GJP.Initialize(do_arg);

  SerialIO::dbl_latt_storemode = dbl_latt_storemode;
  
  LRG.Initialize(); //usually initialised when lattice generated, but I pre-init here so I can load the state from file

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }
  
  GwilsonFdwf* lattice = new GwilsonFdwf;
					       
  if(!load_config){
    printf("Creating gauge field\n");
    lattice->SetGfieldDisOrd(); //unit gauge
  }else{
    ReadLatticeParallel readLat;
    if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
    readLat.read(*lattice,load_config_file);
    if(UniqueID()==0) printf("Config read.\n");
  }

  if(save_config){
    if(UniqueID()==0) printf("Saving config to %s\n",save_config_file);

    QioArg wt_arg(save_config_file,0.001);
    
    wt_arg.ConcurIONumber=32;
    WriteLatticeParallel wl;
    wl.setHeader("disord_id","disord_label",0);
    wl.write(*lattice,wt_arg);
    
    if(!wl.good()) ERR.General("main","()","Failed write lattice %s",save_config_file);

    if(UniqueID()==0) printf("Config written.\n");
  }

  if(gauge_fix){
    lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    lattice->FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }

  #define U_USTAR_EVOLVE_CHECK
#ifdef U_USTAR_EVOLVE_CHECK
  {
    printf("Starting U U* evolve check\n");
    Lattice &lat = *lattice;
    //test the gauge field evolution by calculating the conjugate momentum for the U links
    //and separately for the U* links and compare them
    
    int f_sites = GJP.SnodeSites()*GJP.VolNodeSites() / (lat.FchkbEvl()+1); //half-checkerboard fermion

    int f_vec_count = 2* f_sites * lat.SpinComponents(); //number of 3-vectors in fermion
    size_t f_size = f_vec_count * lat.Colors() * 2;

    int Ncb;
    if(lat.FchkbEvl() == 1) Ncb = 1; 
    else if(lat.FchkbEvl() == 0) Ncb = 2;

    Vector* tmp1 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1","","");
    Vector* tmp2 = (Vector*)smalloc(f_size*sizeof(Float),"tmp2","","");
    Vector* phi = (Vector *) smalloc(f_size*sizeof(Float),"phi","","");
    //~~make a random two-flavour pseudofermion field from the heatbath as usual~~
    LatRanGen LRGbak(LRG);
    lat.RandGaussVector(tmp1, 0.5, Ncb);
    lat.RandGaussVector(tmp2, 0.5, Ncb);
    LRG = LRGbak;

    Float h_init = 0.0;
    // phi <- M_f^\dag (RGV)
    h_init += lat.SetPhi(phi, tmp1, tmp2, 0.2, DAG_YES); //mass is 0.2

    //~~mom field evolve~~

    //invert D^\dag D on phi
    Vector *cg_sol = (Vector*)smalloc(f_size*sizeof(Float),"cg_sol","","");
    CgArg frm_cg_arg_md;
    frm_cg_arg_md.mass = 0.2;
    frm_cg_arg_md.max_num_iter = 5000;
    frm_cg_arg_md.stop_rsd = 1e-06;

    lat.FmatEvlInv(cg_sol, phi, &frm_cg_arg_md, CNV_FRM_NO);

    int g_size = 2*GJP.VolNodeSites() * lat.GsiteSize();

    // //take a copy of the gauge field
     Matrix* gcopy = (Matrix*)smalloc(g_size*sizeof(Float),"gcopy","","");
     memcpy((void*)gcopy,(void*)lat.GaugeField(),g_size*sizeof(Float));

    //evolve U conj mom then copy-conjugate to U* field conj mom
    Matrix* mom = (Matrix*)smalloc(g_size*sizeof(Float),"mom","","");
    for(int i=0;i<2*4*GJP.VolNodeSites();i++) mom[i].ZeroMatrix();

    lat.EvolveMomFforce(mom, cg_sol, 0.5, 1.0); //evolve 

    //restore the gauge field
     memcpy((void*)lat.GaugeField(),(void*)gcopy,g_size*sizeof(Float));
    
    Matrix* mom2 = (Matrix*)smalloc(g_size*sizeof(Float),"mom2","","");
    for(int i=0;i<2*4*GJP.VolNodeSites();i++) mom2[i].ZeroMatrix();
    
    //evolve U* conj mom then copy-conjugate to U field conj mom
    EvolveMomFforceUstarGparity(lat, mom2, cg_sol, 0.5, 1.0);

    //restore the gauge field
     memcpy((void*)lat.GaugeField(),(void*)gcopy,g_size*sizeof(Float));

    int flav_off = GJP.VolNodeSites()*4;

    printf("Checking momentum fields\n");

    //mom and mom2 should be the same
    Float err(0.0);
    
    for(int f=0;f<2;f++){
      for(int t=0;t<GJP.TnodeSites();t++){
	for(int z=0;z<GJP.ZnodeSites();z++){
	  for(int y=0;y<GJP.YnodeSites();y++){
	    for(int x=0;x<GJP.XnodeSites();x++){
	      int pos[4] = {x,y,z,t};
	      for(int mu=0;mu<4;mu++){
		int off = lat.GsiteOffset(pos) + mu + f*flav_off;
		Float* m = (Float*)(mom+off);
		Float* mc = (Float*)(mom2+off);
		if(fabs(*m - *mc) > 1e-08 || fabs(*(m+1) - *(mc+1)) > 1e-08 ){
		  printf("Node %d %d %d %d , Error: f%d (%d %d %d %d), %d: (%f %f), (%f %f)\n",GJP.XnodeCoor(),GJP.YnodeCoor(),GJP.ZnodeCoor(),GJP.TnodeCoor(),f,x,y,z,t,mu,*m,*(m+1),*mc, *(mc+1));
		  err=1.0;
		}
	      }
	    }
	  }
	}
      }
    }
    glb_sum(&err);

    if(err>0.0){
      printf("U U* evolve check failed\n");
      exit(-1);
    }

    sfree(tmp1);
    sfree(tmp2);
    sfree(phi);
    sfree(cg_sol);
    sfree(gcopy);
    sfree(mom);
    sfree(mom2);
    printf("U U* evolve check passed\n");
  }
#endif

  #define GAUGE_TRANSFORM_CHECK
#ifdef GAUGE_TRANSFORM_CHECK
  {
    printf("Starting gauge transform check\n");
    //test the code by performing a gauge transformation and checking the conjugate momentum transforms in the correct way
    //U -> e^{i \pi^a T^a}U
    //V_x U^mu_x V_{x+\mu}^\dagger -> V_x e^{i \pi^a T^a} U^\mu_x  V_{x+\mu}^\dagger = V_x(1 + [i\pi^aT^a] + [i\pi^aT^a][i\pi^bT^b]/2 +....)U^\mu_x V^\dagger_{x+\mu}
    //                                                = V_x(1 + [i\pi^aT^a] + [i\pi^aT^a]V^\dagger V[i\pi^bT^b]/2 +....)V_x^\dagger V_x U^\mu_x V^\dagger_{x+\mu}
    //                                                = e^{iV_x\pi^a T^a V_x^\dagger} V_x U V_{x+\mu}^\dagger
    // so (\pi^a T^a) -> V_x (\pi^a T^a) V_x^\dagger

    //remember that the gauge link on the G-parity boundary transforms as V U V^T
    
    Lattice &lat = *lattice;
    //test the gauge field evolution by calculating the conjugate momentum for the U links
    //and separately for the U* links and compare them
    
    int f_sites = GJP.SnodeSites()*GJP.VolNodeSites() / (lat.FchkbEvl()+1); //half-checkerboard fermion

    int f_vec_count = 2* f_sites * lat.SpinComponents(); //number of 3-vectors in fermion
    size_t f_size = f_vec_count * lat.Colors() * 2;

    int Ncb;
    if(lat.FchkbEvl() == 1) Ncb = 1; 
    else if(lat.FchkbEvl() == 0) Ncb = 2;

    Vector* tmp1 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1","","");
    Vector* tmp2 = (Vector*)smalloc(f_size*sizeof(Float),"tmp2","","");
    Vector* phi = (Vector *) smalloc(f_size*sizeof(Float),"phi","","");
    //~~make a random two-flavour pseudofermion field from the heatbath as usual~~
    LatRanGen LRGbak(LRG);
    lat.RandGaussVector(tmp1, 0.5, Ncb);
    lat.RandGaussVector(tmp2, 0.5, Ncb);
    LRG = LRGbak;

    Float h_init;
    // phi <- M_f^\dag (RGV)
    h_init = lat.SetPhi(phi, tmp1, tmp2, 0.2, DAG_YES); //mass is 0.2
    printf("h_init = %f\n", h_init);

    //~~mom field evolve~~

    //invert D^\dag D on phi
    Vector *cg_sol = (Vector*)smalloc(f_size*sizeof(Float),"cg_sol","","");
    CgArg frm_cg_arg_md;
    frm_cg_arg_md.mass = 0.2;
    frm_cg_arg_md.max_num_iter = 5000;
    frm_cg_arg_md.stop_rsd = 1e-08;

    lat.FmatEvlInv(cg_sol, phi, &frm_cg_arg_md, CNV_FRM_NO);

    int g_size = 2*GJP.VolNodeSites() * lat.GsiteSize();

    //take a copy of the gauge field
    Matrix* gcopy = (Matrix*)smalloc(g_size*sizeof(Float),"gcopy","","");
    memcpy((void*)gcopy,(void*)lat.GaugeField(),g_size*sizeof(Float));

    //evolve U conj mom then copy-conjugate to U* field conj mom
    Matrix* mom = (Matrix*)smalloc(g_size*sizeof(Float),"mom","","");
    for(int i=0;i<2*4*GJP.VolNodeSites();i++) mom[i].ZeroMatrix();

    lat.EvolveMomFforce(mom, cg_sol, 0.5, 1.0); //evolve 

    //restore the gauge field
    memcpy((void*)lat.GaugeField(),(void*)gcopy,g_size*sizeof(Float));
    
    //generate a random gauge transformation
    Matrix *gtrans = (Matrix *) pmalloc(GJP.VolNodeSites()*18*sizeof(Float));
#if 1
    for(int i=0;i<GJP.VolNodeSites();i++){
      LRG.AssignGenerator(i);
      Matrix *m = gtrans+i;
      IFloat * m0 = (IFloat*)m;
      for(int j=0;j<18;j++){
	*(m0+j) = LRG.Urand();
      }
      m->Unitarize(); //make SU(3)
    }
#endif

#if 0
    //unit matrix test
    for(int i=0;i<GJP.VolNodeSites();i++){
      Matrix *m = gtrans+i;
      IFloat * m0 = (IFloat*)m;
      for(int j=0;j<18;j++){
	*(m0+j) = 0.0;
      }
      m0[0]=1.0;
      m0[8]=1.0;
      m0[16]=1.0;
    }
#endif
    GaugeTransformP(gtrans,mom,lat); //for later comparison

    GaugeTransformU(gtrans,lat);
    GaugeTransformF(gtrans,tmp1,lat);
    GaugeTransformF(gtrans,tmp2,lat);

    // redo phi  phi <- M_f^\dag (RGV) and recalculate cg_sol = D^\dag D on phi
    h_init = lat.SetPhi(phi, tmp1, tmp2, 0.2, DAG_YES); //mass is 0.2
    printf("h_init post gauge-transform = %f\n", h_init);

    lat.FmatEvlInv(cg_sol, phi, &frm_cg_arg_md, CNV_FRM_NO);

    //evolve U conj mom then copy-conjugate to U* field conj mom
    Matrix* mom2 = (Matrix*)smalloc(g_size*sizeof(Float),"mom2","","");
    for(int i=0;i<2*4*GJP.VolNodeSites();i++) mom2[i].ZeroMatrix();

    lat.EvolveMomFforce(mom2, cg_sol, 0.5, 1.0); //evolve 

    //restore the gauge field
    memcpy((void*)lat.GaugeField(),(void*)gcopy,g_size*sizeof(Float));

    int flav_off = GJP.VolNodeSites()*4;

    //mom and mom2 should be the same
    Float err(0.0);
    
    for(int f=0;f<2;f++){
      for(int t=0;t<GJP.TnodeSites();t++){
	for(int z=0;z<GJP.ZnodeSites();z++){
	  for(int y=0;y<GJP.YnodeSites();y++){
	    for(int x=0;x<GJP.XnodeSites();x++){
	      int pos[4] = {x,y,z,t};
	      for(int mu=0;mu<4;mu++){
		int off = lat.GsiteOffset(pos) + mu + f*flav_off;
		Float* m = (Float*)(mom+off);
		Float* mc = (Float*)(mom2+off);
		if(fabs(*m - *mc) > 1e-06 || fabs(*(m+1) - *(mc+1)) > 1e-06 ){
		  printf("Node %d %d %d %d , Error: f%d (%d %d %d %d), %d: (%f %f), (%f %f)\n",GJP.XnodeCoor(),GJP.YnodeCoor(),GJP.ZnodeCoor(),GJP.TnodeCoor(),f,x,y,z,t,mu,*m,*(m+1),*mc, *(mc+1));
		  err=1.0;
		}
	      }
	    }
	  }
	}
      }
    }
    glb_sum(&err);

    if(err>0.0){
      printf("Gauge transform check failed\n");
      exit(-1);
    }

    sfree(tmp1);
    sfree(tmp2);
    sfree(phi);
    sfree(cg_sol);
    sfree(gcopy);
    sfree(mom);
    sfree(mom2);
    printf("Gauge transform check passed\n");
  }
#endif

#define _1F_2F_COMPARISON
#ifdef _1F_2F_COMPARISON
  Matrix *_2f_mom;
  int _2f_norm;
  {
    printf("Starting 1f 2f evolve check\n");
    Lattice &lat = *lattice;
    //test the gauge field evolution by calculating the conjugate momentum in both the 2f and 1f models
    //and compare them

    int f_sites = GJP.SnodeSites()*GJP.VolNodeSites() / (lat.FchkbEvl()+1); //half-checkerboard fermion

    int f_vec_count = 2* f_sites * lat.SpinComponents(); //number of 3-vectors in fermion
    size_t f_size = f_vec_count * lat.Colors() * 2;

    int Ncb;
    if(lat.FchkbEvl() == 1) Ncb = 1; 
    else if(lat.FchkbEvl() == 0) Ncb = 2;

    Vector* tmp1 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1","","");
    Vector* tmp2 = (Vector*)smalloc(f_size*sizeof(Float),"tmp2","","");
    Vector* phi = (Vector *) smalloc(f_size*sizeof(Float),"phi","","");
    
    //backup the RNG
    LatRanGen LRGbak(LRG);

    //~~make a random two-flavour pseudofermion field from the heatbath as usual~~
    lat.RandGaussVector(tmp1, 0.5, Ncb);
    lat.RandGaussVector(tmp2, 0.5, Ncb);

    //restore the RNG
    LRG = LRGbak;

    Float h_init = 0.0;
    // phi <- M_f^\dag (RGV)
    h_init += lat.SetPhi(phi, tmp1, tmp2, 0.2, DAG_YES); //mass is 0.2

    //~~mom field evolve~~

    //invert D^\dag D on phi
    Vector *cg_sol = (Vector*)smalloc(f_size*sizeof(Float),"cg_sol","","");
    CgArg frm_cg_arg_md;
    frm_cg_arg_md.mass = 0.2;
    frm_cg_arg_md.max_num_iter = 5000;
    frm_cg_arg_md.stop_rsd = 1e-08;

    lat.FmatEvlInv(cg_sol, phi, &frm_cg_arg_md, CNV_FRM_NO);
    _2f_norm = vect_norm(cg_sol,f_vec_count);

    int g_size = 2*GJP.VolNodeSites() * lat.GsiteSize();

    //evolve U conj mom then copy-conjugate to U* field conj mom
    Matrix* mom = (Matrix*)smalloc(g_size*sizeof(Float),"mom","","");
    for(int i=0;i<2*4*GJP.VolNodeSites();i++) mom[i].ZeroMatrix();

    lat.EvolveMomFforce(mom, cg_sol, 0.5, 1.0); //evolve 

    _2f_mom = mom;

    sfree(tmp1);
    sfree(tmp2);
    sfree(phi);
    sfree(cg_sol);
  }
#endif


  if(gauge_fix) lattice->FixGaugeFree();

  if(UniqueID()==0){ printf("Starting double lattice force\n"); fflush(stdout); }

  
  int array_size = 2*lattice->GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
  Matrix *orig_lattice = (Matrix *) pmalloc(array_size);
  memcpy((void*)orig_lattice, (void*)lattice->GaugeField(), array_size);

  lattice->FreeGauge(); //free memory and reset
  delete lattice; //lattice objects are singleton (scope_lock)

  //setup 1f model. Upon calling GJP.Initialize the lattice size will be doubled in the appropriate directions
  //and the boundary condition set to APRD
  if(gparity_X) do_arg.gparity_1f_X = 1;
  if(gparity_Y) do_arg.gparity_1f_Y = 1;

  GJP.Initialize(do_arg);

  if(GJP.Gparity()){ printf("Que?\n"); exit(-1); }
  if(UniqueID()==0) printf("Doubled lattice : %d %d %d %d\n", GJP.XnodeSites()*GJP.Xnodes(),GJP.YnodeSites()*GJP.Ynodes(),
			   GJP.ZnodeSites()*GJP.Znodes(),GJP.TnodeSites()*GJP.Tnodes());
  
#ifdef HAVE_BFM
  {
    QDP::multi1d<int> nrow(Nd);  
    for(int i = 0;i<Nd;i++) nrow[i] = GJP.Sites(i);
    //  multi1d<LatticeFermion> test(Nd);  
    //  nrow=size;
    QDP::Layout::setLattSize(nrow);
    QDP::Layout::create();
  }
#endif

  GwilsonFdwf* doubled_lattice = new GwilsonFdwf;
  setup_double_latt(*doubled_lattice,orig_lattice,gparity_X,gparity_Y);
  setup_double_rng(gparity_X,gparity_Y);
 
  if(gauge_fix){
    doubled_lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    doubled_lattice->FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }

#ifdef _1F_2F_COMPARISON  
  {
    printf("Starting 1f 2f evolve check part II\n");
    Lattice &lat = *doubled_lattice;

    int g_size = GJP.VolNodeSites() * lat.GsiteSize();
    Matrix *_2f_mom_dbl = (Matrix*)smalloc(g_size*sizeof(Float),"mom","","");
    setup_double_matrixfield(_2f_mom_dbl, _2f_mom, 4, gparity_X,gparity_Y);

    //test the gauge field evolution by calculating the conjugate momentum in both the 2f and 1f models
    //and compare them

    int f_sites = GJP.SnodeSites()*GJP.VolNodeSites() / (lat.FchkbEvl()+1); //half-checkerboard fermion

    int f_vec_count = f_sites * lat.SpinComponents(); //number of 3-vectors in fermion
    size_t f_size = f_vec_count * lat.Colors() * 2;

    int Ncb;
    if(lat.FchkbEvl() == 1) Ncb = 1; 
    else if(lat.FchkbEvl() == 0) Ncb = 2;

    Vector* tmp1 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1","","");
    Vector* tmp2 = (Vector*)smalloc(f_size*sizeof(Float),"tmp2","","");
    Vector* phi = (Vector *) smalloc(f_size*sizeof(Float),"phi","","");
    
    //backup the RNG
    LatRanGen LRGbak(LRG);

    //~~make a random one-flavour pseudofermion field from the heatbath as usual~~
    lat.RandGaussVector(tmp1, 0.5, Ncb);
    lat.RandGaussVector(tmp2, 0.5, Ncb);

    if(GJP.Gparity1fY()){
      //make source on upper-right quadrant negative (RNGs should be correct)
      for(int s=0;s<GJP.SnodeSites();s++){
	for(int t=0;t<GJP.TnodeSites();t++){
	  for(int z=0;z<GJP.ZnodeSites();z++){
	    for(int y=0;y<GJP.YnodeSites();y++){
	      for(int x=0;x<GJP.XnodeSites();x++){
		if( (x+y+z+t+s)%2 == 0) continue; //ferm vect is odd parity

		int gx = x+GJP.XnodeCoor()*GJP.XnodeSites();
		int gy = y+GJP.YnodeCoor()*GJP.YnodeSites();

		//
		if(gx>=GJP.Xnodes()*GJP.XnodeSites()/2 && gy>=GJP.Ynodes()*GJP.YnodeSites()/2){
		  int pos[5] = {x,y,z,t,s};
		  int f_off = lat.FsiteOffsetChkb(pos) * lat.SpinComponents();

		  for(int spn=0;spn<lat.SpinComponents();spn++){
		    *(tmp1+f_off+spn) *=-1;
		    *(tmp2+f_off+spn) *=-1;
		  }

		}
	      }
	    }
	  }
	}
      }
    }

    //restore the RNG
    LRG = LRGbak;

    Float h_init = 0.0;
    // phi <- M_f^\dag (RGV)
    h_init += lat.SetPhi(phi, tmp1, tmp2, 0.2, DAG_YES); //mass is 0.2

    //~~mom field evolve~~

    //invert D^\dag D on phi
    Vector *cg_sol = (Vector*)smalloc(f_size*sizeof(Float),"cg_sol","","");
    CgArg frm_cg_arg_md;
    frm_cg_arg_md.mass = 0.2;
    frm_cg_arg_md.max_num_iter = 5000;
    frm_cg_arg_md.stop_rsd = 1e-08;

    lat.FmatEvlInv(cg_sol, phi, &frm_cg_arg_md, CNV_FRM_NO);
    int _1f_norm = vect_norm(cg_sol,f_vec_count);
    if(GJP.Gparity1fY()) _1f_norm/=2; //quad lattice has redundant degrees of freedom

    if( (_1f_norm - _2f_norm)/(_1f_norm + _2f_norm)*2 > 0.01 ){
      printf("Node %d, 2f and 1f prop norms not equal: %d %d\n", UniqueID(), _2f_norm, _1f_norm); exit(-1);
    }else{
      printf("Node %d, 2f and 1f prop norms are equal: %d %d\n", UniqueID(), _2f_norm, _1f_norm);
    }

    //evolve U conj mom then copy-conjugate to U* field conj mom
    Matrix* mom = (Matrix*)smalloc(g_size*sizeof(Float),"mom","","");
    for(int i=0;i<4*GJP.VolNodeSites();i++) mom[i].ZeroMatrix();

    lat.EvolveMomFforce(mom, cg_sol, 0.5, 1.0); //evolve 

    //check 1f vs 2f mom
    {
      bool err(false);
    
      for(int t=0;t<GJP.TnodeSites();t++){
	for(int z=0;z<GJP.ZnodeSites();z++){
	  for(int y=0;y<GJP.YnodeSites();y++){
	    for(int x=0;x<GJP.XnodeSites();x++){
	      int pos[4] = {x,y,z,t};
	      int gpos[4] = {x+GJP.XnodeCoor()*GJP.XnodeSites(),
			     y+GJP.YnodeCoor()*GJP.YnodeSites(),			     
			     z+GJP.ZnodeCoor()*GJP.ZnodeSites(),
			     t+GJP.TnodeCoor()*GJP.TnodeSites()};

	      for(int mu=0;mu<4;mu++){
		int off = lat.GsiteOffset(pos) + mu;
		Float* m = (Float*)(mom+off);
		Float* mc = (Float*)(_2f_mom_dbl+off);
		if(fabs(*m - *mc) > 1e-05 || fabs(*(m+1) - *(mc+1)) > 1e-05 ){
		  printf("Error: Xnode %d, 1f:2f (%d %d %d %d), %d: (%f %f), (%f %f)\n",GJP.XnodeCoor(),gpos[0],gpos[1],gpos[2],gpos[3],mu,*m,*(m+1),*mc, *(mc+1));
		  err=true;
		}
	      }
	    }
	  }
	}
      }
    
      if(err){
	printf("Failed 1f-2f comparison\n");
	exit(-1);
      }
      printf("Passed 1f-2f comparison\n");
    }


    sfree(tmp1);
    sfree(tmp2);
    sfree(phi);
    sfree(cg_sol);
  }
#endif


  if(gauge_fix) doubled_lattice->FixGaugeFree();
  doubled_lattice->FreeGauge(); //free memory and reset
  delete doubled_lattice; //lattice objects are singleton (scope_lock)

  //#define TEST_U_USTAR_FORCE_RELN
#ifdef TEST_U_USTAR_FORCE_RELN
  {
    if(gparity_X && gparity_Y){ printf("TEST_U_USTAR_FORCE_RELN test not valid for quad lattice\n"); exit(-1); }

    //generate a single lattice with random gauge links and APBC

    do_arg.x_bc = BND_CND_PRD;
    do_arg.x_sites = size[0];
    if(do_arg.x_sites/GJP.Xnodes() < 2){
      ERR.General("main","TEST_U_USTAR_FORCE_RELN","After dividing initial Xnodesites by 2, too few sites on node: %d\n",do_arg.x_sites/GJP.Xnodes());
    }

    GJP.Initialize(do_arg);

#ifdef HAVE_BFM
    {
      QDP::multi1d<int> nrow(Nd);  
      for(int i = 0;i<Nd;i++) nrow[i] = GJP.Sites(i);
      //  multi1d<LatticeFermion> test(Nd);  
      //  nrow=size;
      QDP::Layout::setLattSize(nrow);
      QDP::Layout::create();
    }
#endif
    GwilsonFdwf* lattice = new GwilsonFdwf;
    lattice->SetGfieldDisOrd();
  
    Lattice &lat = *lattice;
    int f_sites = GJP.SnodeSites()*GJP.VolNodeSites() / (lat.FchkbEvl()+1); //half-checkerboard fermion

    int f_vec_count = f_sites * lat.SpinComponents(); //number of 3-vectors in fermion
    size_t f_size = f_vec_count * lat.Colors() * 2;

    int Ncb;
    if(lat.FchkbEvl() == 1) Ncb = 1; 
    else if(lat.FchkbEvl() == 0) Ncb = 2;

    Vector* tmp1 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1","","");
    Vector* tmp2 = (Vector*)smalloc(f_size*sizeof(Float),"tmp2","","");
    Vector* phi = (Vector *) smalloc(f_size*sizeof(Float),"phi","","");
    //~~make a random pseudofermion field from the heatbath as usual~~
    lat.RandGaussVector(tmp1, 0.5, Ncb);
    lat.RandGaussVector(tmp2, 0.5, Ncb);

    Float h_init = 0.0;
    // phi <- M_f^\dag (RGV)
    h_init += lat.SetPhi(phi, tmp1, tmp2, 0.2, DAG_YES); //mass is 0.2

    //~~mom field evolve~~

    //invert D^\dag D on phi
    Vector *cg_sol = (Vector*)smalloc(f_size*sizeof(Float),"cg_sol","","");
    CgArg frm_cg_arg_md;
    frm_cg_arg_md.mass = 0.2;
    frm_cg_arg_md.max_num_iter = 5000;
    frm_cg_arg_md.stop_rsd = 1e-06;

    lat.FmatEvlInv(cg_sol, phi, &frm_cg_arg_md, CNV_FRM_NO);

    //evolve the momentum
    int g_size = GJP.VolNodeSites() * lat.GsiteSize();
    Matrix* mom = (Matrix*)smalloc(g_size*sizeof(Float),"mom","","");
    for(int i=0;i<4*GJP.VolNodeSites();i++) mom[i].ZeroMatrix();

    lat.EvolveMomFforce(mom, cg_sol, 0.5, 1.0); //evolve 

    //complex conjugate all links on lattice U->U*, U*->U
    Matrix *lp = lat.GaugeField();
    for(int t=0;t<GJP.TnodeSites();t++){
      for(int z=0;z<GJP.ZnodeSites();z++){
	for(int y=0;y<GJP.YnodeSites();y++){
	  for(int x=0;x<GJP.XnodeSites();x++){
	    int pos[4] = {x,y,z,t};
	    int off = lat.GsiteOffset(pos);
	    for(int mu=0;mu<4;mu++){
	      Float* m = (Float*)(lp+off+mu);
	      for(int i=1;i<18;i+=2) m[i] = -m[i];
	    }
	  }
	}
      }
    }


    Float *vec_base = (Float*)cg_sol;
    int vec_size = 2 * lat.Colors() * lat.SpinComponents();
    for(int i = 0; i < GJP.VolNodeSites()*GJP.SnodeSites(); i+=2) {
      WilsonVector* wvec = (WilsonVector*)vec_base;
      wvec->conj().ccl(1).gamma(-5);
      vec_base+=vec_size;
    }
    //evolve the momentum
    Matrix* momconj = (Matrix*)smalloc(g_size*sizeof(Float),"momconj","","");
    for(int i=0;i<4*GJP.VolNodeSites();i++) momconj[i].ZeroMatrix();

    lat.EvolveMomFforce(momconj, cg_sol, 0.5, 1.0); //evolve 

    //printf("h_init %f, h_init_conj %f\n",h_init,h_init_conj);

    //mom and momconj should be complex conjugates of each other
    for(int t=0;t<GJP.TnodeSites();t++){
      for(int z=0;z<GJP.ZnodeSites();z++){
	for(int y=0;y<GJP.YnodeSites();y++){
	  for(int x=0;x<GJP.XnodeSites();x++){
	    int pos[4] = {x,y,z,t};
	    for(int mu=0;mu<4;mu++){
	      int off = lat.GsiteOffset(pos) + mu;
	      Float* m = (Float*)(mom+off);
	      Float* mc = (Float*)(momconj+off);
	      printf("%d %d %d %d, %d: (%f %f), (%f %f)\n",x,y,z,t,mu,*m,*(m+1),*mc, *(mc+1));
	    }
	  }
	}
      }
    }


  }
#endif

  //#define TEST_SINGLE_LATT_GTRANS
#ifdef TEST_SINGLE_LATT_GTRANS
  {
    if(gparity_X && gparity_Y){ printf("TEST_U_USTAR_FORCE_RELN test not valid for quad lattice\n"); exit(-1); }
    printf("Starting gauge transform check for single lattice\n");
    //test the code by performing a gauge transformation and checking the conjugate momentum transforms in the correct way
    //U -> e^{i \pi^a T^a}U
    //V_x U^mu_x V_{x+\mu}^\dagger -> V_x e^{i \pi^a T^a} U^\mu_x  V_{x+\mu}^\dagger = V_x(1 + [i\pi^aT^a] + [i\pi^aT^a][i\pi^bT^b]/2 +....)U^\mu_x V^\dagger_{x+\mu}
    //                                                = V_x(1 + [i\pi^aT^a] + [i\pi^aT^a]V^\dagger V[i\pi^bT^b]/2 +....)V_x^\dagger V_x U^\mu_x V^\dagger_{x+\mu}
    //                                                = e^{iV_x\pi^a T^a V_x^\dagger} V_x U V_{x+\mu}^\dagger
    // so (\pi^a T^a) -> V_x (\pi^a T^a) V_x^\dagger

    do_arg.x_bc = BND_CND_PRD;
    do_arg.x_sites = size[0];
    if(do_arg.x_sites/GJP.Xnodes() < 2){
      ERR.General("main","TEST_SINGLE_LATT_GTRANS","After dividing initial Xnodesites by 2, too few sites on node: %d\n",do_arg.x_sites/GJP.Xnodes());
    }
    GJP.Initialize(do_arg);

#ifdef HAVE_BFM
    {
      QDP::multi1d<int> nrow(Nd);  
      for(int i = 0;i<Nd;i++) nrow[i] = GJP.Sites(i);
      //  multi1d<LatticeFermion> test(Nd);  
      //  nrow=size;
      QDP::Layout::setLattSize(nrow);
      QDP::Layout::create();
    }
#endif
    GwilsonFdwf* lattice = new GwilsonFdwf;
    lattice->SetGfieldDisOrd();

    Lattice &lat = *lattice;
    //test the gauge field evolution by calculating the conjugate momentum for the U links
    //and separately for the U* links and compare them
    
    int f_sites = GJP.SnodeSites()*GJP.VolNodeSites() / (lat.FchkbEvl()+1); //half-checkerboard fermion

    int f_vec_count = f_sites * lat.SpinComponents(); //number of 3-vectors in fermion
    size_t f_size = 2*f_vec_count * lat.Colors();

    int Ncb;
    if(lat.FchkbEvl() == 1) Ncb = 1; 
    else if(lat.FchkbEvl() == 0) Ncb = 2;

    Vector* tmp1 = (Vector*)smalloc(f_size*sizeof(Float),"tmp1","","");
    Vector* tmp2 = (Vector*)smalloc(f_size*sizeof(Float),"tmp2","","");
    Vector* phi = (Vector *) smalloc(f_size*sizeof(Float),"phi","","");
    //~~make a random two-flavour pseudofermion field from the heatbath as usual~~
    lat.RandGaussVector(tmp1, 0.5, Ncb);
    lat.RandGaussVector(tmp2, 0.5, Ncb);

    Float h_init;
    // phi <- M_f^\dag (RGV)
    h_init = lat.SetPhi(phi, tmp1, tmp2, 0.2, DAG_YES); //mass is 0.2
    printf("h_init = %f\n", h_init);

    //~~mom field evolve~~

    //invert D^\dag D on phi
    Vector *cg_sol = (Vector*)smalloc(f_size*sizeof(Float),"cg_sol","","");
    CgArg frm_cg_arg_md;
    frm_cg_arg_md.mass = 0.2;
    frm_cg_arg_md.max_num_iter = 5000;
    frm_cg_arg_md.stop_rsd = 1e-08;

    lat.FmatEvlInv(cg_sol, phi, &frm_cg_arg_md, CNV_FRM_NO);

    int g_size = GJP.VolNodeSites() * lat.GsiteSize();

    //calculate U conj mom
    Matrix* mom = (Matrix*)smalloc(g_size*sizeof(Float),"mom","","");
    for(int i=0;i<4*GJP.VolNodeSites();i++) mom[i].ZeroMatrix();

    lat.EvolveMomFforce(mom, cg_sol, 0.5, 1.0); //evolve 
    
    //generate a random gauge transformation
    Matrix *gtrans = (Matrix *) pmalloc(GJP.VolNodeSites()*18*sizeof(Float));
#if 1
    for(int i=0;i<GJP.VolNodeSites();i++){
      LRG.AssignGenerator(i);
      Matrix *m = gtrans+i;
      IFloat * m0 = (IFloat*)m;
      for(int j=0;j<18;j++){
	*(m0+j) = LRG.Urand();
      }
      m->Unitarize(); //make SU(3)
    }
#endif

#if 0
    //unit matrix test
    for(int i=0;i<GJP.VolNodeSites();i++){
      Matrix *m = gtrans+i;
      IFloat * m0 = (IFloat*)m;
      for(int j=0;j<18;j++){
	*(m0+j) = 0.0;
      }
      m0[0]=1.0;
      m0[8]=1.0;
      m0[16]=1.0;
    }
#endif

    //take a copy of the gauge field
    Matrix* ucopy = (Matrix*)smalloc(g_size*sizeof(Float),"ucopy","","");
    memcpy((void*)ucopy,(void*)lat.GaugeField(),g_size*sizeof(Float));

    //take a copy of the original mom field
    Matrix* momcopy = (Matrix*)smalloc(g_size*sizeof(Float),"momcopy","","");
    memcpy((void*)momcopy,(void*)mom,g_size*sizeof(Float));

    GaugeTransformU(gtrans,lat);
    GaugeTransformP(gtrans,mom,lat);

    Matrix* uprimecopy = (Matrix*)smalloc(g_size*sizeof(Float),"uprimecopy","","");
    memcpy((void*)uprimecopy,(void*)lat.GaugeField(),g_size*sizeof(Float));

    //test transform of mom: e^{i\pi}U should give the same as e^{i\pi'}U' where prime indicates gauge transform performed
    Matrix* evolved = (Matrix*)smalloc(g_size*sizeof(Float),"evolved1","","");
    Matrix* evolvedprime = (Matrix*)smalloc(g_size*sizeof(Float),"evolved2","","");

    //do primed first
    lat.EvolveGfield(mom,1);
    memcpy((void*)evolvedprime,(void*)lat.GaugeField(),g_size*sizeof(Float));

    //restore the original gauge field and do unprimed
    memcpy((void*)lat.GaugeField(),(void*)ucopy,g_size*sizeof(Float));
    lat.EvolveGfield(momcopy,1);
    GaugeTransformU(gtrans,lat);
    memcpy((void*)evolved,(void*)lat.GaugeField(),g_size*sizeof(Float));

    {
      bool err(false);
    
      for(int t=0;t<GJP.TnodeSites();t++){
	for(int z=0;z<GJP.ZnodeSites();z++){
	  for(int y=0;y<GJP.YnodeSites();y++){
	    for(int x=0;x<GJP.XnodeSites();x++){
	      int pos[4] = {x,y,z,t};
	      for(int mu=0;mu<4;mu++){
		int off = lat.GsiteOffset(pos) + mu;
		Float* m = (Float*)(evolved+off);
		Float* mc = (Float*)(evolvedprime+off);
		if(fabs(*m - *mc) > 1e-06 || fabs(*(m+1) - *(mc+1)) > 1e-06 ){
		  printf("Error: U_new (%d %d %d %d), %d: (%f %f), (%f %f)\n",x,y,z,t,mu,*m,*(m+1),*mc, *(mc+1));
		  err=true;
		}
	      }
	    }
	  }
	}
      }
    
      if(err) exit(-1);
      printf("Passed test of mom field gtrans relation\n");
      pfree(evolved);
      pfree(evolvedprime);
    }
    //----prepare to recalculated conjugate momentum for transformed field------

    //restore transformed gauge field
    memcpy((void*)lat.GaugeField(),(void*)uprimecopy,g_size*sizeof(Float));

    //gauge transform the fields
    GaugeTransformF(gtrans,tmp1,lat);
    GaugeTransformF(gtrans,tmp2,lat);

    // redo phi  phi <- M_f^\dag (RGV) and recalculate cg_sol = D^\dag D on phi
    h_init = lat.SetPhi(phi, tmp1, tmp2, 0.2, DAG_YES); //mass is 0.2
    printf("h_init post gauge-transform = %f\n", h_init);

    lat.FmatEvlInv(cg_sol, phi, &frm_cg_arg_md, CNV_FRM_NO);

    //evolve U conj mom then copy-conjugate to U* field conj mom
    Matrix* mom2 = (Matrix*)smalloc(g_size*sizeof(Float),"mom2","","");
    for(int i=0;i<4*GJP.VolNodeSites();i++) mom2[i].ZeroMatrix();

    lat.EvolveMomFforce(mom2, cg_sol, 0.5, 1.0); //evolve 

    //mom and mom2 should be the same
    bool err(false);
    
    for(int t=0;t<GJP.TnodeSites();t++){
      for(int z=0;z<GJP.ZnodeSites();z++){
	for(int y=0;y<GJP.YnodeSites();y++){
	  for(int x=0;x<GJP.XnodeSites();x++){
	    int pos[4] = {x,y,z,t};
	    for(int mu=0;mu<4;mu++){
	      int off = lat.GsiteOffset(pos) + mu;
	      Float* m = (Float*)(mom+off);
	      Float* mc = (Float*)(mom2+off);
	      if(fabs(*m - *mc) > 1e-06 || fabs(*(m+1) - *(mc+1)) > 1e-06 ){
		printf("Error: (%d %d %d %d), %d: (%f %f), (%f %f)\n",x,y,z,t,mu,*m,*(m+1),*mc, *(mc+1));
		err=true;
	      }
	    }
	  }
	}
      }
    }
    
    if(err) exit(-1);


    sfree(tmp1);
    sfree(tmp2);
    sfree(phi);
    sfree(cg_sol);
    sfree(ucopy);
    sfree(mom);
    sfree(mom2);
    sfree(gtrans);
    printf("Single lattice Gauge transform check passed\n");
  }
#endif



#ifdef HAVE_BFM
  Chroma::finalize();
#endif

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}

void setup_double_latt(Lattice &double_latt, Matrix* orig_gfield, bool gparity_X, bool gparity_Y){
  //orig latt ( U_0 U_1 ) ( U_2 U_3 ) ( U_4 U_5 ) ( U_6 U_7 )
  //double tatt ( U_0 U_1 U_2 U_3 ) ( U_4 U_5 U_6 U_7 ) ( U_0* U_1* U_2* U_3* ) ( U_4* U_5* U_6* U_7* )

  Matrix *dbl_gfield = double_latt.GaugeField();

  if(!UniqueID()){ printf("Setting up 1f lattice.\n"); fflush(stdout); }
  SingleToDoubleLattice lattdoubler(gparity_X,gparity_Y,orig_gfield,double_latt);
  lattdoubler.Run();
  if(!UniqueID()){ printf("Finished setting up 1f lattice\n"); fflush(stdout); }
}
void setup_double_rng(bool gparity_X, bool gparity_Y){
  //orig 4D rng 2 stacked 4D volumes
  //orig ([R_0 R_1][R'_0 R'_1])([R_2 R_3][R'_2 R'_3])([R_4 R_5][R'_4 R'_5])([R_6 R_7][R'_6 R'_7])
  //double (R_0 R_1 R_2 R_3)(R_4 R_5 R_6 R_7)(R'_0 R'_1 R'_2 R'_3)(R'_4 R'_5 R'_6 R'_7)
  
  //orig 5D rng 2 stacked 4D volumes per ls/2 slice (ls/2 as only one RNG per 2^4 block)

  SingleToDouble4dRNG fourDsetup(gparity_X,gparity_Y);
  SingleToDouble5dRNG fiveDsetup(gparity_X,gparity_Y);
  
  LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form
  
  if(!UniqueID()){ printf("Setting up 1f 4D RNG\n"); fflush(stdout); }
  fourDsetup.Run();      
  if(!UniqueID()){ printf("Setting up 1f 5D RNG\n"); fflush(stdout); }
  fiveDsetup.Run();    
}

void setup_double_matrixfield(Matrix* double_mat, Matrix* orig_mat, int nmat_per_site, bool gparity_X, bool gparity_Y){
  if(!UniqueID()){ printf("Setting up 1f matrix field.\n"); fflush(stdout); }
  SingleToDoubleMatrixField doubler(gparity_X,gparity_Y,nmat_per_site,orig_mat,double_mat);
  doubler.Run();
  if(!UniqueID()){ printf("Finished setting up 1f matrixfield\n"); fflush(stdout); }
}

#include <util/dirac_op.h>


void EvolveMomFforceUstarGparity(Lattice &lat, Matrix *mom, Vector *chi, double mass, double dt){
  printf("Test version\n");

  char* cname = "";
  char* fname = "EvolveMomFforceUstarGparity";

  if(!GJP.Gparity()) ERR.General(cname,fname,"Function executed when G-parity BCs not in use\n");

  Matrix *gauge = lat.GaugeField() ;

  if (lat.Colors() != 3)
    ERR.General(cname,fname,"Wrong nbr of colors.") ;
 
  if (lat.SpinComponents() != 4)
    ERR.General(cname,fname,"Wrong nbr of spin comp.") ;
 
  if (mom == 0)
    ERR.Pointer(cname,fname,"mom") ;
 
  if (chi == 0)
    ERR.Pointer(cname,fname,"chi") ;
 
  //----------------------------------------------------------------
  // allocate space for two CANONICAL fermion fields
  //----------------------------------------------------------------

  size_t f_size = lat.FsiteSize() * GJP.VolNodeSites() ;
  int f_site_size_4d = 2 * lat.Colors() * lat.SpinComponents();
  size_t f_size_4d = f_site_size_4d * GJP.VolNodeSites() ;

  //CK: Need space for both d and C\bar u^T fields stacked
  //Two 4d volumes are stacked on each Ls such that Ls can be split over multiple nodes
  size_t f_size_alloc = f_size *2 *sizeof(double);
  int f_single4dsite_alloc = lat.FsiteSize()*sizeof(double)*2;


  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(f_size_alloc) ;
  if (v1 == 0) ERR.Pointer(cname, fname, str_v1) ;
  VRB.Smalloc(cname, fname, str_v1, v1, f_size_alloc) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(f_size_alloc) ;
  if (v2 == 0) ERR.Pointer(cname, fname, str_v2) ;
  VRB.Smalloc(cname, fname, str_v2, v2, f_size_alloc) ;

  //----------------------------------------------------------------
  // allocate buffer space for two fermion fields that are assoc
  // with only one 4-D site.
  //----------------------------------------------------------------

  char *str_site_v1 = "site_v1" ;
  double *site_v1 = (double *)smalloc(f_single4dsite_alloc) ;
  if (site_v1 == 0) ERR.Pointer(cname, fname, str_site_v1) ;
  VRB.Smalloc(cname, fname, str_site_v1, site_v1, f_single4dsite_alloc) ;

  char *str_site_v2 = "site_v2" ;
  double *site_v2 = (double *)smalloc(f_single4dsite_alloc) ;
  if (site_v2 == 0) ERR.Pointer(cname, fname, str_site_v2) ;
  VRB.Smalloc(cname, fname, str_site_v2, site_v2, f_single4dsite_alloc) ;

  //CK: G-parity site vectors on each Ls are 2 * flavours 4 spin * 3 clr * 2 re/im
  //Sticking with same layout as 5d guys, stack 2 site vectors on each Ls
  //Ls=0 [ d ][ CubarT ] Ls=1 [ d ][ CubarT ] .....

  //----------------------------------------------------------------
  // Calculate v1, v2. Both v1, v2 must be in CANONICAL order after
  // the calculation.
  //----------------------------------------------------------------  

  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;

    DiracOpDwf dwf(lat, v1, v2, &cg_arg, CNV_FRM_YES) ;
    dwf.CalcHmdForceVecs(chi) ;  
  }

  int mu, x, y, z, t, s, lx, ly, lz, lt, ls ;
 
  lx = GJP.XnodeSites() ;
  ly = GJP.YnodeSites() ;
  lz = GJP.ZnodeSites() ;
  lt = GJP.TnodeSites() ;
  ls = GJP.SnodeSites() ;

  Matrix tmp_mat1, tmp_mat2, tmp_mat3;

  for (mu=0; mu<4; mu++){
    for (t=0; t<lt; t++){
    for (z=0; z<lz; z++){
    for (y=0; y<ly; y++){
    for (x=0; x<lx; x++){
      int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
      int vec_offset = f_site_size_4d*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      int pos[] = {x,y,z,t};
      int lattsz[] = {lx,ly,lz,lt};    

      double *v1_plus_mu ;
      double *v2_plus_mu ;

      double *v1_plus_mu_CubarT;
      double *v2_plus_mu_CubarT;      

      int vec_plus_mu_stride ;
      int vec_plus_mu_offset = f_site_size_4d ;

      double d_coeff = -2.0 * dt ; //coefficient for first G-parity field
      double CubarT_coeff = d_coeff; //coefficient for second G-parity field

      /* CK: Replaced nasty 300000 line switch statement with the following:*/
      {
	int pos_p_mu[] = {x,y,z,t};
	pos_p_mu[mu] = (pos_p_mu[mu]+1) % lattsz[mu];
	vec_plus_mu_offset *= pos_p_mu[0] + lx*(pos_p_mu[1] + ly*(pos_p_mu[2] + lz*pos_p_mu[3]));
      }      

      if ((pos[mu]+1) == lattsz[mu]) {
	if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu)==GJP.Nodes(mu)-1){
	  //at global boundary
	  // d <- +CubarT,  CubarT <- -d
	  CubarT_coeff = -CubarT_coeff;
	}else if(GJP.NodeBc(mu)==BND_CND_APRD){
	  // d <- -d,   CubarT <- -CubarT
	  d_coeff = -d_coeff;
	  CubarT_coeff = -CubarT_coeff;
	}

	for (s=0; s<ls; s++) {
	  //if on global lattice boundary we must perform the G-parity flavor swap
	  int d_site_vect_offset = s*f_site_size_4d*2; //where d-fields are stored in site vect
	  int CubarT_site_vect_offset = d_site_vect_offset + f_site_size_4d; //same for CubarT field (offset by 1 4d vect from d)
	      
	  //normally d <- d,  CubarT <- CubarT
	  int d_comms_to_offset = d_site_vect_offset; //place in site vect where comm'd d-field will be stored
	  int CubarT_comms_to_offset = CubarT_site_vect_offset; //....

	  if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu)==GJP.Nodes(mu)-1){
	    //at global boundary
	    // d <- +CubarT,  CubarT <- -d
	    CubarT_comms_to_offset = d_site_vect_offset;
	    d_comms_to_offset = CubarT_site_vect_offset;
	  }

	  //communicate d field
	  getPlusData( (double *)site_v1+d_comms_to_offset,
		       (double *)v1+vec_plus_mu_offset+s*f_size_4d*2,
		       f_site_size_4d, mu) ;
	  getPlusData( (double *)site_v2+d_comms_to_offset,
		       (double *)v2+vec_plus_mu_offset+s*f_size_4d*2,
		       f_site_size_4d, mu) ;
	  //communicate CubarT field
	  getPlusData( (double *)site_v1+CubarT_comms_to_offset,
		       (double *)v1 + vec_plus_mu_offset + f_size_4d + s*f_size_4d*2,
		       f_site_size_4d, mu) ;
	  getPlusData( (double *)site_v2+CubarT_comms_to_offset,
		       (double *)v2 + vec_plus_mu_offset + f_size_4d + s*f_size_4d*2,
		       f_site_size_4d, mu) ;
	} // end for s
	v1_plus_mu = site_v1 ;
	v2_plus_mu = site_v2 ;

	v1_plus_mu_CubarT = site_v1 + f_site_size_4d;
	v2_plus_mu_CubarT = site_v2 + f_site_size_4d;

	vec_plus_mu_stride = f_site_size_4d ; //skip over other stacked field
      } else {
	v1_plus_mu = (double *)v1+vec_plus_mu_offset ;
	v2_plus_mu = (double *)v2+vec_plus_mu_offset ;

	v1_plus_mu_CubarT = (double *)v1+vec_plus_mu_offset + f_size_4d ;
	v2_plus_mu_CubarT = (double *)v2+vec_plus_mu_offset + f_size_4d ;

	vec_plus_mu_stride = 2*f_size_4d - f_site_size_4d ;
      }

      // if(gauge_offset == 4){
      // 	printf("f1:\n");

      // 	printf("v1(x):\n");
      // 	double* m = (double*)v1 + vec_offset;
      // 	for(int i=0;i<12;i++){
      // 	  printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	}
      // 	printf("\n");
	
      // 	printf("v1(x+mu):\n");
      // 	m = v1_plus_mu;
      // 	for(int i=0;i<12;i++){
      // 	  printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	}
      // 	printf("\n");

      // 	printf("v2:\n");
      // 	m = (double*)v2 + vec_offset;
      // 	for(int i=0;i<12;i++){
      // 	  printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	}
      // 	printf("\n");

      // 	printf("v2(x+mu):\n");
      // 	m = v2_plus_mu;
      // 	for(int i=0;i<12;i++){
      // 	  printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	}
      // 	printf("\n");

      // 	printf("f2:\n");

      // 	printf("v1(x):\n");
      // 	m = (double*)v1 + vec_offset +f_size_4d;
      // 	for(int i=0;i<12;i++){
      // 	  printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	}
      // 	printf("\n");
	
      // 	printf("v1(x+mu):\n");
      // 	m = v1_plus_mu_CubarT;
      // 	for(int i=0;i<12;i++){
      // 	  printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	}
      // 	printf("\n");

      // 	printf("v2:\n");
      // 	m = (double*)v2 + vec_offset +f_size_4d;
      // 	for(int i=0;i<12;i++){
      // 	  printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	}
      // 	printf("\n");

      // 	printf("v2(x+mu):\n");
      // 	m = v2_plus_mu_CubarT;
      // 	for(int i=0;i<12;i++){
      // 	  printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	}
      // 	printf("\n");

      // }

      GwilsonFdwf* lattice = dynamic_cast<GwilsonFdwf*>(&lat);

      // if(gauge_offset == 4){ //site is x=y=z=t=mu=0
      // 	double v[24];
      // 	double w[24];
      // 	for(int i=0;i<24;i++){
      // 	  v[i] = LRG.Urand();
      // 	  w[i] = LRG.Urand();
      // 	}
      // 	Matrix test1;
      // 	lattice->sproj_tr[mu]( (double *)&test1,
      // 			       v,w,
      // 			       1, 0,0) ;
	

      // 	printf("test a:\n");
      // 	double *m = (double *)&test1;
      // 	for(int i=0;i<3;i++){
      // 	  for(int j=0;j<3;j++){ 
      // 	    printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	  }
      // 	  printf("\n");
      // 	}

      // 	lattice->sproj_tr[mu]( (double *)&test1,
      // 			       w,v,
      // 			       1, 0,0) ;


      // 	printf("test b:\n");
      // 	m = (double *)&test1;
      // 	for(int i=0;i<3;i++){
      // 	  for(int j=0;j<3;j++){ 
      // 	    printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	  }
      // 	  printf("\n");
      // 	}


	// //tmp_mat1_{ij} = tr_s ( P_+^mu v1(x+mu)_i v2^\dagger(x)_j 
	// lattice->sproj_tr[mu]( (double *)&test1,
	// 	      (double *)v1_plus_mu,
	// 	      (double *)v2+vec_offset,
	// 	      1, vec_plus_mu_stride, 2*f_size_4d-f_site_size_4d) ;

	// if(gauge_offset == 4){
	//   printf("f1 contrib a normal s1:\n");
	//   double *m = (double *)&test1;
	//   for(int i=0;i<3;i++){
	//     for(int j=0;j<3;j++){ 
	//       printf("(%f,%f) ",*m,*(m+1)); m+=2;
	//     }
	//     printf("\n");
	//   }
	// }

	// lattice->sproj_tr[mu]( (double *)&test1,
	// 		       (double *)v2+vec_offset,
	// 		       (double *)v1_plus_mu,
	// 		       1, vec_plus_mu_stride, 2*f_size_4d-f_site_size_4d) ;

	// if(gauge_offset == 4){
	//   printf("f1 contrib a s1:\n");
	//   double *m = (double *)&tmp_mat1;
	//   for(int i=0;i<3;i++){
	//     for(int j=0;j<3;j++){ 
	//       printf("(%f,%f) ",*m,*(m+1)); m+=2;
	//     }
	//     printf("\n");
	//   }
	// }
      // }

      
      //Calculate contributions to force at site
      //from d field

      //d field part 1
      //normally tmp_mat1_{ij} = tr_s ( P_+^mu v1(x+mu)_i v2^\dagger(x)_j    with v1 = X,  v2 = M X
      //for conjugate instead we want tmp_mat1_{ij} = tr_s ( P_+^mu v2(x)_i v1^\dagger(x+mu)_j

      lattice->sproj_tr[mu]( (double *)&tmp_mat1,
			     (double *)v2+vec_offset,
			     (double *)v1_plus_mu,
			     ls, 2*f_size_4d-f_site_size_4d, vec_plus_mu_stride) ;

      // if(gauge_offset == 4){
      // 	printf("f1 contrib a:\n");
      // 	double *m = (double *)&tmp_mat1;
      // 	for(int i=0;i<3;i++){
      // 	  for(int j=0;j<3;j++){ 
      // 	    printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	  }
      // 	  printf("\n");
      // 	}
      // }

      //d field part 2
      //normally tmp_mat2_{ij} = tr_s ( P_-^mu v2(x+mu)_i v1^\dagger(x)_j 
      //for conjugate instead we want   tmp_mat2_{ij} = tr_s ( P_-^mu v1(x)_i v2^\dagger(x+mu)_j 

      lattice->sproj_tr[mu+4]( (double *)&tmp_mat2,
			       (double *)v1+vec_offset,
			       (double *)v2_plus_mu,
			       ls, 2*f_size_4d-f_site_size_4d, vec_plus_mu_stride) ;

      // if(gauge_offset == 4){
      // 	printf("f1 contrib b:\n");
      // 	double *m = (double *)&tmp_mat2;
      // 	for(int i=0;i<3;i++){
      // 	  for(int j=0;j<3;j++){ 
      // 	    printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	  }
      // 	  printf("\n");
      // 	}
      // }

      tmp_mat1 += tmp_mat2;
      tmp_mat2.Trans(tmp_mat1);//must transpose outer product on colour index
      tmp_mat1 = tmp_mat2;
      tmp_mat1 *=d_coeff; //different coefficient for the two flavours

      // if(gauge_offset == 4){
      // 	printf("f1 contrib:\n");
      // 	double *m = (double *)&tmp_mat1;
      // 	for(int i=0;i<3;i++){
      // 	  for(int j=0;j<3;j++){ 
      // 	    printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	  }
      // 	  printf("\n");
      // 	}
      // }
      

      //Cubar^T field part 1
      //normally tmp_mat2_{ij} = tr_s ( P_-^mu v1(x)_i v2^\dagger(x+mu)_j 
      //for conjugate we want   tmp_mat2_{ij} = tr_s ( P_-^mu v2(x+mu)_i v1^\dagger(x)_j
      
      lattice->sproj_tr[mu+4]( (double *)&tmp_mat2,
			     (double *)v2_plus_mu_CubarT,
			     (double *)v1+vec_offset+f_size_4d,
			     ls, vec_plus_mu_stride, 2*f_size_4d-f_site_size_4d) ;

      // if(gauge_offset == 4){
      // 	printf("f2 contrib a:\n");
      // 	double *m = (double *)&tmp_mat2;
      // 	for(int i=0;i<3;i++){
      // 	  for(int j=0;j<3;j++){ 
      // 	    printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	  }
      // 	  printf("\n");
      // 	}
      // }

      //Cubar^T field part 2
      //tmp_mat3_{ij} = tr_s ( P_+^mu v2(x)_i v1^\dagger(x+mu)_j
      //for conjugate we want   tmp_mat3_{ij} = tr_s ( P_+^mu v1(x+mu)_i v2^\dagger(x)_j
      
      lattice->sproj_tr[mu]( (double *)&tmp_mat3,
			     (double *)v1_plus_mu_CubarT,
			     (double *)v2+vec_offset+f_size_4d,
			     ls, vec_plus_mu_stride, 2*f_size_4d-f_site_size_4d) ;

      // if(gauge_offset == 4){
      // 	printf("f2 contrib b:\n");
      // 	double *m = (double *)&tmp_mat3;
      // 	for(int i=0;i<3;i++){
      // 	  for(int j=0;j<3;j++){ 
      // 	    printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	  }
      // 	  printf("\n");
      // 	}
      // }

      tmp_mat2+= tmp_mat3;
      tmp_mat2*=CubarT_coeff;

      tmp_mat1 += tmp_mat2;

      // if(gauge_offset == 4){
      // 	printf("f2 contrib:\n");
      // 	double *m = (double *)&tmp_mat2;
      // 	for(int i=0;i<3;i++){
      // 	  for(int j=0;j<3;j++){ 
      // 	    printf("(%f,%f) ",*m,*(m+1)); m+=2;
      // 	  }
      // 	  printf("\n");
      // 	}
      // }

      // If GJP.Snodes > 1 sum up contributions from all s nodes
      if(GJP.Snodes() > 1) {
	glb_sum_multi_dir((double *)&tmp_mat1,4,sizeof(Matrix)/sizeof(double));
      }

      tmp_mat2.DotMEqual(*(gauge+gauge_offset+GJP.VolNodeSites()*4), tmp_mat1) ; //multiply by U* field rather than U field

      //tmp_mat1.Trans(tmp_mat2) ;
      //tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;
      //tmp_mat2.Conj(tmp_mat2); //these 3 steps give 0

      tmp_mat1.Dagger(tmp_mat2) ;
      tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;
      

      *(mom+gauge_offset+GJP.VolNodeSites()*4) += tmp_mat2 ;

      //set force for U link
      Matrix* mom_ustar = mom+gauge_offset+GJP.VolNodeSites()*4;
      Matrix* mom_u = mom+gauge_offset;
      mom_u->Conj((double*)mom_ustar);


    } } } } // end for x,y,z,t
  } // end for mu
 
//------------------------------------------------------------------
// deallocate smalloc'd space
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_site_v2, site_v2) ;
  sfree(site_v2) ;
 
  VRB.Sfree(cname, fname, str_site_v1, site_v1) ;
  sfree(site_v1) ;
 
  VRB.Sfree(cname, fname, str_v2, v2) ;
  sfree(v2) ;
 
  VRB.Sfree(cname, fname, str_v1, v1) ;
  sfree(v1) ;
}


void GaugeTransformU(Matrix *gtrans, Lattice &lat){
  Matrix recv_buf;
  Matrix tmp;
  //apply the gauge transformation to U
  for(int t=0;t<GJP.TnodeSites();t++){
    for(int z=0;z<GJP.ZnodeSites();z++){
      for(int y=0;y<GJP.YnodeSites();y++){
	for(int x=0;x<GJP.XnodeSites();x++){
	  int pos[4] = {x,y,z,t};
	  int v_x_off = x + GJP.XnodeSites()*(y+GJP.YnodeSites()*(z+GJP.ZnodeSites()*t));
	  Matrix &v_x = *(gtrans + v_x_off);
	  Matrix vdag_x; vdag_x.Dagger(v_x);

	  for(int mu=0;mu<4;mu++){
	    int u_x_off = lat.GsiteOffset(pos) + mu;
	    Matrix &u_x = *(lat.GaugeField() + u_x_off);

	    //get V_x+mu
	    int posp[4] = {x,y,z,t};
	    posp[mu] = (posp[mu]+1)%GJP.NodeSites(mu);
	    Matrix *v_xpmu_ptr = gtrans + posp[0] + GJP.XnodeSites()*(posp[1]+GJP.YnodeSites()*(posp[2]+GJP.ZnodeSites()*posp[3]));

	    if(pos[mu] == GJP.NodeSites(mu)-1 && GJP.Nodes(mu)!=1){
	      //doesnt need to be fast!
	      getPlusData((double *)&recv_buf, (double *)v_xpmu_ptr, 18, mu);
	      v_xpmu_ptr = &recv_buf; 
	    }
	    //dagger/transpose it
	    Matrix vdag_xpmu;
	    if(GJP.Bc(mu) == BND_CND_GPARITY && pos[mu] == GJP.NodeSites(mu)-1 && GJP.NodeCoor(mu) == GJP.Nodes(mu)-1){
	      vdag_xpmu.Trans(*v_xpmu_ptr);
	    }else{
	      vdag_xpmu.Dagger(*v_xpmu_ptr);
	    }

	    //gauge transform link
	    tmp.DotMEqual(v_x,u_x);
	    u_x.DotMEqual(tmp,vdag_xpmu);

	    if(GJP.Gparity()){
	      u_x_off += 4*GJP.VolNodeSites();
	      Matrix &ustar_x = *(lat.GaugeField() + u_x_off);
	      ustar_x.Conj(u_x);
	    }

	  }
	}
      }
    }
  }
}

void GaugeTransformP(Matrix *gtrans,Matrix* mom,Lattice &lat){
  Matrix tmp;
  //apply the gauge transformation to U and also to mom for later comparison
  for(int t=0;t<GJP.TnodeSites();t++){
    for(int z=0;z<GJP.ZnodeSites();z++){
      for(int y=0;y<GJP.YnodeSites();y++){
	for(int x=0;x<GJP.XnodeSites();x++){
	  int pos[4] = {x,y,z,t};
	  int v_x_off = x + GJP.XnodeSites()*(y+GJP.YnodeSites()*(z+GJP.ZnodeSites()*t));
	  Matrix &v_x = *(gtrans + v_x_off);
	  Matrix vdag_x; vdag_x.Dagger(v_x);

	  for(int mu=0;mu<4;mu++){
	    int u_x_off = lat.GsiteOffset(pos) + mu;
	    Matrix &mom_x = *(mom + u_x_off);

	    //gauge transform mom
	    tmp.DotMEqual(v_x,mom_x);
	    mom_x.DotMEqual(tmp,vdag_x);

	    if(GJP.Gparity()){
	      u_x_off += 4*GJP.VolNodeSites();
	      Matrix &momstar_x = *(mom + u_x_off);
	      momstar_x.Conj(mom_x);
	    }
	  }
	}
      }
    }
  }
}

void GaugeTransformF(Matrix *gtrans,Vector* ferm,Lattice &lat){
  Vector tmp;
  //apply the gauge transformation to the odd checkerboarded fermion
  for(int s=0;s<GJP.SnodeSites();s++){
    for(int t=0;t<GJP.TnodeSites();t++){
      for(int z=0;z<GJP.ZnodeSites();z++){
	for(int y=0;y<GJP.YnodeSites();y++){
	  for(int x=0;x<GJP.XnodeSites();x++){
	    if( (x+y+z+t+s)%2 == 0) continue; //skip even sites

	    int v_off = x + GJP.XnodeSites()*(y+GJP.YnodeSites()*(z+GJP.ZnodeSites()*t));
	    Matrix &v_x = *(gtrans + v_off);

	    int pos[5] = {x,y,z,t,s};
	    int f_off = lat.FsiteOffsetChkb(pos) * lat.SpinComponents();

	    for(int spn=0;spn<lat.SpinComponents();spn++){
	      Vector &f_x = *(ferm+f_off+spn);
	      tmp.DotXEqual(v_x,f_x);
	      f_x = tmp;
	    }

	    if(GJP.Gparity()){
	      Matrix vstar_x; vstar_x.Conj(v_x);
	      f_off += GJP.VolNodeSites()/2 * lat.SpinComponents(); //one half-checkerboard offset for second flav
	      for(int spn=0;spn<lat.SpinComponents();spn++){
		Vector &f2_x = *(ferm+f_off+spn);
		tmp.DotXEqual(vstar_x,f2_x);
		f2_x = tmp;
	      }
	    }
	  }
	}
      }
    }
  }
}





















// void setup_double_latt(Lattice &double_latt, Matrix* orig_gfield, bool gparity_X, bool gparity_Y){
//   //orig latt ( U_0 U_1 ) ( U_2 U_3 ) ( U_4 U_5 ) ( U_6 U_7 )
//   //double tatt ( U_0 U_1 U_2 U_3 ) ( U_4 U_5 U_6 U_7 ) ( U_0* U_1* U_2* U_3* ) ( U_4* U_5* U_6* U_7* )
  
//   if(gparity_Y){ printf("setup_double_latt for quad latt not yet implemented\n"); exit(-1); }

//   //Matrix *orig_gfield = orig_latt.GaugeField();
//   Matrix *dbl_gfield = double_latt.GaugeField();

//   if(gparity_X && !gparity_Y){
//     if(GJP.Xnodes()>1){
//       if(!UniqueID()){ printf("Setting up doubled lattice. sizeof(Float) %d sizeof(IFloat) %d\n",sizeof(Float), sizeof(IFloat)); fflush(stdout); }

//       SingleToDoubleLattice lattdoubler(orig_gfield,double_latt);
//       lattdoubler.Run(gparity_X,gparity_Y);

//       if(!UniqueID()){ printf("Finished setting up doubled lattice\n"); fflush(stdout); }
//     }else{
//       //only one node in X-direction
//       //copy data from orig_latt stored on this node
//       int pos[4];
//       for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]++){
// 	for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]++){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]++){
// 	    for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]++){
// 	      for(int mu=0;mu<4;mu++){
// 		int dbl_site_idx = 4*(pos[0]+GJP.XnodeSites() *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		if(pos[0] < GJP.XnodeSites()/2){
// 		  //site is stored on-node in orig_gfield
// 		  int orig_site_idx = 4*(pos[0]+GJP.XnodeSites()/2 *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		  dbl_gfield[dbl_site_idx] = orig_gfield[orig_site_idx];
// 		}else{
// 		  //site is stored on-node in orig_gfield but needs complex conjugating
// 		  int orig_site_idx = 4*(pos[0]-GJP.XnodeSites()/2 +GJP.XnodeSites()/2 *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		  dbl_gfield[dbl_site_idx].Conj(orig_gfield[orig_site_idx]);
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }
// #if 0
//       printf("Converted single to double lattice. Scan across X-direction of original:\n");
//       for(int mu=0;mu<4;mu++){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]++){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]++){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]++){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites()/2;pos[0]++){
// 		int orig_site_idx = 4*(pos[0]+GJP.XnodeSites()/2 *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		printf("U(%d,%d,%d,%d)_%d = (%f,%f) ",pos[0],pos[1],pos[2],pos[3],mu,*((IFloat*)(orig_gfield+orig_site_idx)),*((IFloat*)(orig_gfield+orig_site_idx)+1));
// 	      }
// 	      printf("\n");
// 	    }
// 	  }
// 	}
//       }
//       printf("\nDoubled lattice:\n");
//       for(int mu=0;mu<4;mu++){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]++){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]++){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]++){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]++){
// 		int dbl_site_idx = 4*(pos[0]+GJP.XnodeSites() *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		printf("U(%d,%d,%d,%d)_%d = (%f,%f) ",pos[0],pos[1],pos[2],pos[3],mu,*((IFloat*)(dbl_gfield+dbl_site_idx)),*((IFloat*)(dbl_gfield+dbl_site_idx)+1));
// 	      }
// 	      printf("\n");
// 	    }
// 	  }
// 	}
// 	printf("\n");
//       }
// #endif
//     }


//   }//if(gparity_X && ! gparity_Y)....
// }



// void setup_double_rng(bool gparity_X, bool gparity_Y){
//   //orig 4D rng 2 stacked 4D volumes
//   //orig ([R_0 R_1][R'_0 R'_1])([R_2 R_3][R'_2 R'_3])([R_4 R_5][R'_4 R'_5])([R_6 R_7][R'_6 R'_7])
//   //double (R_0 R_1 R_2 R_3)(R_4 R_5 R_6 R_7)(R'_0 R'_1 R'_2 R'_3)(R'_4 R'_5 R'_6 R'_7)
  
//   //orig 5D rng 2 stacked 4D volumes per ls/2 slice (ls/2 as only one RNG per 2^4 block)

//   if(gparity_Y){ printf("setup_double_latt for quad latt not yet implemented\n"); exit(-1); }

//   if(gparity_X && !gparity_Y){
//     if(!UniqueID()) printf("Setting up RNG from original stacked version\n");


//   //orig 4D rng 2 stacked 4D volumes
//   //orig   ([R_0 R_1][R'_0 R'_1]) ([R_2 R_3][R'_2 R'_3]) ([R_4 R_5][R'_4 R'_5]) ([R_6 R_7][R'_6 R'_7]) ([R_8 R_9][R'_8 R'_9]) ([R_10 R_11][R'_10 R'_11]) ([R_12 R_13][R'_12 R'_13]) ([R_14 R_15][R'_14 R'_15])
//   //double (R_0 R_1 R_2 R_3)      (R_4 R_5 R_6 R_7)      (R_8 R_9 R_10 R_11)    (R_12 R_13 R_13 R_15)  (R'_0 R'_1 R'_2 R'_3)  (R'_4 R'_5 R'_6 R'_7)      (R'_8 R'_9 R'_10 R'_11)    (R'_12 R'_13 R'_14 R'_15)

//     if(GJP.Xnodes()>1){
//       SingleToDouble4dRNG fourDsetup;
//       SingleToDouble5dRNG fiveDsetup;

//       LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form
      
//       if(!UniqueID()){ printf("Setting up 4D RNG\n"); fflush(stdout); }
//       fourDsetup.Run(gparity_X,gparity_Y);      
//       if(!UniqueID()){ printf("Setting up 5D RNG\n"); fflush(stdout); }
//       fiveDsetup.Run(gparity_X,gparity_Y);    
//     }else{    
//       int n_rgen_4d = GJP.VolNodeSites()/16; //applies both to original and doubled latt
//       int n_rgen = n_rgen_4d;
//       if (GJP.SnodeSites()>=2)
// 	n_rgen = GJP.VolNodeSites()*GJP.SnodeSites() / 32;

//       int stk_index_4d_off = n_rgen_4d/2; //offset for R' on 4D orig latt
//       int blocks_per_s_layer = n_rgen /( GJP.SnodeSites() / 2 ); //also same for original and doubled latt
//       int stk_index_5d_off = blocks_per_s_layer/2; //offset for R' on 5D orig latt

//       //copy the originals
//       UGrandomGenerator *ugran_4d_orig = new UGrandomGenerator[n_rgen_4d];
//       for(int i=0;i<n_rgen_4d;i++) ugran_4d_orig[i] = LRG.UGrandGen4D(i);

//       UGrandomGenerator *ugran_orig = new UGrandomGenerator[n_rgen];
//       for(int i=0;i<n_rgen;i++) ugran_orig[i] = LRG.UGrandGen(i);

//       LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form

  
//       int pos[5];

//       for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
// 		//do the 4D RNG
// 		if(pos[4]==0){
// 		  if(pos[0]>=GJP.XnodeSites()/2){
// 		    int orig_idx = (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + stk_index_4d_off;
// 		    int new_idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		    LRG.UGrandGen4D(new_idx) = ugran_4d_orig[orig_idx];
// 		  }else{
// 		    int orig_idx = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		    int new_idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		    LRG.UGrandGen4D(new_idx) = ugran_4d_orig[orig_idx];
// 		  }
// 		}
// 		//do the 5D RNG
// 		if(pos[0]>=GJP.XnodeSites()/2){
// 		  int orig_idx = pos[4]/2*blocks_per_s_layer + (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + stk_index_5d_off;
// 		  int new_idx = pos[4]/2*blocks_per_s_layer + pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		  LRG.UGrandGen(new_idx) = ugran_orig[orig_idx];
// 		}else{
// 		  int orig_idx = pos[4]/2*blocks_per_s_layer + pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		  int new_idx = pos[4]/2*blocks_per_s_layer + pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
// 		  LRG.UGrandGen(new_idx) = ugran_orig[orig_idx];
// 		}

// 	      }
// 	    }
// 	  }
// 	}
//       }
//       delete[] ugran_4d_orig;
//       delete[] ugran_orig;
//     }//single node

//   }//gpx and gpy


// }









