//CK: In this test we check that the gauge force term is evaluated correctly

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
#include <util/pt.h>
#ifdef HAVE_BFM
#include <chroma.h>
#endif

#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

void setup_double_mom(Matrix *double_mom, Matrix* orig_mom, bool gparity_X, bool gparity_Y);
void setup_double_latt(Lattice &double_latt, Matrix* orig_gfield, bool gparity_X, bool gparity_Y);
void setup_double_rng(bool gparity_X, bool gparity_Y);

void GaugeTransformU(Matrix *gtrans, Lattice &lat);
void GaugeTransformP(Matrix *gtrans, Matrix* mom, Lattice &lat);
void GaugeTransformF(Matrix *gtrans,Vector* ferm,Lattice &lat);
ForceArg EvolveMomGforceTest(Lattice &lattice,Matrix *mom, double dt);

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

  printf("Lattice size is %d %d %d %d %d\n",size[0],size[1],size[2],size[3],size[4]);

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

#if TARGET == BGQ
  omp_set_num_threads(64);
#endif

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }
  
  GimprRectFdwf* lattice = new GimprRectFdwf;
					       
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
#define GAUGE_FORCE_CHECK

#ifdef GAUGE_FORCE_CHECK
  Matrix* mom;
  {
    Lattice &lat = *lattice;
    int g_size = 2*GJP.VolNodeSites() * lat.GsiteSize();
    mom = (Matrix*)smalloc(g_size*sizeof(double),"mom","","");
    for(int i=0;i<2*4*GJP.VolNodeSites();i++) mom[i].ZeroMatrix();

    lat.EvolveMomGforce(mom,1);
    //EvolveMomGforceTest(lat,mom,1);

    //above only does mom for U field (by design) but for comparison with double latt, we want both
    double* momfl = (double*)mom;
    for(int i=g_size/2;i<g_size;i++){
      if(i%2==0) momfl[i] = momfl[i-g_size/2];
      else momfl[i] = -momfl[i-g_size/2];
    }
  }
#endif

  //#define PT_CHECK

#ifdef PT_CHECK
  if(GJP.Xnodes()>1){
    Lattice &lat = *lattice;
    ParTransGauge pt(lat);

    int N = 4;
    int vol = GJP.VolNodeSites();  
    Matrix *tmp1[N];
    for(int i = 0;i<N;i++){
      tmp1[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
    }

    int dirs_f[] = {1,3,5,7}; //mu = 0,1,2,3  psi(x) = U_mu^dag(x-mu)psi(x-mu) (forward transport), both U^dag(x-mu) and psi(x-mu) are communicated.
    int dirs_b[] = {0,2,4,6}; //mu = 0,1,2,3  psi(x) = U_mu(x)psi(x+mu) (backward transport), psi(x+mu) is communicated.

    Matrix *Unit = (Matrix *) fmalloc(vol*sizeof(Matrix));
    for(int i = 0; i < vol;i++)
      Unit[i] = 1.;

    Matrix *Units[4];
    for(int i = 0; i < N;i++)
      Units[i] = Unit;
    
    pt.run(N,tmp1,Units,dirs_m);
  }


#endif


  if(gauge_fix) lattice->FixGaugeFree();

  if(UniqueID()==0) printf("Starting double lattice force\n");
 
  if(!UniqueID()) printf("Backing up old gauge field\n"); fflush(stdout);
  int array_size = 2*lattice->GsiteSize() * GJP.VolNodeSites() * sizeof(double);
  Matrix *orig_lattice = (Matrix *) smalloc(array_size,"","","");
  for(int i=0;i<2*4*GJP.VolNodeSites();i++) orig_lattice[i] = lattice->GaugeField()[i];
  //memcpy((void*)orig_lattice, (void*)(lattice->GaugeField()), array_size);  (causes SEGV for some reason...)

  lattice->FreeGauge(); //free memory and reset
  delete lattice; //lattice objects are singleton (scope_lock)

  //setup 1f model. Upon calling GJP.Initialize the lattice size will be doubled in the appropriate directions
  //and the boundary condition set to APRD
  if(gparity_X) do_arg.gparity_1f_X = 1;
  if(gparity_Y) do_arg.gparity_1f_Y = 1;

  if(!UniqueID())printf("Reinitializing GJP\n"); fflush(stdout);
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

  GimprRectFdwf* doubled_lattice = new GimprRectFdwf;
  if(!UniqueID())printf("Setting up doubled latt\n"); fflush(stdout);
  setup_double_latt(*doubled_lattice,orig_lattice,gparity_X,gparity_Y);
  if(!UniqueID())printf("Setting up doubled RNG\n"); fflush(stdout);
  setup_double_rng(gparity_X,gparity_Y);
  if(!UniqueID())printf("Finished doubling things\n"); fflush(stdout); 

  if(gauge_fix){
    doubled_lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    doubled_lattice->FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }

#ifdef GAUGE_FORCE_CHECK
  {
    Lattice &lat = *doubled_lattice;
    int g_size = GJP.VolNodeSites() * lat.GsiteSize();

    Matrix* mom_doubled = (Matrix*)smalloc(g_size*sizeof(double),"mom_doubled","","");
    setup_double_mom(mom_doubled,mom,gparity_X,gparity_Y); //put 2f mom in same format as 1f doubled-latt mom


    Matrix* mom2 = (Matrix*)smalloc(g_size*sizeof(double),"mom2","","");
    for(int i=0;i<4*GJP.VolNodeSites();i++) mom2[i].ZeroMatrix();

    lat.EvolveMomGforce(mom2,1);

    //mom2 and mom_doubled should be the same
    {
      bool err(false);
    
      for(int t=0;t<GJP.TnodeSites();t++){
	for(int z=0;z<GJP.ZnodeSites();z++){
	  for(int y=0;y<GJP.YnodeSites();y++){
	    for(int x=0;x<GJP.XnodeSites();x++){
	      int pos[4] = {x,y,z,t};
	      for(int mu=0;mu<4;mu++){
		int off = lat.GsiteOffset(pos) + mu;
		double* m = (double*)(mom_doubled+off);
		double* mc = (double*)(mom2+off);
		if(fabs(*m - *mc) > 1e-06 || fabs(*(m+1) - *(mc+1)) > 1e-06 ){
		  printf("Error: 2f mom: 1f mom (%d %d %d %d), %d: (%f %f), (%f %f)\n",x,y,z,t,mu,*m,*(m+1),*mc, *(mc+1));
		  err=true;
		}
	      }
	    }
	  }
	}
      }
    
      if(err) exit(-1);
      if(!UniqueID()) printf("Passed gauge force test\n");
    }
    sfree(mom);
    sfree(mom_doubled);
    sfree(mom2);
  }
#endif
  if(gauge_fix) doubled_lattice->FixGaugeFree();
  doubled_lattice->FreeGauge(); //free memory and reset
  delete doubled_lattice; //lattice objects are singleton (scope_lock)

#ifdef HAVE_BFM
  Chroma::finalize();
#endif

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}



void setup_double_mom(Matrix *double_mom, Matrix* orig_mom, bool gparity_X, bool gparity_Y){
  //orig latt ( U_0 U_1 ) ( U_2 U_3 ) ( U_4 U_5 ) ( U_6 U_7 )
  //double tatt ( U_0 U_1 U_2 U_3 ) ( U_4 U_5 U_6 U_7 ) ( U_0* U_1* U_2* U_3* ) ( U_4* U_5* U_6* U_7* )
  
  if(!UniqueID()){ printf("Setting up doubled mom.\n");  fflush(stdout); }
  SingleToDoubleMatrixField momdoubler(gparity_X,gparity_Y,4,orig_mom, double_mom);
  momdoubler.Run();

  if(!UniqueID()){ printf("Finished setting up doubled mom\n"); fflush(stdout); }
}





void setup_double_latt(Lattice &double_latt, Matrix* orig_gfield, bool gparity_X, bool gparity_Y){
  //orig latt ( U_0 U_1 ) ( U_2 U_3 ) ( U_4 U_5 ) ( U_6 U_7 )
  //double tatt ( U_0 U_1 U_2 U_3 ) ( U_4 U_5 U_6 U_7 ) ( U_0* U_1* U_2* U_3* ) ( U_4* U_5* U_6* U_7* )
  if(!UniqueID()){ printf("Setting up doubled lattice.\n"); }

  SingleToDoubleLattice lattdoubler(gparity_X,gparity_Y,orig_gfield,double_latt);
  lattdoubler.Run();
  
  if(!UniqueID()){ printf("Finished setting up doubled lattice\n"); fflush(stdout); }
}



void setup_double_rng(bool gparity_X, bool gparity_Y){
  //orig 4D rng 2 stacked 4D volumes
  //orig ([R_0 R_1][R'_0 R'_1])([R_2 R_3][R'_2 R'_3])([R_4 R_5][R'_4 R'_5])([R_6 R_7][R'_6 R'_7])
  //double (R_0 R_1 R_2 R_3)(R_4 R_5 R_6 R_7)(R'_0 R'_1 R'_2 R'_3)(R'_4 R'_5 R'_6 R'_7)
  
  //orig 5D rng 2 stacked 4D volumes per ls/2 slice (ls/2 as only one RNG per 2^4 block)
  if(!UniqueID()) printf("Setting up RNG from original stacked version\n");


  //orig 4D rng 2 stacked 4D volumes
  //orig   ([R_0 R_1][R'_0 R'_1]) ([R_2 R_3][R'_2 R'_3]) ([R_4 R_5][R'_4 R'_5]) ([R_6 R_7][R'_6 R'_7]) ([R_8 R_9][R'_8 R'_9]) ([R_10 R_11][R'_10 R'_11]) ([R_12 R_13][R'_12 R'_13]) ([R_14 R_15][R'_14 R'_15])
  //double (R_0 R_1 R_2 R_3)      (R_4 R_5 R_6 R_7)      (R_8 R_9 R_10 R_11)    (R_12 R_13 R_13 R_15)  (R'_0 R'_1 R'_2 R'_3)  (R'_4 R'_5 R'_6 R'_7)      (R'_8 R'_9 R'_10 R'_11)    (R'_12 R'_13 R'_14 R'_15)

  SingleToDouble4dRNG fourDsetup(gparity_X,gparity_Y);
  SingleToDouble5dRNG fiveDsetup(gparity_X,gparity_Y);

  LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form
      
  if(!UniqueID()){ printf("Setting up 4D RNG\n"); fflush(stdout); }
  fourDsetup.Run();      
  if(!UniqueID()){ printf("Setting up 5D RNG\n"); fflush(stdout); }
  fiveDsetup.Run();
  if(!UniqueID()){ printf("Finished setting up RNGs\n"); fflush(stdout); }
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


ForceArg EvolveMomGforceTest(Lattice &lattice,Matrix *mom, double dt)
{
  double L1=0.0;
  double L2=0.0;
  double Linf=0.0;

  static Matrix mt0;
  static Matrix *mp0 = &mt0;

  double plaq_coeff = - GJP.Beta() * ( 1.0 - 8.0 * GJP.C1() ) / 3.0 ;
  double rect_coeff = - GJP.Beta() * (             GJP.C1() ) / 3.0 ;
  

  static int vol = GJP.VolNodeSites();  //Local lattice volume
  const int N = 4;                      //Num of dimensions
  double tmp_plaq = plaq_coeff;          //1-8*c_1
  double tmp_rect = rect_coeff;          //c_1

  //Pointer to block of unit matrices for each lattice site
  Matrix *Unit = (Matrix *) fmalloc(vol*sizeof(Matrix));

  //Set all these matrices to the identity
  for(int i = 0; i < vol;i++)
    Unit[i] = 1.;

  //Temporary matrices for use in calculation of the staples
  Matrix *tmp1[N];
  Matrix *tmp2[N];

  for(int i = 0;i<N;i++)
    {
      tmp1[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
      tmp2[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
      //Set all bytes in tmp2 to zero
      bzero((char*)tmp2[i],vol*sizeof(Matrix));
    }

  //Holds the sum of staples associated with each lattice site
  Matrix *result[4];
  for(int i = 0;i<4;i++){
    result[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
  }

  //Array of four pointers to a block of unit matrices
  Matrix *Units[4];
  for(int i = 0; i < N;i++)
    Units[i] = Unit;

  int mu,nu;
  {
    int dirs_p[] = {0,2,4,6,0,2,4};   //Positive directions
    int dirs_m[] = {1,3,5,7,1,3,5};   //Negative directions

    //in pt.C : int mu = dir[d]; bool backwards = mu%2 ? true : false; mu /= 2;
    //hence dirs_p corresponds destination offsets of {+x, +y, +z, +t, +x, +y, +z}
    //      dirs_m                                    {-x, -y, -z, -t, -x, -y, -z}

    //Instantiate parallel transporter
    ParTransGauge pt(lattice);

    //Take a vector field V(x).  Suppose we wish to parallel transport
    //this vector field from the point x to the point x-mu.  That is,
    //we wish to find a new field V'(x-mu) in terms of V(x) such that
    //the new field V'(x-mu) changes appropriately under gauge transformation.
    //
    //The parallel transport that satisfies this property is
    //V'(x-mu) = U_mu(x-mu)V(x).  This combination transforms like
    //a vector field.
    //
    //When using the parallel transport class, one must specify an intial
    //field, a direction, and a final field.  However, the direction
    //that one must specify is the direction of the link matrix, not the
    //direction of the parallel transport.  Thus, to calculate V' as above,
    //use:
    //
    //pt.run(1,V',V,dir_Plus_mu)
    //
    //The new vector field will be indexed according to its new position.

    //CK: Above explanation is very confusing for me!
    //heres what the code actually does:
    //if direction is 1,3,5,7 V'(x) = U_mu(x-mu)V(x-mu)  where mu = x,y,z,t  
    //if direction is 0,2,4,6 V'(x) = U_mu^dag(x)V(x+mu) where mu = x,y,z,t

      for(nu = 1;nu<4;nu++){

	//First calculate the staple in the positive nu direction
	//position is 'p'
        pt.run(N,tmp1,Units,dirs_m+nu); //tmp1(p) = U_nu(p-\nu)
  	pt.run(N,result,tmp1,dirs_m); //result(p) = U_x(p-x) tmp1(p-x) = U_x(p-x)U_nu(p-\nu-x)
	pt.run(N,tmp1,result,dirs_p+nu); //tmp1(p) = U_nu^dag(p) result(p+\nu) = U_nu^dag(p)U_x(p-x+nu)U_nu(p-x)

	//tmp2 contains the sum of the staples for a given link
	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_plaq,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_plaq,tmp1[i],tmp2[i],vol);

	//Calculating one rectangular staple
	pt.run(N,tmp1,Units,dirs_p); //tmp1(p) = U_x^dag(p)
	pt.run(N,result,tmp1,dirs_m+nu); //result(p) = U_mu(p-nu)tmp(p-nu) = U_nu(p-nu)U_x^dag(p-nu)
	pt.run(N,tmp1,result,dirs_m); //tmp1(p) = U_x(p-x)result(p-x) = U_x(p-x)U_nu(p-nu-x)U_x^dag(p-nu-x)
	pt.run(N,result,tmp1,dirs_m); //result(p) = U_x(p-x)tmp1(p-x) = U_x(p-x)U_x(p-2x)U_nu(p-nu-2x)U_x^dag(p-nu-2x)
	pt.run(N,tmp1,result,dirs_p+nu); //tmp1(p) = U_nu^dag(p)result(p+nu) = U_nu^dag(p)U_x(p-x+nu)U_x(p-2x+nu)U_nu(p-2x)U_x^dag(p-2x)

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_m+nu);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_p);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_m+nu);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_p+nu);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Calculating the staple in the negative nu direction
	pt.run(N,tmp1,Units,dirs_p+nu);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m+nu);

	//Add this result into tmp2
	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_plaq,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_plaq,tmp1[i],tmp2[i],vol);

	//Calculating one rectangular staple
	pt.run(N,tmp1,Units,dirs_p);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m+nu);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_p+nu);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_p);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_p+nu);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_m+nu);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);
      }
      //Multiply on the left by our original link matrix to get force term
      pt.run(N,result,tmp2,dirs_p); //result(p) = U_x^dag(p)tmp2(p+x)
  }

#if 1
#pragma omp parallel for default(shared) private(mu) reduction(+:L1,L2)
  for(int index=0;index<4*vol;index++){
    Matrix mp1;
    int i = index%vol;
    mu = index/vol;
    Matrix *mtmp = (result[mu]+i);
    mp1.Dagger((double *)mtmp);
    mtmp->TrLessAntiHermMatrix(mp1);
    double *ihp = (double *)(mom+i*4+mu);  //The gauge momentum
    double *dotp = (double *)mp0;
    double *dotp2 = (double *) (result[mu]+(i));
    fTimesV1PlusV2Single(ihp, dt, dotp2, ihp, 18);  //Update the gauge momentum
    double norm = ((Matrix*)dotp2)->norm();
    double tmp = sqrt(norm);
    L1 += tmp;
    L2 += norm;
  }
#else

      Matrix mp1;
      for(mu = 0; mu<4;mu++)
	{
	  Matrix *mtmp = result[mu];
	  // Takes TrLessAntiHerm part to get the force
	  for(int i = 0; i<vol;i++)
	    {
	      #if 1
	      mp1.Dagger((double *)mtmp);
	      mtmp->TrLessAntiHermMatrix(mp1);
	      #else
	      mtmp->TrLessAntiHermMatrix();
	      #endif
	      mtmp++;
	    }
	}
      ForceFlops += vol*60;

  int x[4];
  
  for(x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0]) {
    for(x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1]) {
      for(x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2]) {
	for(x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {
	  
	  //Calculates the array index offset for the gauge links
	  //located at lattice point x
	  int uoff = lattice.GsiteOffset(x);
	  
	  for (int mu = 0; mu < 4; ++mu) {
//	    GforceSite(*mp0, x, mu);   //Old force calculation
	    
	    double *ihp = (double *)(mom+uoff+mu);  //The gauge momentum
	    double *dotp = (double *)mp0;
	    double *dotp2 = (double *) (result[mu]+(uoff/4));
	    fTimesV1PlusV2(ihp, dt, dotp2, ihp, 18);  //Update the gauge momentum
	    double norm = ((Matrix*)dotp2)->norm();
	    double tmp = sqrt(norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp>Linf ? tmp : Linf);
	  }
	}
      }
    }
  }
}
#endif

  //Free some memory
  ffree(Unit);
  for(int i = 0;i<N;i++){
  ffree(tmp1[i]);
  ffree(tmp2[i]);
  }
  for(int i = 0;i<4;i++) 
  ffree(result[i]);

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);

}
