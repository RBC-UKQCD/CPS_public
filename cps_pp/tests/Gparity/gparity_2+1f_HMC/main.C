#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include <string>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include<util/lattice.h>
#include<util/random.h>
#include<util/time_cps.h>

#include<alg/alg_hmc.h>
#include<alg/common_arg.h>
#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>

#include<alg/alg_int.h>
#include<alg/int_arg.h>

#include<alg/no_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_plaq.h>
#include<alg/alg_pbp.h>
#include<alg/pbp_arg.h>
#include<alg/alg_remez.h>
#include<alg/alg_wline.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>
#include<util/command_line.h>
//#include <sys/bgl/bgl_sys_all.h>

#undef USE_SCU_CHECKSUMS
#ifdef USE_SCU_CHECKSUMS
#include <qcdocos/scu_checksum.h>
#endif
#include<util/rcomplex.h>
#include<util/dirac_op.h>

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include <alg/alg_fix_gauge.h>
#include <alg/fix_gauge_arg.h>

#include <util/data_shift.h>

#include <alg/lanc_arg.h>
#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>

#include <util/gparity_singletodouble.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif
//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

HmcArg hmc_arg;
HmcArg hmc_arg_pass;

ActionGaugeArg gauge_arg;
ActionRationalQuotientArg rat_quo_l_arg;
ActionRationalQuotientArg rat_quo_h_arg;
IntABArg ab1_arg;
IntABArg ab2_arg;
IntABArg ab3_arg;
DoArg do_arg;

void checkpoint(int traj);
void setupRatQuoArg(ActionRationalQuotientArg &into, const int &ndet, Float* bsn_masses, Float* frm_masses, int *pwr_num, int *pwr_den);
void setup_double_latt(Lattice &double_latt, Matrix* orig_gfield, bool gparity_X, bool gparity_Y);
void setup_double_rng(bool gparity_X, bool gparity_Y);
void setup_double_matrixfield(Matrix* double_mat, Matrix* orig_mat, int nmat_per_site, bool gparity_X, bool gparity_Y);

char * strconcat(char *a,char*b){
  char *tmp = new char[strlen(a)+strlen(b)+4];
  sprintf(tmp,"%s%s",a,b);
  return tmp;
}

int main(int argc, char *argv[])
{ 
  Start(&argc,&argv);
  CommandLine::is(argc,argv);

  char *cname=argv[0];
  char *fname="main()";
  Float dtime;

  bool gparity_X(false);
  bool gparity_Y(false);

  int arg0 = CommandLine::arg_as_int(0);
  printf("Arg0 is %d\n",arg0);
  if(arg0==1){
    gparity_X=true;
    printf("Doing G-parity HMC test in X direction\n");
  }else if(arg0==2){
    printf("Doing G-parity HMC test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
  }else{
    ERR.General("main","","Invalid test index");
  }

  bool dbl_latt_storemode(false);
  bool load_config(false);
  char *load_config_file;
  int size[] = {4,4,4,4,4};
  bool save_lrg(false);
  char *save_lrg_file;
  bool load_lrg(false);
  char *load_lrg_file;
  int nsteps = 4;

  bool readwrite_remez_poles(false);
  char* remez_file_stub;

  int i=2;
  while(i<argc){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-load_config",15) == 0){
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
    }else if( strncmp(cmd,"-nsteps",10) == 0){
      if(i==argc-1){ printf("-nsteps requires an argument\n"); exit(-1); }
      nsteps = atoi(argv[i+1]);
      i+=2; 
    }else if( strncmp(cmd,"-save_lrg",15) == 0){
      if(i==argc-1){ printf("-save_lrg requires an argument\n"); exit(-1); }
      save_lrg=true;
      save_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-load_lrg",15) == 0){
      if(i==argc-1){ printf("-load_lrg requires an argument\n"); exit(-1); }
      load_lrg=true;
      load_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-readwrite_remez_poles",15) == 0){
      if(i==argc-1){ printf("-readwrite_remez_poles requires an argument\n"); exit(-1); }
      readwrite_remez_poles=true;
      remez_file_stub = argv[i+1];
      i+=2;
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }
  SerialIO::dbl_latt_storemode = dbl_latt_storemode;

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
  if(load_config){
    do_arg.start_conf_kind = START_CONF_FILE;
    do_arg.start_conf_filename = load_config_file;
  }else{
    do_arg.start_conf_kind = START_CONF_ORD;
  }
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

  if(gparity_X) do_arg.x_bc = BND_CND_GPARITY;
  if(gparity_Y) do_arg.y_bc = BND_CND_GPARITY;

  CommonArg common_arg_hmc;
 
  // do_arg.verbose_level=VERBOSE_RESULT_LEVEL;
  GJP.Initialize(do_arg);
  // VRB.Level(VERBOSE_RESULT_LEVEL);
  LRG.Initialize();
  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }

  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }

  // Outer config loop

  gauge_arg.gluon = G_CLASS_IMPR_RECT;
  gauge_arg.action_arg.force_measure = FORCE_MEASURE_YES;
  gauge_arg.action_arg.force_label = "Gauge";
  
  {
    Float bsn_mass = 1.0;
    Float frm_mass = 0.4;
    int pwr_num = 1;
    int pwr_den = 2;
    setupRatQuoArg(rat_quo_l_arg,1, &bsn_mass, &frm_mass, &pwr_num, &pwr_den);
  }    

  {
    Float bsn_mass = 1.0;
    Float frm_mass = 0.5;
    int pwr_num = 1;
    int pwr_den = 4;
    setupRatQuoArg(rat_quo_h_arg,1, &bsn_mass, &frm_mass, &pwr_num, &pwr_den);
  }  

  hmc_arg.steps_per_traj = 1;
  hmc_arg.step_size =   1.2500000000000000e-01;
  hmc_arg.metropolis = METROPOLIS_NO;//YES;
  hmc_arg.reunitarize = REUNITARIZE_YES;
  hmc_arg.reverse = REVERSE_NO;
  hmc_arg.reproduce = REPRODUCE_NO;
  hmc_arg.reproduce_attempt_limit = 1;
  hmc_arg.wfm_md_sloppy = 1;

  hmc_arg_pass = hmc_arg;

  //!< Construct numerical integrators
  ab1_arg.type = INT_OMELYAN;
  ab1_arg.A_steps = 1;
  ab1_arg.B_steps = 1;
  ab1_arg.level = EMBEDDED_INTEGRATOR;
  ab1_arg.lambda =   2.2000000000000000e-01;
    
  ab2_arg.type = INT_OMELYAN;
  ab2_arg.A_steps = 1;
  ab2_arg.B_steps = 1;
  ab2_arg.level = EMBEDDED_INTEGRATOR;
  //ab2_arg.level = TOP_LEVEL_INTEGRATOR;
  ab2_arg.lambda =   2.2000000000000000e-01;

  ab3_arg.type = INT_OMELYAN;
  ab3_arg.A_steps = 1;
  ab3_arg.B_steps = 1;
  ab3_arg.level = TOP_LEVEL_INTEGRATOR;
  ab3_arg.lambda =   2.2000000000000000e-01;

  //!< Create fictitous Hamiltonian (mom + action)
  AlgMomentum mom;
  AlgActionGauge gauge(mom, gauge_arg);

  if(readwrite_remez_poles) rat_quo_h_arg.rat_poles_file = strconcat(remez_file_stub,"_2f_h.remez");   
  if(readwrite_remez_poles) rat_quo_l_arg.rat_poles_file = strconcat(remez_file_stub,"_2f_l.remez");   

  AlgActionRationalQuotient rat_quo_h(mom, rat_quo_h_arg);
  AlgActionRationalQuotient rat_quo_l(mom, rat_quo_l_arg);

  AlgIntOmelyan ab1(mom, gauge, ab1_arg);
  AlgIntOmelyan ab2(ab1, rat_quo_h, ab2_arg);
  AlgIntOmelyan ab3(ab2, rat_quo_l, ab3_arg);

  //switch to 1f mode to setup integrators
  if(!UniqueID()){ printf("Switching back to 1f mode\n"); fflush(stdout); }

  {
    Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    //deallocate lattice
    lattice.FreeGauge();
    LatticeFactory::Destroy();
  }

  if(gparity_X) do_arg.gparity_1f_X = 1;
  if(gparity_Y) do_arg.gparity_1f_Y = 1;

  GJP.Initialize(do_arg);

  AlgMomentum mom_1f;
  AlgActionGauge gauge_1f(mom_1f, gauge_arg);

  if(readwrite_remez_poles) rat_quo_h_arg.rat_poles_file = strconcat(remez_file_stub,"_1f_h.remez");   
  if(readwrite_remez_poles) rat_quo_l_arg.rat_poles_file = strconcat(remez_file_stub,"_1f_l.remez");   

  AlgActionRationalQuotient rat_quo_h_1f(mom_1f, rat_quo_h_arg);
  AlgActionRationalQuotient rat_quo_l_1f(mom_1f, rat_quo_l_arg);

  AlgIntOmelyan ab1_1f(mom_1f, gauge_1f, ab1_arg);
  AlgIntOmelyan ab2_1f(ab1_1f, rat_quo_h_1f, ab2_arg);
  AlgIntOmelyan ab3_1f(ab2_1f, rat_quo_l_1f, ab3_arg);

  if(!UniqueID()){ printf("Switching back to 2f mode\n"); fflush(stdout); }

  //switch back to 2f mode
  {
    Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    //deallocate lattice
    lattice.FreeGauge();
    LatticeFactory::Destroy();
  }

  if(gparity_X){ do_arg.gparity_1f_X = 0; do_arg.x_bc = BND_CND_GPARITY; }
  if(gparity_Y){ do_arg.gparity_1f_Y = 0; do_arg.y_bc = BND_CND_GPARITY; }

  GJP.Initialize(do_arg);
  GJP.EnableGparity1f2fComparisonCode();
  //start config loop

  if(!UniqueID()){ printf("Starting trajectory loop\n"); fflush(stdout); }

  //FILE *fp = fopen("links.dat","w");

  for(int conf=0; conf< 1; conf ++ ) {
    for(int traj=0;traj< nsteps;traj++) {
      int array_size;
      //backup initial lattice and RNG
      {
	Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE); 
	array_size = 2*lattice.GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
	LatticeFactory::Destroy();
      }

      Matrix *_2f_lat_prehmc = (Matrix *) pmalloc(array_size);
      {
	Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	memcpy((void*)_2f_lat_prehmc, (void*)lat.GaugeField(), array_size);
	LatticeFactory::Destroy();
      }
      LatRanGen _2f_LRG_prehmc(LRG);

      //run hybrid Monte Carlo for 2f G-parity
      VRB.Result("","main()","Running traj %d without reproduction\n",traj);
      hmc_arg_pass.reproduce = REPRODUCE_NO;

      {
	if(!UniqueID()){ printf("Starting 2f evolution for traj %d\n",traj); fflush(stdout); }
	AlgHmc hmc(ab3, common_arg_hmc, hmc_arg_pass);
	Float time = -dclock();
	hmc.run();
	time += dclock();
	print_flops("AlgHmc 2f G-parity","run()",0,time);
      }

#define GP_2F_SYMMETRY_CHECK
#ifdef GP_2F_SYMMETRY_CHECK
      if(GJP.Gparity()){
	Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	int base_loc[4];
	for(int i=0;i<4;i++) base_loc[i] = GJP.NodeCoor(i)*GJP.NodeSites(i); //absolute pos of start point
    
	//printf("Lattice StrOrd is %s\n",StrOrdType_map[(int)lattice.StrOrd()].name);

	int pos[4];
	
	for(int t=0;t<GJP.NodeSites(3);t++){
	  pos[3] = t;
	  for(int z=0;z<GJP.NodeSites(2);z++){
	    pos[2] = z;
	    for(int y=0;y<GJP.NodeSites(1);y++){
	      pos[1] = y;
	      for(int x=0;x<GJP.NodeSites(0);x++){
		pos[0] = x;
		
		for(int mu=0;mu<4;mu++){
		  Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu));
		  Matrix *Ustar = const_cast<Matrix *>(lattice.GetLink(pos, mu,1));

		  Matrix tmp;
		  tmp.Conj((IFloat*)U);
		    
		  for(int ii=0;ii<18;ii++){
		    if(tmp.elem(ii)!=Ustar->elem(ii) ){
		      printf("Symmetry broken on traj idx %d at location %d %d %d %d, %d (elem %d)\n",traj,x,y,z,t,mu,ii); 
		      printf("%.16f != %.16f\n",tmp.elem(ii),Ustar->elem(ii));
		      printf("Matrix should be:\n");
		      for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
			  const Rcomplex &c = tmp(i,j);
			  printf("[%.3f %.3f] ",c.real(), c.imag());
			}
			printf("\n");
		      }
		      printf("\n");
		      
		      printf("But got:\n");
		      for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
			  const Rcomplex &c = (*Ustar)(i,j);
			  printf("[%.3f %.3f] ",c.real(), c.imag());
			}
			printf("\n");
		      }
		      printf("\n");
			
		      exit(-1);
		    }
		  }
		  
		}
		
	      }
	      
	      
	    }
	  }
	}
	
       	LatticeFactory::Destroy();
      }
#endif

      //store 2f G-parity results
      Matrix *_2f_lat_posthmc = (Matrix *) pmalloc(array_size);
      {
	Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	memcpy((void*)_2f_lat_posthmc, (void*)lat.GaugeField(), array_size);
	LatticeFactory::Destroy();
      }
      LatRanGen _2f_LRG_posthmc(LRG);
      
      //restore LRG back to prior to HMC 
      LRG = _2f_LRG_prehmc;

      //deallocate lattice
      {
	Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	//deallocate lattice
	lattice.FreeGauge();
	LatticeFactory::Destroy();
      }

      //repeat evolution from start with 1f G-parity

      if(gparity_X) do_arg.gparity_1f_X = 1;
      if(gparity_Y) do_arg.gparity_1f_Y = 1;

      GJP.Initialize(do_arg);
      GJP.EnableGparity1f2fComparisonCode();

#ifdef HAVE_BFM
      {
	QDP::multi1d<int> nrow(Nd);  
	for(int i = 0;i<Nd;i++) nrow[i] = GJP.Sites(i);
	QDP::Layout::setLattSize(nrow);
	QDP::Layout::create();
      }
#endif
      
      {
	Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE); //reallocates lattice
	setup_double_latt(lattice,_2f_lat_prehmc,gparity_X,gparity_Y);	
	LatticeFactory::Destroy();
      }
      setup_double_rng(gparity_X,gparity_Y);

      //run hybrid Monte Carlo for 1f G-parity
      {
	if(!UniqueID()){ printf("Starting 1f evolution for traj %d\n",traj); fflush(stdout); }
	AlgHmc hmc(ab3_1f, common_arg_hmc, hmc_arg_pass);
	Float time = -dclock();
	hmc.run();
	time += dclock();
	print_flops("AlgHmc 1f G-parity","run()",0,time);
      }

      //boost 2f result to 1f lattice setup
      {
	Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE); 
	array_size = lattice.GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
	LatticeFactory::Destroy();
      }
      Matrix *_2f_lat_posthmc_1fsetup = (Matrix *) pmalloc(array_size);
      setup_double_matrixfield(_2f_lat_posthmc_1fsetup,_2f_lat_posthmc,4,gparity_X,gparity_Y);

      //compare
      {
	bool err(false);
	Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	for(int t=0;t<GJP.TnodeSites();t++){
	  for(int z=0;z<GJP.ZnodeSites();z++){
	    for(int y=0;y<GJP.YnodeSites();y++){
	      for(int x=0;x<GJP.XnodeSites();x++){
		int pos[4] = {x,y,z,t};
		for(int mu=0;mu<4;mu++){
		  int off = lattice.GsiteOffset(pos) + mu;
		  Float* m = (Float*)(lattice.GaugeField()+off);
		  Float* mc = (Float*)(_2f_lat_posthmc_1fsetup+off);
		  if(fabs(*m - *mc) > 1e-06 || fabs(*(m+1) - *(mc+1)) > 1e-06 ){
		    printf("Error: 1f:2f (%d %d %d %d), %d: (%f %f), (%f %f)\n",x,y,z,t,mu,*m,*(m+1),*mc, *(mc+1));
		    err=true;
		  }
		}
	      }
	    }
	  }
	}
	LatticeFactory::Destroy();

	if(err){
	  printf("Failed evolve test on traj step %d\n",traj);
	  exit(-1);
	}
	printf("Passed evolve test on traj step %d\n",traj);
      }

      //restore 2f setup
      
      //deallocate lattice
      {
	Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	//deallocate lattice
	lattice.FreeGauge();
	LatticeFactory::Destroy();
      }

      //repeat evolution from start with 1f G-parity

      if(gparity_X){ do_arg.gparity_1f_X = 0; do_arg.x_bc = BND_CND_GPARITY; }
      if(gparity_Y){ do_arg.gparity_1f_Y = 0; do_arg.y_bc = BND_CND_GPARITY; }

      GJP.Initialize(do_arg);
      GJP.EnableGparity1f2fComparisonCode();
#ifdef HAVE_BFM
      {
	QDP::multi1d<int> nrow(Nd);  
	for(int i = 0;i<Nd;i++) nrow[i] = GJP.Sites(i);
	QDP::Layout::setLattSize(nrow);
	QDP::Layout::create();
      }
#endif

      {
	Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE); 
	array_size = 2*lattice.GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
	LatticeFactory::Destroy();
      }

      //restore LRG and lattice
      LRG = _2f_LRG_posthmc;
      {
	Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE); //reallocates lattice
	memcpy((void*)lattice.GaugeField(),(void*)_2f_lat_posthmc, array_size);
	LatticeFactory::Destroy();
      }

      pfree(_2f_lat_prehmc);
      pfree(_2f_lat_posthmc);
      pfree(_2f_lat_posthmc_1fsetup);
	      
    }//End of inter-cfg sweep

    //#define WRITE_LATT
#ifdef WRITE_LATT
    char lat_file[256];
    // Save this config to disk
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    
    sprintf(lat_file,"gauge.%d",conf);
    QioArg wt_arg(lat_file,0.001);
    
    wt_arg.ConcurIONumber=32;
    WriteLatticeParallel wl;
    wl.setHeader("ens_idx","ens_label",conf);
    wl.write(lat,wt_arg);
    
    if(!wl.good()) 
      ERR.General(cname,fname,"Failed write lattice %s",lat_file);
    LatticeFactory::Destroy();
#endif

  } 

  End();
  printf("End of program\n");
 return(0);
}

void setupRatQuoArg(ActionRationalQuotientArg &into, const int &ndet, Float* bsn_masses, Float* frm_masses, int *pwr_num, int *pwr_den){
  //bi_arg
  into.bi_arg.fermion = F_CLASS_DWF;
  into.bi_arg.bilinears.bilinears_len = ndet;
  into.bi_arg.bilinears.bilinears_val = new BilinearDescr[ndet];
  for(int i=0;i<ndet;i++){
    into.bi_arg.bilinears.bilinears_val[i].mass = 0.0;
    into.bi_arg.bilinears.bilinears_val[i].max_num_iter = 5000;
  }
  into.bi_arg.action_arg.force_measure = FORCE_MEASURE_YES;
  into.bi_arg.action_arg.force_label = "RationalQuotient";
  
  into.spread = 0.0;
  into.remez_generate= 0;
  into.rat_poles_file = "";

  //bsn_mass
  into.bsn_mass.bsn_mass_len = ndet;
  into.bsn_mass.bsn_mass_val = new Float[ndet];
  for(int i=0;i<ndet;i++){
    into.bsn_mass.bsn_mass_val[i] = bsn_masses[i];
  }
  //frm_mass
  into.frm_mass.frm_mass_len = ndet;
  into.frm_mass.frm_mass_val = new Float[ndet];
  for(int i=0;i<ndet;i++){
    into.frm_mass.frm_mass_val[i] = frm_masses[i];
  }
  //bsn_mass_epsilon
  into.bsn_mass_epsilon.bsn_mass_epsilon_len = ndet;
  into.bsn_mass_epsilon.bsn_mass_epsilon_val = new Float[ndet];
  for(int i=0;i<ndet;i++){
    into.bsn_mass_epsilon.bsn_mass_epsilon_val[i] = 0;
  }
  //frm_mass_epsilon
  into.frm_mass_epsilon.frm_mass_epsilon_len = ndet;
  into.frm_mass_epsilon.frm_mass_epsilon_val = new Float[ndet];
  for(int i=0;i<ndet;i++){
    into.frm_mass_epsilon.frm_mass_epsilon_val[i] = 0;
  }

  //bosons
  into.bosons.bosons_len = ndet;
  into.bosons.bosons_val = new RationalDescr[ndet];
  for(int i=0;i<ndet;i++){
    into.bosons.bosons_val[i].field_type = BOSON;
    into.bosons.bosons_val[i].power_num = pwr_num[i];
    into.bosons.bosons_val[i].power_den = pwr_den[i];
    into.bosons.bosons_val[i].precision = 40;
    into.bosons.bosons_val[i].stop_rsd_fg_mult = 1.0;
    
    ApproxDescr *approx = &into.bosons.bosons_val[i].md_approx;
    approx->approx_type = RATIONAL_APPROX_POWER;
    approx->bounds_type = RATIONAL_BOUNDS_MANUAL;
    approx->lambda_low =   4.0000000000000002e-04;
    approx->lambda_high =   2.5000000000000000e+00;
    approx->stop_rsd.stop_rsd_len = 9;
    approx->stop_rsd.stop_rsd_val = new Float[9];
    approx->stop_rsd.stop_rsd_val[0] =   2.0000000000000002e-05;
    approx->stop_rsd.stop_rsd_val[1] =   1.9999999999999999e-06;
    approx->stop_rsd.stop_rsd_val[2] =   1.9999999999999999e-07;
    approx->stop_rsd.stop_rsd_val[3] =   9.9999999999999995e-08;
    approx->stop_rsd.stop_rsd_val[4] =   9.9999999999999995e-08;
    approx->stop_rsd.stop_rsd_val[5] =   9.9999999999999995e-08;
    approx->stop_rsd.stop_rsd_val[6] =   9.9999999999999995e-08;
    approx->stop_rsd.stop_rsd_val[7] =   9.9999999999999995e-08;
    approx->stop_rsd.stop_rsd_val[8] =   9.9999999999999995e-08;
    
    approx = &into.bosons.bosons_val[i].mc_approx;
    approx->approx_type = RATIONAL_APPROX_POWER;
    approx->bounds_type = RATIONAL_BOUNDS_MANUAL;
    approx->lambda_low =   4.0000000000000002e-04;
    approx->lambda_high =   2.5000000000000000e+00;
    approx->stop_rsd.stop_rsd_len = 15;
    approx->stop_rsd.stop_rsd_val = new Float[15];
    for(int j=0;j<15;j++) approx->stop_rsd.stop_rsd_val[j]= 1.0000000000000000e-10; 

    into.bosons.bosons_val[i].stag_bsn_mass = 0.0;
  }
  //fermions
  into.fermions.fermions_len = ndet;
  into.fermions.fermions_val = new RationalDescr[ndet];
  for(int i=0;i<ndet;i++){
    into.fermions.fermions_val[i].field_type = FERMION;
    into.fermions.fermions_val[i].power_num = pwr_num[i];
    into.fermions.fermions_val[i].power_den = pwr_den[i];
    into.fermions.fermions_val[i].precision = 40;
    into.fermions.fermions_val[i].stop_rsd_fg_mult = 1.0;
    
    ApproxDescr *approx = &into.fermions.fermions_val[i].md_approx;
    approx->approx_type = RATIONAL_APPROX_POWER;
    approx->bounds_type = RATIONAL_BOUNDS_MANUAL;
    approx->lambda_low =   4.0000000000000002e-04;
    approx->lambda_high =   2.5000000000000000e+00;
    approx->stop_rsd.stop_rsd_len = 9;
    approx->stop_rsd.stop_rsd_val = new Float[9];
    approx->stop_rsd.stop_rsd_val[0] =   2.0000000000000002e-05;
    approx->stop_rsd.stop_rsd_val[1] =   1.9999999999999999e-06;
    approx->stop_rsd.stop_rsd_val[2] =   1.9999999999999999e-07;
    approx->stop_rsd.stop_rsd_val[3] =   9.9999999999999995e-08;
    approx->stop_rsd.stop_rsd_val[4] =   9.9999999999999995e-08;
    approx->stop_rsd.stop_rsd_val[5] =   9.9999999999999995e-08;
    approx->stop_rsd.stop_rsd_val[6] =   9.9999999999999995e-08;
    approx->stop_rsd.stop_rsd_val[7] =   9.9999999999999995e-08;
    approx->stop_rsd.stop_rsd_val[8] =   9.9999999999999995e-08;

    approx = &into.fermions.fermions_val[i].mc_approx;
    approx->approx_type = RATIONAL_APPROX_POWER;
    approx->bounds_type = RATIONAL_BOUNDS_MANUAL;
    approx->lambda_low =   4.0000000000000002e-04;
    approx->lambda_high =   2.5000000000000000e+00;
    approx->stop_rsd.stop_rsd_len = 15;
    approx->stop_rsd.stop_rsd_val = new Float[15];
    for(int j=0;j<15;j++) approx->stop_rsd.stop_rsd_val[j]= 1.0000000000000000e-10; 

    into.fermions.fermions_val[i].stag_bsn_mass = 0.0;
  }

  into.eigen.eigen_measure = EIGEN_MEASURE_YES;
  into.eigen.stop_rsd = 0.0005;
  into.eigen.max_num_iter = 10000;
  into.eigen.eig_lo_stem = "eig_low";
  into.eigen.eig_hi_stem = "eig_hi";
}







void checkpoint(int traj)
{
  char *cname="cps::";
  char *fname="checkpoint()";


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
