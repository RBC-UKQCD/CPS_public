#if defined(USE_TBC_INPUT)
#warning "Using TBC specified in do_arg"
#elif defined(USE_TBC_FB)
#warning "Using F=P+A and B=P-A temporal boundary condition combinations"
#elif defined(DOUBLE_TLATT)
#warning "Doubling time direction and using specified in do_arg"
#define USE_TBC_INPUT
#else
#error "Must specify TBC flag USE_TBC_INPUT or USE_TBC_FB"
#endif

#ifdef TESTING
#warning "Compiling for TEST mode"
#endif

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
#include<util/lattice/fbfm.h>
#include <util/time_cps.h>
#include <util/smalloc.h>

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <util/command_line.h>

#include<unistd.h>
#include<config.h>

#include <util/data_shift.h>

#include <alg/propmanager.h>
#include <alg/alg_fix_gauge.h>
#include <alg/gparity_contract_arg.h>
#include <alg/alg_gparitycontract.h>
#include <alg/prop_dft.h>
#include <alg/eigen/Krylov_5d.h>

#ifdef USE_EIGEN
//Stupid header guard conflict between eigen and chroma
#undef EIGEN_IO_H
#endif
#include <chroma.h>

#include <omp.h>
#include <pthread.h>

#include <alg/bfm_arg.h>

#include <alg/a2a/threemomentum.h>
#include <alg/a2a/utils.h>
#include <alg/a2a/fmatrix.h>

#include "propmomcontainer.h"
#include "mesonmomenta.h"
#include "wallsinkprop.h"
#include "meas.h"
#include "gparity.h"
#include "lattice_doubleT.h"
#include "cshift.h"
#include "tests.h"


USING_NAMESPACE_CPS

//#define TESTING

int main(int argc,char *argv[])
{
#ifdef TESTING
  return run_tests(argc,argv);
#endif

  Start(&argc,&argv);
  if(argc<2) ERR.General("","main()","Arguments: Require a directory containing the vml files");

  BfmArg bfm_arg; //note, only the solver and mobius scale are used here
  DoArg do_arg;
  LancArg lanc_arg_l;
  LancArg lanc_arg_h;
  GparityAMAarg2 ama_arg;

  decode_vml_all(do_arg, bfm_arg, lanc_arg_l, lanc_arg_h, ama_arg, argv[argc-1]);
  if(ama_arg.conf_start >= ama_arg.conf_lessthan || ama_arg.conf_incr == 0) ERR.General("","main()","Invalid configuration args");
  if(lanc_arg_l.mass != ama_arg.ml) ERR.General("","main()","Light lanczos mass differs from value in AMA args");
  if(lanc_arg_h.mass != ama_arg.mh) ERR.General("","main()","Heavy lanczos mass differs from value in AMA args");

  bool lanczos_tune_l = false;
  bool lanczos_tune_h = false;
  bool dbl_latt_storemode = false;
  Fbfm::use_mixed_solver = true;
  bool mres_do_flavor_project = true;
  bool bk_do_flavor_project = true;
  bool do_alternative_mom = true;
  bool disable_lanczos = false;
  bool random_prop_solns = false; //don't invert, just make the solutions random spin-color-flavor matrices
  bool skip_gauge_fix = false;
  bool random_exact_tsrc_offset = true; //randomly shift the full set of exact src timeslices
  bool rng_test = false; //just load config, generate some uniform random numbers then move onto next config
  int tshift = 0;
  {
    int i = 1;
    while(i<argc-1){
      if( std::string(argv[i]) == "-lanczos_tune_l" ){
	lanczos_tune_l = true;
	i++;
      }else if( std::string(argv[i]) == "-lanczos_tune_h" ){
	lanczos_tune_h = true;
	i++;
      }else if( std::string(argv[i]) == "-load_dbl_latt" ){
	if(!UniqueID()) printf("Loading double latt\n");
	dbl_latt_storemode = true;
	i++;
      }else if( std::string(argv[i]) == "-disable_random_exact_tsrc_offset" ){
	if(!UniqueID()) printf("Disabling random offset of exact tsrc\n");
	random_exact_tsrc_offset = false;
	i++;
      }else if( std::string(argv[i]) == "-disable_mixed_solver" ){
	if(!UniqueID()) printf("Disabling mixed solver\n");
	Fbfm::use_mixed_solver = false;
	i++;
      }else if( std::string(argv[i]) == "-disable_lanczos" ){
	if(!UniqueID()) printf("Not computing or using low-modes\n");
	disable_lanczos = true;
	i++;
      }else if( std::string(argv[i]) == "-disable_mres_flavor_project" ){ //for comparison with old code
	if(!UniqueID()) printf("Disabling mres flavor project\n");
	mres_do_flavor_project = false;
	i++;
      }else if( std::string(argv[i]) == "-disable_bk_flavor_project" ){ //for comparison with old code
	if(!UniqueID()) printf("Disabling BK flavor project\n");
	bk_do_flavor_project = false;
	i++;
      }else if( std::string(argv[i]) == "-disable_use_alternate_mom" ){ 
	if(!UniqueID()) printf("Disabling use of alternative momentum combinations\n");
	do_alternative_mom = false;
	i++;
      }else if( std::string(argv[i]) == "-random_prop_solns"){
	if(!UniqueID()) printf("Not inverting, just using random propagator solutions\n");
	random_prop_solns = true;
	i++;
      }else if( std::string(argv[i]) == "-skip_gauge_fix"){
	skip_gauge_fix = true;
	if(!UniqueID()){ printf("Skipping gauge fixing\n"); fflush(stdout); }
	i++;
      }else if( std::string(argv[i]) == "-rng_test"){
	rng_test = true;
	if(!UniqueID()){ printf("Doing RNG test\n"); fflush(stdout); }
	i++;
      }else if( std::string(argv[i]) == "-tshift_gauge"){
	std::stringstream ss; ss << argv[i+1]; ss >> tshift;
	if(!UniqueID()){ printf("Shifting gauge field by %d in time direction\n",tshift); fflush(stdout); }
	i+=2;
      }else{
	ERR.General("","main","Unknown argument: %s",argv[i]);
      }
    }
  }

  if(UniqueID()==0) printf("Configuration format is '%s'\n",ama_arg.config_fmt);
  if(!contains_pctd(ama_arg.config_fmt)) ERR.General("","main()","GparityAMAarg config format '%s' does not contain a %%d",ama_arg.config_fmt);

  if(UniqueID()==0) printf("RNG format is '%s'\n",ama_arg.rng_fmt);
  if(!contains_pctd(ama_arg.rng_fmt)) ERR.General("","main()","GparityAMAarg rng format '%s' does not contain a %%d",ama_arg.rng_fmt);
  
  GJP.Initialize(do_arg);
  GJP.StartConfKind(START_CONF_MEM); //we will handle the gauge field read/write thankyou!
  GJP.StartSeedKind(START_SEED_INPUT); //and the LRG load!

  if(do_arg.start_seed_kind != START_SEED_FILE) printf("WARNING: Using fixed input start seed. You might want to use START_SEED_FILE to aid reproducibility!\n");

#if TARGET == BGQ
  LRG.setSerial();
#endif

  LRG.Initialize();
  
  init_fbfm(&argc,&argv,bfm_arg);

  if(UniqueID()==0)
    if(Fbfm::use_mixed_solver) printf("Using Fbfm mixed precision solver\n");
    else printf("Using Fbfm double precision solver\n");

  const int Lt = GJP.Tnodes()*GJP.TnodeSites();

  check_bk_tsources(ama_arg); //Check the time seps for BK before we have to do any work

  GnoneFbfm lattice;
  CommonArg carg("label","filename");
  char load_config_file[1000];
  char load_rng_file[1000];

  //Decide on the meson and quark momenta we wish to compute
  MesonMomenta pion_momenta;
  PionMomenta::setup(pion_momenta,do_alternative_mom);
  pion_momenta.printAllCombs("pion_momenta");
    
  MesonMomenta su2_singlet_momenta;
  LightFlavorSingletMomenta::setup(su2_singlet_momenta);
  su2_singlet_momenta.printAllCombs("su2_singlet_momenta");

  MesonMomenta kaon_momenta;
  KaonMomenta::setup(kaon_momenta);
  kaon_momenta.printAllCombs("kaon_momenta");

  //Determine the quark momenta we will need
  QuarkMomenta light_quark_momenta;
  QuarkMomenta heavy_quark_momenta;
    
  pion_momenta.appendQuarkMomenta(Light, light_quark_momenta); //adds the quark momenta it needs
  su2_singlet_momenta.appendQuarkMomenta(Light, light_quark_momenta);
  kaon_momenta.appendQuarkMomenta(Light, light_quark_momenta); //each momentum is unique
  kaon_momenta.appendQuarkMomenta(Heavy, heavy_quark_momenta);

  if(!UniqueID()){
    printf("Light quark momenta to be computed:\n");
    for(int i=0;i<light_quark_momenta.nMom();i++)
      std::cout << light_quark_momenta.getMom(i).str() << '\n';
    printf("Heavy quark momenta to be computed:\n");
    for(int i=0;i<heavy_quark_momenta.nMom();i++)
      std::cout << heavy_quark_momenta.getMom(i).str() << '\n';
  }

  for(int conf=ama_arg.conf_start; conf < ama_arg.conf_lessthan; conf += ama_arg.conf_incr){
    Float conf_start_time = dclock();

    //Read/generate the gauge configuration 
    if(do_arg.start_conf_kind == START_CONF_ORD || rng_test){
      if(!UniqueID()) printf("Using unit gauge links\n");
      lattice.SetGfieldOrd();
    }else if(do_arg.start_conf_kind == START_CONF_DISORD){
      if(!UniqueID()) printf("Using random gauge links\n");
      lattice.SetGfieldDisOrd();
      printf("Gauge checksum = %d\n", lattice.CheckSum());
    }else if(do_arg.start_conf_kind == START_CONF_FILE){    
      if(sprintf(load_config_file,ama_arg.config_fmt,conf) < 0){
	ERR.General("","main()","Configuration filename creation problem : %s | %s",load_config_file,ama_arg.config_fmt);
      }
      //load the configuration
      ReadLatticeParallel readLat;
      if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
      if(dbl_latt_storemode){
	if(!UniqueID()) printf("Disabling U* field reconstruction\n");
	readLat.disableGparityReconstructUstarField();
      }
      readLat.read(lattice,load_config_file);
    }else{
      ERR.General("","main()","Invalid do_arg.start_conf_kind\n");
    }

    if(do_arg.start_seed_kind == START_SEED_FILE){   
      if(sprintf(load_rng_file,ama_arg.rng_fmt,conf) < 0){
	ERR.General("","main()","RNG filename creation problem : %s | %s",load_rng_file,ama_arg.rng_fmt);
      }
      if(UniqueID()==0) printf("Loading RNG state from %s\n",load_rng_file);
      int default_concur=0;
#if TARGET==BGQ
      default_concur=1;
#endif
      LRG.Read(load_rng_file,default_concur);
      if(UniqueID()==0) printf("RNG read.\n");
    }

#ifdef DOUBLE_TLATT
    if(!UniqueID()) printf("Doubling lattice temporal size\n");
    LatticeTimeDoubler doubler;
    doubler.doubleLattice(lattice,do_arg);
#endif


    if(rng_test){ //just load config, generate some uniform random numbers then move onto next config
      if(!UniqueID()) printf("Random offsets in range 0..Lt\n");
      for(int i=0;i<50;i++){
	int offset = int(floor( LRG.Lrand(Lt,0) )) % Lt;
	if(!UniqueID()) printf("%d\n",offset);
      }
      continue;
    }

    if(tshift != 0)
      Tshift4D( (Float*)lattice.GaugeField(), 4*3*3*2, tshift); //do optional temporal shift

    //Gauge fix lattice if required. Do this before the fermion time BCs are applied to the gauge field
    if(ama_arg.fix_gauge.fix_gauge_kind != FIX_GAUGE_NONE){
      Float time = -dclock();
      AlgFixGauge fix_gauge(lattice,&carg,&ama_arg.fix_gauge);
      if(skip_gauge_fix){
	if(!UniqueID()) printf("Skipping gauge fix -> Setting all GF matrices to unity\n");
	gaugeFixUnity(lattice,ama_arg.fix_gauge);      
      }else{
	fix_gauge.run();
      }
      print_time("main","Gauge fix",time + dclock());
    }

    if(lanczos_tune_l || lanczos_tune_h){
      Float time = -dclock();
      doLanczos(lattice, lanczos_tune_l ? lanc_arg_l : lanc_arg_h, GJP.Tbc());
      print_time("main","Lanczos tune",time + dclock());
      if(UniqueID()==0){
	printf("Main job complete\n"); 
	fflush(stdout);
      }
      return 0;
    }

    typedef std::auto_ptr<BFM_Krylov::Lanczos_5d<double> > LanczosPtrType;

    //Generate eigenvectors
#ifdef USE_TBC_INPUT
    LanczosPtrType lanc_l(NULL), lanc_h(NULL);
    Float time = -dclock();
    if(!disable_lanczos) lanc_l = doLanczos(lattice,lanc_arg_l,GJP.Tbc());
    time += dclock();    
    print_time("main","Light quark Lanczos",time);

    time = -dclock();
    if(!disable_lanczos) lanc_h = doLanczos(lattice,lanc_arg_h,GJP.Tbc());
    time += dclock();    
    print_time("main","Heavy quark Lanczos",time);
#else  //USE_TBC_FB
    LanczosPtrType lanc_l_P(NULL), lanc_l_A(NULL), lanc_h_P(NULL), lanc_h_A(NULL);

    Float time = -dclock();
    if(!disable_lanczos) lanc_l_P = doLanczos(lattice,lanc_arg_l,BND_CND_PRD);
    time += dclock();    
    print_time("main","Light quark Lanczos PRD",time);

    time = -dclock();
    if(!disable_lanczos) lanc_l_A = doLanczos(lattice,lanc_arg_l,BND_CND_APRD);
    time += dclock();    
    print_time("main","Light quark Lanczos APRD",time);

    time = -dclock();
    if(!disable_lanczos) lanc_h_P = doLanczos(lattice,lanc_arg_h,BND_CND_PRD);
    time += dclock();    
    print_time("main","Heavy quark Lanczos PRD",time);

    time = -dclock();
    if(!disable_lanczos) lanc_h_A = doLanczos(lattice,lanc_arg_h,BND_CND_APRD);
    time += dclock();    
    print_time("main","Heavy quark Lanczos APRD",time);
#endif

    

    //We want stationary mesons and moving mesons. For GPBC there are two inequivalent directions: along the G-parity axis and perpendicular to it. 
    PropMomContainer props; //stores generated propagators by tag

    std::string results_dir(ama_arg.results_dir);

    std::vector<int> tslice_sloppy(ama_arg.sloppy_solve_timeslices.sloppy_solve_timeslices_val, 
				   ama_arg.sloppy_solve_timeslices.sloppy_solve_timeslices_val + ama_arg.sloppy_solve_timeslices.sloppy_solve_timeslices_len);
    std::vector<int> tslice_exact(ama_arg.exact_solve_timeslices.exact_solve_timeslices_val, 
				  ama_arg.exact_solve_timeslices.exact_solve_timeslices_val + ama_arg.exact_solve_timeslices.exact_solve_timeslices_len);
    
    std::vector<int> bk_tseps(ama_arg.bk_tseps.bk_tseps_val, ama_arg.bk_tseps.bk_tseps_val + ama_arg.bk_tseps.bk_tseps_len);

    if(random_exact_tsrc_offset){
      int offset = int(floor( LRG.Lrand(Lt,0) )) % Lt;
      if(!UniqueID()) printf("Shifting exact src timeslices by offset %d\n",offset);
      for(int i=0;i<tslice_exact.size();i++){
	int nval = (tslice_exact[i]+offset) % Lt;
	if(!UniqueID()) printf("Exact src timeslice %d -> %d\n",tslice_exact[i], nval);
	tslice_exact[i] = nval;
      }
    }

#ifdef USE_TBC_INPUT
    BndCndType tbcs[1] = { GJP.Tbc() };
    LanczosPtrType lanc_l_bcs[1] = { lanc_l }; //note this invalidates the old auto_ptrs
    LanczosPtrType lanc_h_bcs[1] = { lanc_h };
    int ntbc = 1;
    TbcStatus pi_k_tbcuse[1] = { TbcStatus(GJP.Tbc()) };
    int npiktbc = 1;
    TbcStatus bk_tbcuse[2] = { TbcStatus(GJP.Tbc()), TbcStatus(GJP.Tbc()) };
    TbcStatus mres_tbcuse(GJP.Tbc());
    bool store_midprop_l[1] = { true };
#else //USE_TBC_FB
    BndCndType tbcs[2] = { BND_CND_PRD, BND_CND_APRD };
    LanczosPtrType lanc_l_bcs[2] = { lanc_l_P, lanc_l_A };
    LanczosPtrType lanc_h_bcs[2] = { lanc_h_P, lanc_h_A };
    int ntbc = 2;
    TbcStatus pi_k_tbcuse[2] = { TbcStatus(CombinationF), TbcStatus(CombinationB) };
    int npiktbc = 2;
    TbcStatus bk_tbcuse[2] = { TbcStatus(CombinationF), TbcStatus(CombinationB) };
    TbcStatus mres_tbcuse(BND_CND_APRD);
    bool store_midprop_l[2] = { false, true };
#endif

    for(int status = 0; status < 2; status++){ //sloppy, exact
      Float status_start_time = dclock();
      PropPrecision pp = status == 0 ? Sloppy : Exact;
      std::string status_str = status == 0 ? "sloppy" : "exact";

      const std::vector<int> &tslices = status == 0 ? tslice_sloppy : tslice_exact;
      double prec = status == 0 ? ama_arg.sloppy_precision : ama_arg.exact_precision;
      
      if(tslices.size() == 0){
	if(!UniqueID()) printf("Skipping %s contractions because 0 source timeslices given\n",status_str.c_str());
	continue;
      }
      if(!UniqueID()) printf("Starting %s contractions\n",status_str.c_str());

      //Quark inversions
      for(int bci=0;bci<ntbc;bci++){
	std::string bc_name(BndCndType_map[(int)tbcs[bci]].name);
	std::string bc_info = std::string("Propagators Light ") + status_str + std::string(" ") + bc_name;
	time = -dclock();
	quarkInvert(props, Light, pp, prec, ama_arg.ml,tbcs[bci],tslices,light_quark_momenta,store_midprop_l[bci],lattice,lanc_l_bcs[bci].get(),random_prop_solns);
	print_time("main",bc_info.c_str(),time+dclock());
	
	bc_info = std::string("Propagators Heavy ") + status_str + std::string(" ") + bc_name;
	time = -dclock();
	quarkInvert(props, Heavy, pp, prec, ama_arg.mh,tbcs[bci],tslices,heavy_quark_momenta,false,lattice,lanc_h_bcs[bci].get(),random_prop_solns); //no need to store the midpoint prop for the heavy quarks
	print_time("main",bc_info.c_str(),time+dclock());
      }
#ifdef USE_TBC_FB  //Make F and B combinations from P and A computed above
      quarkCombine(props,Light,pp,tslices,light_quark_momenta);
      quarkCombine(props,Heavy,pp,tslices,heavy_quark_momenta);
#endif

      props.printAllTags();
      props.writePropNormTdep(results_dir);

      for(int piktbci = 0; piktbci < npiktbc; piktbci++){
	const TbcStatus & tbs = pi_k_tbcuse[piktbci];

	//Pion 2pt LW functions pseudoscalar and axial sinks	     
	measurePion2ptLW(props,pp,tbs,tslices,pion_momenta,results_dir,conf);
	
	//Pion 2pt WW function pseudoscalar sink
	measurePion2ptPPWW(props,pp,tbs,tslices,pion_momenta,lattice,results_dir,conf);
	
	//SU(2) flavor singlet
	if(GJP.Gparity()) measureLightFlavorSingletLW(props,pp,tbs,tslices,su2_singlet_momenta,results_dir,conf);
	
	//Kaon 2pt LW functions pseudoscalar and axial sinks
	measureKaon2ptLW(props,pp,tbs,tslices,kaon_momenta,results_dir,conf);
      
	//Kaon 2pt WW function pseudoscalar sink
	measureKaon2ptPPWW(props,pp,tbs,tslices,kaon_momenta,lattice,results_dir,conf);
      }

      //BK O_{VV+AA} 3pt contractions
      //Need to ensure that props exist on t0 and t0+tsep for all tseps
      measureBK(props,pp,tslices,bk_tseps,kaon_momenta,bk_tbcuse[0],results_dir,conf,bk_do_flavor_project);

      //J5 and J5q for mres
      measureMres(props,pp,mres_tbcuse,tslices,pion_momenta,results_dir,conf, mres_do_flavor_project);

      props.clear(); //delete all propagators thus far computed
      print_time("main",status_str.c_str(),dclock()-status_start_time);
    }
    lattice.FixGaugeFree();
    print_time("main","Configuration time",dclock()-conf_start_time);
  }//end of conf loop

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  return 0;
}
