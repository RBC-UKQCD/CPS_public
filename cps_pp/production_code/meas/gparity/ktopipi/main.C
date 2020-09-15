//We can save some time by efficiently precomputing SpinColorFlavorMatrix parts for the K->pipi before the inner site loop. However storing these takes a lot of memory and, at least for the DD1 with a 16^3 job
//causes the machine to run out of memory. Use the below options to disable the precompute.
//#define DISABLE_TYPE1_PRECOMPUTE
//#define DISABLE_TYPE2_PRECOMPUTE
//#define DISABLE_TYPE3_PRECOMPUTE
//#define DISABLE_TYPE3_SPLIT_VMV //also disables precompute
//#define DISABLE_TYPE4_PRECOMPUTE

//This option disables the majority of the compute but keeps everything else intact allowing you to test the memory usage without doing a full run
//#define MEMTEST_MODE

#define NODE_DISTRIBUTE_MESONFIELDS //Save memory by keeping meson fields only on single node until needed

//bfm headers
#ifdef USE_BFM
#include<chroma.h>
#include<bfm.h>
#include<util/lattice/bfm_eigcg.h> // This is for the Krylov.h function "matrix_dgemm"
#include<util/lattice/bfm_evo.h>
#endif

//cps headers
#include<alg/common_arg.h>
#include<alg/fix_gauge_arg.h>
#include<alg/do_arg.h>
#include<alg/meas_arg.h>
#include<alg/a2a_arg.h>
#include<alg/lanc_arg.h>
#include<alg/ktopipi_jobparams.h>
#include<util/qioarg.h>
#include<util/ReadLatticePar.h>
#include<alg/alg_fix_gauge.h>
#include<util/flavormatrix.h>
#include<alg/wilson_matrix.h>
#include<util/spincolorflavormatrix.h>

#if defined(USE_GRID) && !defined(DISABLE_GRID_A2A)
#include<util/lattice/fgrid.h>
#endif

#ifdef USE_MPI
//mpi headers
#warning "WARNING : USING MPI"
#include<mpi.h>
#endif

//c++ classes
#include<sys/stat.h>
#include<unistd.h>

#ifdef USE_BFM
using namespace Chroma;
#endif
using namespace cps;

#include <alg/a2a/a2a.h>
#include <alg/a2a/mesonfield.h>
#include <alg/a2a/compute_kaon.h>
#include <alg/a2a/compute_pion.h>
#include <alg/a2a/compute_sigma.h>
#include <alg/a2a/compute_pipi.h>
#include <alg/a2a/compute_ktopipi.h>
#include <alg/a2a/main.h>

#ifdef USE_DESTRUCTIVE_FFT

#ifdef A2A_PREC_DOUBLE
typedef A2ApoliciesDoubleManualAlloc A2Apolicies;
#elif defined(A2A_PREC_SINGLE)
typedef A2ApoliciesSingleManualAlloc A2Apolicies;
#elif defined(A2A_PREC_SIMD_DOUBLE)
typedef A2ApoliciesSIMDdoubleManualAlloc A2Apolicies;
#elif defined(A2A_PREC_SIMD_SINGLE)
typedef A2ApoliciesSIMDsingleManualAlloc A2Apolicies;
#else
#error "Must provide an A2A precision"
#endif

#else

#ifdef A2A_PREC_DOUBLE
typedef A2ApoliciesDoubleAutoAlloc A2Apolicies;
#elif defined(A2A_PREC_SINGLE)
typedef A2ApoliciesSingleAutoAlloc A2Apolicies;
#elif defined(A2A_PREC_SIMD_DOUBLE)
typedef A2ApoliciesSIMDdoubleAutoAlloc A2Apolicies;
#elif defined(A2A_PREC_SIMD_SINGLE)
typedef A2ApoliciesSIMDsingleAutoAlloc A2Apolicies;
#else
#error "Must provide an A2A precision"
#endif

#endif


int main (int argc,char **argv )
{
  Start(&argc, &argv);
  if(!UniqueID()){ printf("Arguments:\n"); fflush(stdout); }
  for(int i=0;i<argc;i++){
    if(!UniqueID()){ printf("%d \"%s\"\n",i,argv[i]); fflush(stdout); }
  }
#ifdef NODE_DISTRIBUTE_MESONFIELDS
  if(!UniqueID()) printf("Using node distribution of meson fields\n");
#endif
#ifdef MEMTEST_MODE
  if(!UniqueID()) printf("Running in MEMTEST MODE (so don't expect useful results)\n");
#endif

  const char *cname=argv[0];
  const int TrajStart = atoi(argv[2]);
  const int LessThanLimit = atoi(argv[3]);

  int nthreads = 1;
#if TARGET == BGQ
  nthreads = 64;
#endif
  bool randomize_vw = false; //rather than doing the Lanczos and inverting the propagators, etc, just use random vectors for V and W. This should only be used after you have tested
  bool tune_lanczos_light = false; //just run the light lanczos on first config then exit
  bool tune_lanczos_heavy = false; //just run the heavy lanczos on first config then exit
  bool skip_gauge_fix = false;
  bool double_latt = false; //most ancient 8^4 quenched lattices stored both U and U*. Enable this to read those configs
  bool mixed_solve = true; //do high mode inversions using mixed precision solves. Is disabled if we turn off the single-precision conversion of eigenvectors (because internal single-prec inversion needs singleprec eigenvectors)
  bool evecs_single_prec = true; //convert the eigenvectors to single precision to save memory

  const int ngrid_arg = 7;
  const std::string grid_args[ngrid_arg] = { "--debug-signals", "--dslash-generic", "--dslash-unroll", "--dslash-asm", "--shm", "--lebesgue", "--cacheblocking" };
  const int grid_args_skip[ngrid_arg] = { 1, 1, 1, 1, 2, 1, 2 };
  
  int arg = 4;
  while(arg < argc){
    char* cmd = argv[arg];
    if( strncmp(cmd,"-nthread",8) == 0){
      if(arg == argc-1){ 
	if(!UniqueID()){ printf("-nthread must be followed by a number!\n"); fflush(stdout); }
	exit(-1);
      }
      std::stringstream ss; ss << argv[arg+1];
      ss >> nthreads;
      if(!UniqueID()){ printf("Setting number of threads to %d\n",nthreads); }
      arg+=2;
    }else if( strncmp(cmd,"-randomize_vw",15) == 0){
      randomize_vw = true;
      if(!UniqueID()){ printf("Using random vectors for V and W, skipping Lanczos and inversion stages\n"); fflush(stdout); }
      arg++;
    }else if( strncmp(cmd,"-tune_lanczos_light",15) == 0){
      tune_lanczos_light = true;
      if(!UniqueID()){ printf("Just tuning light lanczos on first config\n"); fflush(stdout); }
      arg++;
    }else if( strncmp(cmd,"-tune_lanczos_heavy",15) == 0){
      tune_lanczos_heavy = true;
      if(!UniqueID()){ printf("Just tuning heavy lanczos on first config\n"); fflush(stdout); }
      arg++;
    }else if( strncmp(cmd,"-double_latt",15) == 0){
      double_latt = true;
      if(!UniqueID()){ printf("Loading doubled lattices\n"); fflush(stdout); }
      arg++;
    }else if( strncmp(cmd,"-skip_gauge_fix",20) == 0){
      skip_gauge_fix = true;
      if(!UniqueID()){ printf("Skipping gauge fixing\n"); fflush(stdout); }
      arg++;
    }else if( strncmp(cmd,"-disable_evec_singleprec_convert",30) == 0){
      evecs_single_prec = false;
      mixed_solve = false;
      if(!UniqueID()){ printf("Disabling single precision conversion of evecs\n"); fflush(stdout); }
      arg++;
    }else if( strncmp(cmd,"-disable_mixed_prec_CG",30) == 0){
      mixed_solve = false;
      if(!UniqueID()){ printf("Disabling mixed-precision CG\n"); fflush(stdout); }
      arg++;
    }else{
      bool is_grid_arg = false;
      for(int i=0;i<ngrid_arg;i++){
	if( std::string(cmd) == grid_args[i] ){
	  if(!UniqueID()){ printf("main.C: Ignoring Grid argument %s\n",cmd); fflush(stdout); }
	  arg += grid_args_skip[i];
	  is_grid_arg = true;
	  break;
	}
      }
      if(!is_grid_arg){
	if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
	exit(-1);
      }
    }
  }
  const char *fname="main(int,char**)";
  
#ifdef A2A_LANCZOS_SINGLE
  if(!evecs_single_prec) ERR.General("",fname,"Must use single-prec eigenvectors when doing Lanczos in single precision\n");
#endif
  
  if(chdir(argv[1])!=0) ERR.General("",fname,"Unable to switch to directory '%s'\n",argv[1]);
  CommonArg common_arg("",""), common_arg2("","");
  DoArg do_arg;
  JobParams jp;
  MeasArg meas_arg;
  FixGaugeArg fix_gauge_arg;
  A2AArg a2a_arg, a2a_arg_s;
  LancArg lanc_arg, lanc_arg_s;

  if(!do_arg.Decode("do_arg.vml","do_arg")){
    do_arg.Encode("do_arg.templ","do_arg");
    VRB.Result(cname,fname,"Can't open do_arg.vml!\n");exit(1);
  }
  if(!jp.Decode("job_params.vml","job_params")){
    jp.Encode("job_params.templ","job_params");
    VRB.Result(cname,fname,"Can't open job_params.vml!\n");exit(1);
  }
  if(!meas_arg.Decode("meas_arg.vml","meas_arg")){
    meas_arg.Encode("meas_arg.templ","meas_arg");
    std::cout<<"Can't open meas_arg!"<<std::endl;exit(1);
  }
  if(!a2a_arg.Decode("a2a_arg.vml","a2a_arg")){
    a2a_arg.Encode("a2a_arg.templ","a2a_arg");
    VRB.Result(cname,fname,"Can't open a2a_arg.vml!\n");exit(1);
  }
  if(!a2a_arg_s.Decode("a2a_arg_s.vml","a2a_arg_s")){
    a2a_arg_s.Encode("a2a_arg_s.templ","a2a_arg_s");
    VRB.Result(cname,fname,"Can't open a2a_arg_s.vml!\n");exit(1);
  }
  if(!lanc_arg.Decode("lanc_arg.vml","lanc_arg")){
    lanc_arg.Encode("lanc_arg.templ","lanc_arg");
    VRB.Result(cname,fname,"Can't open lanc_arg.vml!\n");exit(1);
  }
  if(!lanc_arg_s.Decode("lanc_arg_s.vml","lanc_arg_s")){
    lanc_arg_s.Encode("lanc_arg_s.templ","lanc_arg_s");
    VRB.Result(cname,fname,"Can't open lanc_arg_s.vml!\n");exit(1);
  }
  if(!fix_gauge_arg.Decode("fix_gauge_arg.vml","fix_gauge_arg")){
    fix_gauge_arg.Encode("fix_gauge_arg.templ","fix_gauge_arg");
    VRB.Result(cname,fname,"Can't open fix_gauge_arg.vml!\n");exit(1);
  }

  //Chroma::initialize(&argc,&argv);
  GJP.Initialize(do_arg);
  LRG.Initialize();

  if(double_latt) SerialIO::dbl_latt_storemode = true;

  if(!UniqueID()) printf("Initial memory post-initialize:\n");
  printMem();


  if(!tune_lanczos_light && !tune_lanczos_heavy){ 
    assert(a2a_arg.nl <= lanc_arg.N_true_get);
    assert(a2a_arg_s.nl <= lanc_arg_s.N_true_get);
  }
#ifdef USE_BFM
  cps_qdp_init(&argc,&argv);
  Chroma::initialize(&argc,&argv);
#endif

  int Lt = GJP.Tnodes()*GJP.TnodeSites();
  
  LatticeSolvers solvers(jp,nthreads); //for BFM holds a double and single precision bfm instance, nothing for Grid
  
  if(chdir(meas_arg.WorkDirectory)!=0) ERR.General("",fname,"Unable to switch to work directory '%s'\n",meas_arg.WorkDirectory);
  double time;

  if(!UniqueID()) printf("Memory prior to config loop:\n");
  printMem();

#if defined(USE_GRID_A2A) || defined(USE_GRID_LANCZOS)
  typedef A2Apolicies LanczosPolicies;
  if(GJP.Gparity()){
#ifndef USE_GRID_GPARITY
    ERR.General("","","Must compile main program with flag USE_GRID_GPARITY to enable G-parity\n");
#endif
  }else{
#ifdef USE_GRID_GPARITY
    ERR.General("","","Must compile main program with flag USE_GRID_GPARITY off to disable G-parity\n");
#endif
  }      
#else
  typedef void LanczosPolicies;
#endif

  //Setup parameters of fields
  typedef typename A2Apolicies::ComplexType mf_Complex;
  typedef typename A2Apolicies::LatticeType LatticeType;
  typedef typename A2Apolicies::SourcePolicies::DimensionPolicy::ParamType Field3DparamType;
  typedef typename A2Apolicies::FermionFieldType::InputParamType Field4DparamType;
  Field4DparamType field4dparams; setupFieldParams<mf_Complex>(field4dparams);
  Field3DparamType field3dparams; setupFieldParams<mf_Complex>(field3dparams);
  
  //-------------------- Main Loop Begin! -------------------- //
  for(int conf = TrajStart; conf < LessThanLimit; conf += meas_arg.TrajIncrement) {
    double conf_time = -dclock();
    if(!UniqueID()) std::cout<<"Starting configuration "<<conf<< std::endl;

    meas_arg.TrajCur = conf;

    std::string dir(meas_arg.WorkDirectory);
    common_arg.set_filename(dir.c_str());

    //-------------------- Read gauge field --------------------//
    ReadGaugeField(meas_arg,double_latt); 
    ReadRngFile(meas_arg,double_latt); 

    if(!UniqueID()) printf("Memory after gauge and RNG read:\n");
    printMem();

    LatticeSetup<LatticeType> lattice_setup(jp,solvers); //for BFM this creates a lattice object and imports the gauge field into the bfm instances, for Grid a lattice object and import of the gauge field
    LatticeType &lat = lattice_setup.getLattice();

    if(tune_lanczos_light || tune_lanczos_heavy){
      if(!UniqueID()) printf("Tuning lanczos %s with mass %f\n", tune_lanczos_light ? "light": "heavy", tune_lanczos_light ? lanc_arg.mass : lanc_arg_s.mass);
      time = -dclock();
      Lanczos<LanczosPolicies> eig;
      eig.compute(tune_lanczos_light ? lanc_arg : lanc_arg_s, solvers, lat);
      time += dclock();
      print_time("main","Lanczos",time);
      exit(0);
    }

    //-------------------- Light quark v and w --------------------//
    Lanczos<LanczosPolicies> eig;

    if(!randomize_vw){
      if(!UniqueID()) printf("Running light quark Lanczos\n");
      time = -dclock();
      eig.compute(lanc_arg, solvers, lat);
      time += dclock();
      print_time("main","Light quark Lanczos",time);

      if(!UniqueID()) printf("Memory after light quark Lanczos:\n");
      printMem();      

#ifndef A2A_LANCZOS_SINGLE
      if(evecs_single_prec){
	eig.toSingle();
	if(!UniqueID()) printf("Memory after single-prec conversion of light quark evecs:\n");
	printMem();
      }
#endif
    }

    if(!UniqueID()) printf("Computing light quark A2A vectors\n");
    time = -dclock();

    if(!UniqueID()) printf("V vector requires %f MB, W vector %f MB of memory\n", 
			   A2AvectorV<A2Apolicies>::Mbyte_size(a2a_arg,field4dparams), A2AvectorW<A2Apolicies>::Mbyte_size(a2a_arg,field4dparams) );
    
    A2AvectorV<A2Apolicies> V(a2a_arg, field4dparams);
    A2AvectorW<A2Apolicies> W(a2a_arg, field4dparams);

#ifdef USE_DESTRUCTIVE_FFT
    V.allocModes(); W.allocModes();
#endif
    
    if(!randomize_vw){
      computeA2Avectors<A2Apolicies,LanczosPolicies>::compute(V,W,mixed_solve,evecs_single_prec, lat, eig, solvers);
      //W.computeVW(V, lat, *eig.eig, evecs_single_prec, solvers.dwf_d, mixed_solve ? & solvers.dwf_f : NULL);
    }else randomizeVW<A2Apolicies>(V,W);    

    if(!UniqueID()) printf("Memory after light A2A vector computation:\n");
    printMem();

    eig.freeEvecs();

    if(!UniqueID()) printf("Memory after light evec free:\n");
    printMem();

    time += dclock();
    print_time("main","Light quark A2A vectors",time);

    //-------------------- Strange quark v and w --------------------//
    Lanczos<LanczosPolicies> eig_s;

    if(!randomize_vw){
      if(!UniqueID()) printf("Running strange quark Lanczos\n");
      time = -dclock();
      eig_s.compute(lanc_arg_s,solvers,lat);
      time += dclock();
      print_time("main","Strange quark Lanczos",time);

      if(!UniqueID()) printf("Memory after heavy quark Lanczos:\n");
      printMem();

#ifndef A2A_LANCZOS_SINGLE
      if(evecs_single_prec){
	eig_s.toSingle();
	if(!UniqueID()) printf("Memory after single-prec conversion of heavy quark evecs:\n");
	printMem();
      }
#endif
    }

    if(!UniqueID()) printf("Computing strange quark A2A vectors\n");
    time = -dclock();

    if(!UniqueID()) printf("V_s vector requires %f MB, W_s vector %f MB of memory\n", 
			   A2AvectorV<A2Apolicies>::Mbyte_size(a2a_arg_s,field4dparams), A2AvectorW<A2Apolicies>::Mbyte_size(a2a_arg_s,field4dparams) );

    A2AvectorV<A2Apolicies> V_s(a2a_arg_s,field4dparams);
    A2AvectorW<A2Apolicies> W_s(a2a_arg_s,field4dparams);

#ifdef USE_DESTRUCTIVE_FFT
    V_s.allocModes(); W_s.allocModes();
#endif
    
    if(!randomize_vw){
      computeA2Avectors<A2Apolicies,LanczosPolicies>::compute(V_s,W_s,mixed_solve,evecs_single_prec, lat, eig_s, solvers);
      //W_s.computeVW(V_s, lat, *eig_s.eig, evecs_single_prec, solvers.dwf_d, mixed_solve ? & solvers.dwf_f : NULL);
    }else randomizeVW<A2Apolicies>(V_s,W_s);      

    if(!UniqueID()) printf("Memory after heavy A2A vector computation:\n");
    printMem();

    eig_s.freeEvecs();

    if(!UniqueID()) printf("Memory after heavy evec free:\n");
    printMem();

    time += dclock();
    print_time("main","Strange quark A2A vectors",time);

    //-------------------Fix gauge----------------------------
    AlgFixGauge fix_gauge(lat,&common_arg,&fix_gauge_arg);
    if(skip_gauge_fix){
      if(!UniqueID()) printf("Skipping gauge fix -> Setting all GF matrices to unity\n");
      gaugeFixUnity(lat,fix_gauge_arg);      
    }else{
      if(!UniqueID()){ printf("Gauge fixing\n"); fflush(stdout); }
      time = -dclock();
#ifndef MEMTEST_MODE
      fix_gauge.run();
#else
      gaugeFixUnity(lat,fix_gauge_arg);
#endif      
      time += dclock();
      print_time("main","Gauge fix",time);
    }

    if(!UniqueID()) printf("Memory after gauge fix:\n");
    printMem();

    //-------------------------Compute the kaon two-point function---------------------------------
    {
      if(!UniqueID()) printf("Computing kaon 2pt function\n");
      time = -dclock();
      fMatrix<typename A2Apolicies::ScalarComplexType> kaon(Lt,Lt);
      ComputeKaon<A2Apolicies>::compute(kaon,
				     W, V, W_s, V_s,
					jp.kaon_rad, lat, field3dparams);
      std::ostringstream os; os << meas_arg.WorkDirectory << "/traj_" << conf << "_kaoncorr";
      kaon.write(os.str());
      time += dclock();
      print_time("main","Kaon 2pt function",time);
    }

    if(!UniqueID()) printf("Memory after kaon 2pt function computation:\n");
    printMem();

    //The pion two-point function and pipi/k->pipi all utilize the same meson fields. Generate those here
    //For convenience pointers to the meson fields are collected into a single object that is passed to the compute methods
    RequiredMomentum<StandardPionMomentaPolicy> pion_mom; //these are the W and V momentum combinations

    MesonFieldMomentumContainer<A2Apolicies> mf_ll_con; //stores light-light meson fields, accessible by momentum
    MesonFieldMomentumContainer<A2Apolicies> mf_ll_con_2s; //Gparity only
    
    if(!UniqueID()) printf("Computing light-light meson fields\n");
    time = -dclock();
    if(!GJP.Gparity()) ComputePion<A2Apolicies>::computeMesonFields(mf_ll_con, meas_arg.WorkDirectory,conf, pion_mom, W, V, jp.pion_rad, lat, field3dparams);
    else ComputePion<A2Apolicies>::computeGparityMesonFields(mf_ll_con, mf_ll_con_2s, meas_arg.WorkDirectory,conf, pion_mom, W, V, jp.pion_rad, lat, field3dparams);
    time += dclock();
    print_time("main","Light-light meson fields",time);

    if(!UniqueID()) printf("Memory after light-light meson field computation:\n");
    printMem();

    //----------------------------Compute the pion two-point function---------------------------------
    int nmom = pion_mom.nMom();

    if(!UniqueID()) printf("Computing pion 2pt function\n");
    time = -dclock();
    for(int p=0;p<nmom;p+=2){ //note odd indices 1,3,5 etc have equal and opposite momenta to 0,2,4... 
      if(!UniqueID()) printf("Starting pidx %d\n",p);
      fMatrix<typename A2Apolicies::ScalarComplexType> pion(Lt,Lt);
      ComputePion<A2Apolicies>::compute(pion, mf_ll_con, pion_mom, p);
      //Note it seems Daiqian's pion momenta are opposite what they should be for 'conventional' Fourier transform phase conventions:
      //f'(p) = \sum_{x,y}e^{ip(x-y)}f(x,y)  [conventional]
      //f'(p) = \sum_{x,y}e^{-ip(x-y)}f(x,y) [Daiqian]
      //This may have been a mistake as it only manifests in the difference between the labelling of the pion momenta and the sign of 
      //the individual quark momenta.
      //However it doesn't really make any difference. If you enable DAIQIAN_PION_PHASE_CONVENTION
      //the output files will be labelled in Daiqian's convention
#define DAIQIAN_PION_PHASE_CONVENTION

      std::ostringstream os; os << meas_arg.WorkDirectory << "/traj_" << conf << "_pioncorr_mom";
#ifndef DAIQIAN_PION_PHASE_CONVENTION
      os << pion_mom.getMesonMomentum(p).file_str(2);  //note the divisor of 2 is to put the momenta in units of pi/L and not pi/2L
#else
      os << (-pion_mom.getMesonMomentum(p)).file_str(2);
#endif
      pion.write(os.str());
    }
    time += dclock();
    print_time("main","Pion 2pt function",time);

    if(!UniqueID()) printf("Memory after pion 2pt function computation:\n");
    printMem();

    time = -dclock();
    if(!UniqueID()) printf("Computing sigma mesonfield computation\n");
    ComputeSigma<A2Apolicies>::computeAndWrite(meas_arg.WorkDirectory,conf,W,V, jp.pion_rad, lat, field3dparams);
    time += dclock();
    print_time("main","Sigma meson fields ",time);
    
    //------------------------------I=0 and I=2 PiPi two-point function---------------------------------
    if(!UniqueID()) printf("Computing pi-pi 2pt function\n");
    double timeC(0), timeD(0), timeR(0), timeV(0);
    double* timeCDR[3] = {&timeC, &timeD, &timeR};

    for(int psrcidx=0; psrcidx < nmom; psrcidx++){
      ThreeMomentum p_pi1_src = pion_mom.getMesonMomentum(psrcidx);

      for(int psnkidx=0; psnkidx < nmom; psnkidx++){	
	fMatrix<typename A2Apolicies::ScalarComplexType> pipi(Lt,Lt);
	ThreeMomentum p_pi1_snk = pion_mom.getMesonMomentum(psnkidx);
	
	MesonFieldProductStore<A2Apolicies> products; //try to reuse products of meson fields wherever possible

	char diag[3] = {'C','D','R'};
	for(int d = 0; d < 3; d++){
	  if(!UniqueID()){ printf("Doing pipi figure %c, psrcidx=%d psnkidx=%d\n",diag[d],psrcidx,psnkidx); fflush(stdout); }

	  time = -dclock();
	  ComputePiPiGparity<A2Apolicies>::compute(pipi, diag[d], p_pi1_src, p_pi1_snk, jp.pipi_separation, jp.tstep_pipi, mf_ll_con, products);
	  std::ostringstream os; os << meas_arg.WorkDirectory << "/traj_" << conf << "_Figure" << diag[d] << "_sep" << jp.pipi_separation;
#ifndef DAIQIAN_PION_PHASE_CONVENTION
	  os << "_mom" << p_pi1_src.file_str(2) << "_mom" << p_pi1_snk.file_str(2);
#else
	  os << "_mom" << (-p_pi1_src).file_str(2) << "_mom" << (-p_pi1_snk).file_str(2);
#endif
	  pipi.write(os.str());
	  time += dclock();
	  *timeCDR[d] += time;
	}
      }

      { //V diagram
	if(!UniqueID()){ printf("Doing pipi figure V, pidx=%d\n",psrcidx); fflush(stdout); }
	time = -dclock();
	fVector<typename A2Apolicies::ScalarComplexType> figVdis(Lt);
	ComputePiPiGparity<A2Apolicies>::computeFigureVdis(figVdis,p_pi1_src,jp.pipi_separation,mf_ll_con);
	std::ostringstream os; os << meas_arg.WorkDirectory << "/traj_" << conf << "_FigureVdis_sep" << jp.pipi_separation;
#ifndef DAIQIAN_PION_PHASE_CONVENTION
	os << "_mom" << p_pi1_src.file_str(2);
#else
	os << "_mom" << (-p_pi1_src).file_str(2);
#endif
	figVdis.write(os.str());
	time += dclock();
	timeV += time;
      }
    }//end of psrcidx loop

    print_time("main","Pi-pi figure C",timeC);
    print_time("main","Pi-pi figure D",timeD);
    print_time("main","Pi-pi figure R",timeR);
    print_time("main","Pi-pi figure V",timeV);

    if(!UniqueID()) printf("Memory after pi-pi 2pt function computation:\n");
    printMem();

    //--------------------------------------K->pipi contractions--------------------------------------------------------
    //We first need to generate the light-strange W*W contraction
    std::vector<A2AmesonField<A2Apolicies,A2AvectorWfftw,A2AvectorWfftw> > mf_ls_ww;
    ComputeKtoPiPiGparity<A2Apolicies>::generatelsWWmesonfields(mf_ls_ww,W,W_s,jp.kaon_rad,lat, field3dparams);

    std::vector<int> k_pi_separation(jp.k_pi_separation.k_pi_separation_len);
    for(int i=0;i<jp.k_pi_separation.k_pi_separation_len;i++) k_pi_separation[i] = jp.k_pi_separation.k_pi_separation_val[i];

    if(!UniqueID()) printf("Memory after computing W*W meson fields:\n");
    printMem();

    typedef ComputeKtoPiPiGparity<A2Apolicies>::ResultsContainerType ResultsContainerType;
    typedef ComputeKtoPiPiGparity<A2Apolicies>::MixDiagResultsContainerType MixDiagResultsContainerType;

    MesonFieldMomentumContainer<A2Apolicies>* ll_meson_field_ptrs[2] = { &mf_ll_con, &mf_ll_con_2s };
    const int nsource = GJP.Gparity() ? 2 : 1;
    const std::string src_str[2] = { "", "_src2s" };
    
    //For type1 loop over momentum of pi1 (conventionally the pion closest to the kaon)
    int ngp = 0; for(int i=0;i<3;i++) if(GJP.Bc(i)==BND_CND_GPARITY) ngp++;
#define TYPE1_DO_ASSUME_ROTINVAR_GP3  //For GPBC in 3 directions we can assume rotational invariance around the G-parity diagonal vector (1,1,1) and therefore calculate only one off-diagonal momentum

    if(!UniqueID()) printf("Starting type 1 contractions, nmom = %d\n",nmom);
    time = -dclock();
    
    for(int pidx=0; pidx < nmom; pidx++){
#ifdef TYPE1_DO_ASSUME_ROTINVAR_GP3
      if(ngp == 3 && pidx >= 4) continue; // p_pi1 = (-1,-1,-1), (1,1,1) [diag] (1,-1,-1), (-1,1,1) [orth] only
#endif
      for(int sidx=0; sidx<nsource;sidx++){
      
	if(!UniqueID()) printf("Starting type 1 contractions with pidx=%d and source idx %d\n",pidx,sidx);
	if(!UniqueID()) printf("Memory status before type1 K->pipi:\n");
	printMem();

	ThreeMomentum p_pi1 = pion_mom.getMesonMomentum(pidx);
	std::vector<ResultsContainerType> type1;
	ComputeKtoPiPiGparity<A2Apolicies>::type1(type1,
						  k_pi_separation, jp.pipi_separation, jp.tstep_type12, jp.xyzstep_type1, p_pi1,
						  mf_ls_ww, *ll_meson_field_ptrs[sidx],
						  V, V_s,
						  W, W_s);
	for(int kpi_idx=0;kpi_idx<k_pi_separation.size();kpi_idx++){
	  std::ostringstream os; os << meas_arg.WorkDirectory << "/traj_" << conf << "_type1_deltat_" << k_pi_separation[kpi_idx] << src_str[sidx] << "_sep_" << jp.pipi_separation;
#ifndef DAIQIAN_PION_PHASE_CONVENTION
	  os << "_mom" << p_pi1.file_str(2);
#else
	  os << "_mom" << (-p_pi1).file_str(2);
#endif
	  type1[kpi_idx].write(os.str());
	}
	if(!UniqueID()) printf("Memory status after type1 K->pipi:\n");
	printMem();
      }
    }

    
    time += dclock();
    print_time("main","K->pipi type 1",time);

    if(!UniqueID()) printf("Memory after type1 K->pipi:\n");
    printMem();

    //Type 2 and 3 are optimized by performing the sum over pipi momentum orientations within the contraction
    time = -dclock();    
    for(int sidx=0; sidx< nsource; sidx++){
      if(!UniqueID()) printf("Starting type 2 contractions with source idx %d\n", sidx);
      std::vector<ResultsContainerType> type2;
      ComputeKtoPiPiGparity<A2Apolicies>::type2(type2,
					     k_pi_separation, jp.pipi_separation, jp.tstep_type12, pion_mom,
					     mf_ls_ww, *ll_meson_field_ptrs[sidx],
					     V, V_s,
					     W, W_s);
      for(int kpi_idx=0;kpi_idx<k_pi_separation.size();kpi_idx++){
	std::ostringstream os; os << meas_arg.WorkDirectory << "/traj_" << conf << "_type2_deltat_" << k_pi_separation[kpi_idx] << src_str[sidx] << "_sep_" << jp.pipi_separation;
	type2[kpi_idx].write(os.str());
      }
    }
    time += dclock();
    print_time("main","K->pipi type 2",time);
    
    if(!UniqueID()) printf("Memory after type2 K->pipi:\n");
    printMem();
    

    time = -dclock();
    for(int sidx=0; sidx< nsource; sidx++){
      if(!UniqueID()) printf("Starting type 3 contractions with source idx %d\n", sidx);
      std::vector<ResultsContainerType> type3;
      std::vector<MixDiagResultsContainerType> mix3;
      ComputeKtoPiPiGparity<A2Apolicies>::type3(type3,mix3,
					     k_pi_separation, jp.pipi_separation, 1, pion_mom,
					     mf_ls_ww, *ll_meson_field_ptrs[sidx],
					     V, V_s,
					     W, W_s);
      for(int kpi_idx=0;kpi_idx<k_pi_separation.size();kpi_idx++){
	std::ostringstream os; os << meas_arg.WorkDirectory << "/traj_" << conf << "_type3_deltat_" << k_pi_separation[kpi_idx] << src_str[sidx] << "_sep_" << jp.pipi_separation;
	write(os.str(),type3[kpi_idx],mix3[kpi_idx]);
      }
    }
    time += dclock();
    print_time("main","K->pipi type 3",time);
    
    if(!UniqueID()) printf("Memory after type3 K->pipi:\n");
    printMem();
    

    {
      //Type 4 has no momentum loop as the pion disconnected part is computed as part of the pipi 2pt function calculation
      time = -dclock();
      if(!UniqueID()) printf("Starting type 4 contractions\n");
      ResultsContainerType type4;
      MixDiagResultsContainerType mix4;
      
      ComputeKtoPiPiGparity<A2Apolicies>::type4(type4, mix4,
					     1,
					     mf_ls_ww,
					     V, V_s,
					     W, W_s);
      
      {
	std::ostringstream os; os << meas_arg.WorkDirectory << "/traj_" << conf << "_type4";
	write(os.str(),type4,mix4);
      }
      time += dclock();
      print_time("main","K->pipi type 4",time);
    
      if(!UniqueID()) printf("Memory after type4 K->pipi and end of config loop:\n");
      printMem();
    }

    conf_time += dclock();
    print_time("main","Configuration total",conf_time);
  }//end of config loop

  if(!UniqueID()) printf("Done\n");
  End();
}

