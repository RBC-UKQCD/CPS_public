//#define USE_GRID
#define USE_GRID_A2A
#define USE_GRID_LANCZOS
#include<chroma.h>

#ifdef USE_CALLGRIND
#include<valgrind/callgrind.h>
#else
#define CALLGRIND_START_INSTRUMENTATION ;
#define CALLGRIND_STOP_INSTRUMENTATION ;
#define CALLGRIND_TOGGLE_COLLECT ;
#endif

#ifdef USE_VTUNE
#include<ittnotify.h>
#else
void __itt_pause(){}
void __itt_resume(){}
void __itt_detach(){}
#endif

//bfm headers
#ifdef USE_BFM
#include<bfm.h>
#include<util/lattice/bfm_eigcg.h> // This is for the Krylov.h function "matrix_dgemm"
#include<util/lattice/bfm_evo.h>
#endif

//cps headers
#include<util/time_cps.h>
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

//using namespace Chroma;
using namespace cps;

#include <alg/a2a/template_wizardry.h>
#include <alg/a2a/spin_color_matrices.h>
#include <alg/a2a/a2a.h>
#include <alg/a2a/mesonfield.h>
#include <alg/a2a/compute_ktopipi_base.h>

#include "benchmark_mesonfield.h"

typedef A2ApoliciesSIMDdoubleAutoAlloc GridA2Apolicies;
typedef A2ApoliciesDoubleAutoAlloc ScalarA2Apolicies;

typedef typename GridA2Apolicies::ComplexType grid_Complex;
typedef typename ScalarA2Apolicies::ComplexType mf_Complex;
typedef typename mf_Complex::value_type mf_Float;


int main(int argc,char *argv[])
{
  Start(&argc, &argv);
  int ngp;
  { std::stringstream ss; ss << argv[1]; ss >> ngp; }

  if(!UniqueID()) printf("Doing G-parity in %d directions\n",ngp);

  bool save_config(false);
  bool load_config(false);
  bool load_lrg(false);
  bool save_lrg(false);
  char *load_config_file;
  char *save_config_file;
  char *save_lrg_file;
  char *load_lrg_file;
  bool verbose(false);
  bool unit_gauge(false);

  int size[] = {2,2,2,2,2};
  int nthreads = 1;
  int ntests = 10;
  
  double tol = 1e-8;
  int nlowmodes = 100;
  printf("Argc is %d\n",argc);
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
      for(int d=0;d<5;d++)
	size[d] = toInt(argv[i+1+d]);
      i+=6;
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
    }else if( strncmp(cmd,"-verbose",15) == 0){
      verbose=true;
      i++;
    }else if( strncmp(cmd,"-nthread",15) == 0){
      nthreads = toInt(argv[i+1]);
      printf("Set nthreads to %d\n", nthreads);
      i+=2;
    }else if( strncmp(cmd,"-ntest",15) == 0){
      ntests = toInt(argv[i+1]);
      printf("Set ntests to %d\n", ntests);
      i+=2;
    }else if( strncmp(cmd,"-unit_gauge",15) == 0){
      unit_gauge=true;
      i++;
    }else if( strncmp(cmd,"-tolerance",15) == 0){
      std::stringstream ss; ss << argv[i+1];
      ss >> tol;
      if(!UniqueID()) printf("Set tolerance to %g\n",tol);
      i+=2;
    }else if( strncmp(cmd,"-nl",15) == 0){
      std::stringstream ss; ss << argv[i+1];
      ss >> nlowmodes;
      if(!UniqueID()) printf("Set nl to %d\n",nlowmodes);
      i+=2;
    }else if( strncmp(cmd,"-mf_outerblocking",15) == 0){
      int* b[3] = { &BlockedMesonFieldArgs::bi, &BlockedMesonFieldArgs::bj, &BlockedMesonFieldArgs::bp };
      for(int a=0;a<3;a++){
	std::stringstream ss; ss << argv[i+1+a];
	ss >> *b[a];
      }
      i+=4;
    }else if( strncmp(cmd,"-mf_innerblocking",15) == 0){
      int* b[3] = { &BlockedMesonFieldArgs::bii, &BlockedMesonFieldArgs::bjj, &BlockedMesonFieldArgs::bpp };
      for(int a=0;a<3;a++){
	std::stringstream ss; ss << argv[i+1+a];
	ss >> *b[a];
      }
      i+=4;
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }

  printf("Lattice size is %d %d %d %d\n",size[0],size[1],size[2],size[3],size[4]);

  CommonArg common_arg;
  DoArg do_arg;  setupDoArg(do_arg,size,ngp,verbose);

  GJP.Initialize(do_arg);
  GJP.SetNthreads(nthreads);
  
#if TARGET == BGQ
  LRG.setSerial();
#endif
  LRG.Initialize(); //usually initialised when lattice generated, but I pre-init here so I can load the state from file

  GnoneFnone lattice;

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }					       
  if(!load_config){
    printf("Creating gauge field\n");
    if(!unit_gauge) lattice.SetGfieldDisOrd();
    else lattice.SetGfieldOrd();
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

  A2AArg a2a_args;
  a2a_args.nl = nlowmodes;
  a2a_args.nhits = 1;
  a2a_args.rand_type = UONE;
  a2a_args.src_width = 1;

  if(0) testCyclicPermute();
  
  if(0) demonstrateFFTreln<ScalarA2Apolicies>(a2a_args);


  if(0) testA2AvectorFFTrelnGparity<ScalarA2Apolicies>(a2a_args,lattice);
  if(0) testA2AvectorFFTrelnGparity<GridA2Apolicies>(a2a_args,lattice);

  if(0) testMultiSource<ScalarA2Apolicies>(a2a_args,lattice);

  if(0) testMfFFTreln<ScalarA2Apolicies>(a2a_args,lattice);
  if(0) testMfFFTreln<GridA2Apolicies>(a2a_args,lattice);

  if(0) testFFTopt<ScalarA2Apolicies>();
  if(0) testFFTopt<GridA2Apolicies>();

  if(0) testA2AFFTinv<ScalarA2Apolicies>(a2a_args,lattice);
  
  if(0) testVVdag<ScalarA2Apolicies>(lattice);
  if(0) testVVdag<GridA2Apolicies>(lattice);

  if(0) testDestructiveFFT<A2ApoliciesDoubleManualAlloc>(a2a_args,lattice);
  
  if(0) testA2AallocFree(a2a_args,lattice);


  if(0) benchmarkMFcontractKernel<GridA2Apolicies>(ntests,nthreads);
  
  if(1) testGridg5Contract<Grid::vComplexD>(); //Keep active because its very quick and checks that the kernel is correct

  if(0) benchmarkTrace(ntests,tol);
  if(0) benchmarkSpinFlavorTrace(ntests,tol);
  if(0) benchmarkTraceProd(ntests,tol);
  if(0) benchmarkColorTranspose(ntests,tol);
  if(0) benchmarkmultGammaLeft(ntests, tol);
  
  if(1) testMFcontract<ScalarA2Apolicies,GridA2Apolicies>(a2a_args, nthreads,tol);
  if(1) benchmarkMFcontract<ScalarA2Apolicies,GridA2Apolicies>(a2a_args, ntests, nthreads);


  printf("Finished\n"); fflush(stdout);
  
  return 0;
}
