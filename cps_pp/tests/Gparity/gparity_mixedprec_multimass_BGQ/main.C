#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/lattice/fbfm.h>
#include<util/random.h>
#include<util/time_cps.h>

#include<alg/alg_hmc.h>
#include<alg/bfm_arg.h>
#include<alg/common_arg.h>
#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>

#include<alg/alg_int.h>
#include<alg/int_arg.h>
#include<alg/alg_wline.h>
#include<alg/no_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_plaq.h>
#include<alg/alg_pbp.h>
#include<alg/alg_tcharge.h>
#include<alg/alg_smear.h>
#include<alg/ape_smear_arg.h>
#include<alg/pbp_arg.h>
#include<alg/alg_remez.h>
#include<alg/alg_w_spect.h>
#include<alg/array_arg.h>
#include<alg/alg_fix_gauge.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>

#include <util/lat_cont.h>

#include <chroma.h>
#include <omp.h>
#include <pthread.h>
#include <sys/time.h>
//-------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

const char *cname = "";

ActionRationalQuotientArg rat_quo_arg;
DoArg do_arg;

#define decode_vml(arg_name) \
  if(!UniqueID()) printf("Decoding %s.vml\n",#arg_name); fflush(stdout); \
  do{									\
        if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
            ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
    } while(0)

void decode_vml_all(void)
{
    char *fname = "decode_vml_all()";

    //decode_vml(b_arg);
    decode_vml(do_arg);
    decode_vml(rat_quo_arg);
}

void init_bfm(int *argc, char **argv[])
{
    Chroma::initialize(argc, argv);
    multi1d<int> nrow(Nd);

    for(int i = 0; i< Nd; ++i)
        nrow[i] = GJP.Sites(i);

    Layout::setLattSize(nrow);
    Layout::create();

    //Set for Mobius fermions
    Fbfm::bfm_args[0].solver = HmCayleyTanh;
    Fbfm::bfm_args[0].precon_5d = 0;

    // mixed-precision CG *based on environment variable*, *true by default*
    char* use_mixed_solver_env = getenv( "use_mixed_solver" );
    Fbfm::use_mixed_solver = true;
    if ( use_mixed_solver_env && strcmp( use_mixed_solver_env, "false" ) == 0 )
      Fbfm::use_mixed_solver = false;
    VRB.Result( "cps", "init_bfm", "Fbfm::use_mixed_solver: %d\n", Fbfm::use_mixed_solver );

    Fbfm::bfm_args[0].Ls = GJP.SnodeSites();
    Fbfm::bfm_args[0].M5 = GJP.DwfHeight();
    Fbfm::bfm_args[0].mass = 0.1;
    Fbfm::bfm_args[0].residual = 1e-8;
    Fbfm::bfm_args[0].max_iter = 100000;
    Fbfm::bfm_args[0].Csw = 0.0;

    Fbfm::bfm_args[0].node_latt[0] = QDP::Layout::subgridLattSize()[0];
    Fbfm::bfm_args[0].node_latt[1] = QDP::Layout::subgridLattSize()[1];
    Fbfm::bfm_args[0].node_latt[2] = QDP::Layout::subgridLattSize()[2];
    Fbfm::bfm_args[0].node_latt[3] = QDP::Layout::subgridLattSize()[3];

    multi1d<int> procs = QDP::Layout::logicalSize();
    multi1d<int> ncoor = QDP::Layout::nodeCoord();

    Fbfm::bfm_args[0].local_comm[0] = procs[0] > 1 ? 0 : 1;
    Fbfm::bfm_args[0].local_comm[1] = procs[1] > 1 ? 0 : 1;
    Fbfm::bfm_args[0].local_comm[2] = procs[2] > 1 ? 0 : 1;
    Fbfm::bfm_args[0].local_comm[3] = procs[3] > 1 ? 0 : 1;

    for(int i=0;i<4;i++) Fbfm::bfm_args[0].ncoor[i] = ncoor[i];

    if(GJP.Gparity()){
      Fbfm::bfm_args[0].gparity = 1;
      printf("G-parity directions: ");
      for(int d=0;d<3;d++)
	if(GJP.Bc(d) == BND_CND_GPARITY){ Fbfm::bfm_args[0].gparity_dir[d] = 1; printf("%d ",d); }
	else Fbfm::bfm_args[0].gparity_dir[d] = 0;
      for(int d=0;d<4;d++){
	Fbfm::bfm_args[0].nodes[d] = procs[d];
      }
      printf("\n");
    }

    // mobius_scale = b + c in Andrew's notation
    bfmarg::mobius_scale = 32.0/12.0; 

    Fbfm::current_arg_idx = 0;

#if TARGET == BGQ  
    bfmarg::Threads(64); //32
#else
    bfmarg::Threads(1);
#endif
    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(1);
}

void setup(int *argc, char ***argv)
{
    const char *fname = "setup()";

    Start(argc, argv);

    if(*argc < 2) {
        ERR.General(cname, fname, "Must provide VML directory.\n");
    }

    if(chdir( (*argv)[1]) != 0) {
      ERR.General(cname, fname, "Changing directory to %s failed.\n", (*argv)[1]);
    }

    decode_vml_all();
 
    VRB.Result(cname, fname, "Read VML files successfully.\n");

#if TARGET == BGQ
  LRG.setSerial();
#endif

    GJP.Initialize(do_arg);
    LRG.Initialize();

    VRB.Result(cname, fname, "VRB.Level(%d)\n", do_arg.verbose_level);
    VRB.Level(do_arg.verbose_level);

    init_bfm(argc, argv);
    VRB.Result(cname, fname, "Mobius scale (2c+1) = %.10f\n", bfmarg::mobius_scale);
}

int main(int argc, char *argv[])
{
    const char *fname = "main()";

    setup(&argc, &argv);

    if(argc<3) ERR.General("","main","Not enough arguments. Usage: <executable> <directory> <mode> <options>\n");
      
    int mode;
    { std::stringstream ss; ss << argv[2]; ss>>mode; }

    if(mode >= MultiShiftCGcontroller::NMultiShiftCGMode ) ERR.General("","main","Unknown multi-shift mode\n");
    MultiShiftController.setEnvironmentMode(MultiShiftCGcontroller::Generic, (MultiShiftCGcontroller::Mode)mode );
    MultiShiftController.setEnvironment(MultiShiftCGcontroller::Generic);

    int i=3;
    while(i<argc){
      char* cmd = argv[i];  
      if( strncmp(cmd,"-min_fp_resid",15) == 0){
	std::stringstream ss;
	ss << argv[i+1];
	double into;
	ss >> into;
	MultiShiftController.setMinimumSinglePrecResidual(into);
	if(UniqueID()) printf("Set minimum floating point residual for mixed-prec multi-mass shift to %e\n",into);
	i+=2;
      }else if( strncmp(cmd,"-rel_up_freq",15) == 0){
	std::stringstream ss;
	ss << argv[i+1];
	int into;
	ss >> into;
	MultiShiftController.setReliableUpdateFreq(into);
	if(UniqueID()) printf("Set reliable update frequency to %d\n",into);
	i+=2;
      }else ERR.General("","main","Unknown argument %s\n",cmd);
    }

    VRB.Result( "CPS", "main", "omp_get_num_threads[1] -> %d", omp_get_num_threads() );
    
    #pragma omp parallel
    {
      if ( UniqueID() == 0 && omp_get_thread_num() == 0 ) {
        VRB.Result( "CPS", "main", "omp_get_num_threads[2] -> %d", omp_get_num_threads() );
      }
    }

    VRB.Result( "CPS", "main", "omp_get_num_threads[3] -> %d", omp_get_num_threads() );
    
    AlgMomentum mom;
    AlgActionRationalQuotient rat_quo(mom, rat_quo_arg);
    
    struct timeval start,stop,diff;
    gettimeofday(&start,NULL);
    rat_quo.heatbath();
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff);      

    if(!UniqueID()) printf("Time %d.%6.6d s\n",diff.tv_sec,diff.tv_usec);
    VRB.Result(cname, fname, "Program ended normally.\n");
    
    return 0;
}
