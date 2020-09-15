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
#include<util/lattice/bfm_mixed_solver.h>

#include <util/lat_cont.h>

#include <chroma.h>
#include <omp.h>
#include <pthread.h>

#include <sys/types.h>
#include <dirent.h>
#include <sstream>
#include <fstream>
//-------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

const char *cname = "";

HmcArg hmc_arg;

ActionGaugeArg gauge_arg;
ActionRationalQuotientArg rat_quo_l_hsb1_arg;
ActionRationalQuotientArg rat_quo_l_hsb2_arg;
ActionRationalQuotientArg rat_quo_tm_arg;
ActionRationalQuotientArg rat_quo_h_arg;

// define all integrators; vml file arguments will not change
IntABArg ab1_arg;
IntABArg ab2_arg;
IntABArg ab3_arg;
IntABArg ab4_arg;

//BfmArg b_arg;
EvoArg evo_arg;
DoArg do_arg;
PbpArg pbp_arg;
NoArg no_arg;
//ApeSmearArg ape_arg;

void checkpoint(int traj);

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
    decode_vml(hmc_arg);
    decode_vml(evo_arg);
    decode_vml(gauge_arg);
    decode_vml(rat_quo_tm_arg);
    decode_vml(rat_quo_l_hsb1_arg);
    decode_vml(rat_quo_l_hsb2_arg);
    decode_vml(rat_quo_h_arg);
    decode_vml(ab1_arg);
    decode_vml(ab2_arg);
    decode_vml(ab3_arg);
    decode_vml(ab4_arg);
    decode_vml(pbp_arg);
    // decode_vml(ape_arg);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj);
void measure_plaq(CommonArg &common_arg);
//void measure_wline(CommonArg &common_arg);
//void measure_tc(CommonArg &common_arg, int cycle);
void measure_pbp(CommonArg &common_arg, int traj);
void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab);

void config_repro(const char* repro_config, const Float &tolerance);

void init_bfm(int *argc, char **argv[])
{
    cps_qdp_init(argc,argv);
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
    if(!UniqueID()) printf( "cps::init_bfm : Fbfm::use_mixed_solver: %d\n", Fbfm::use_mixed_solver );

    //Use mixed precision multi-shift
    char* mixed_prec_multishift_env = getenv( "use_mixedprec_multishift" );
    MultiShiftController.setEnvironmentMode(MultiShiftCGcontroller::MolecularDynamics, MultiShiftCGcontroller::SINGLE_PREC_RELIABLE_UPDATE_PLUS_OUTER_DEFECT_CORRECTION_LOOP);
    //MultiShiftController.setEnvironmentMode(MultiShiftCGcontroller::EnergyCalculation, MultiShiftCGcontroller::SINGLE_PREC_RELIABLE_UPDATE_PLUS_OUTER_DEFECT_CORRECTION_LOOP);
    //MultiShiftController.setEnvironmentMode(MultiShiftCGcontroller::Heatbath, MultiShiftCGcontroller::SINGLE_PREC_RELIABLE_UPDATE_PLUS_OUTER_DEFECT_CORRECTION_LOOP);

    MultiShiftController.setEnvironmentMode(MultiShiftCGcontroller::EnergyCalculation,MultiShiftCGcontroller::DOUBLE_PREC); 
    MultiShiftController.setEnvironmentMode(MultiShiftCGcontroller::Heatbath,MultiShiftCGcontroller::DOUBLE_PREC);

 
    MultiShiftController.setReliableUpdateFreq(200);
    MultiShiftController.setMaximumDefectCorrectionCycles(20);

    if ( mixed_prec_multishift_env && strcmp( mixed_prec_multishift_env, "false" ) == 0 ){
      MultiShiftController.setEnvironmentMode(MultiShiftCGcontroller::MolecularDynamics,MultiShiftCGcontroller::DOUBLE_PREC);
      MultiShiftController.setEnvironmentMode(MultiShiftCGcontroller::EnergyCalculation,MultiShiftCGcontroller::DOUBLE_PREC); 
      MultiShiftController.setEnvironmentMode(MultiShiftCGcontroller::Heatbath,MultiShiftCGcontroller::DOUBLE_PREC);
    } 

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

    bfmarg::Threads(GJP.Nthreads());
    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(1);
}



//By setting quo_tm_arg.bi_arg.fermion = F_CLASS_BFM_TYPE2 we enable using Fbfm to do the Twisted Mass fermions also
void init_bfm_wilsontm(){
  if(!UniqueID()) printf("Doing DSDR term with Fbfm also\n");

  Fbfm::bfm_args[1].solver = WilsonTM;
  Fbfm::bfm_args[1].precon_5d = 0;

  Fbfm::bfm_args[1].Ls = 1;
  Fbfm::bfm_args[1].M5 = 0;
  Fbfm::bfm_args[1].mass = 0.1; //These are overwriten with the correct values within the evolution code
  Fbfm::bfm_args[1].twistedmass = 0.1;
  Fbfm::bfm_args[1].residual = 1e-8;
  Fbfm::bfm_args[1].max_iter = 100000;
  Fbfm::bfm_args[1].Csw = 0.0;
  Fbfm::bfm_args[1].time_report_iter = -16;

  Fbfm::bfm_args[1].node_latt[0] = QDP::Layout::subgridLattSize()[0];
  Fbfm::bfm_args[1].node_latt[1] = QDP::Layout::subgridLattSize()[1];
  Fbfm::bfm_args[1].node_latt[2] = QDP::Layout::subgridLattSize()[2];
  Fbfm::bfm_args[1].node_latt[3] = QDP::Layout::subgridLattSize()[3];

  multi1d<int> procs = QDP::Layout::logicalSize();
  multi1d<int> ncoor = QDP::Layout::nodeCoord();

  Fbfm::bfm_args[1].local_comm[0] = procs[0] > 1 ? 0 : 1;
  Fbfm::bfm_args[1].local_comm[1] = procs[1] > 1 ? 0 : 1;
  Fbfm::bfm_args[1].local_comm[2] = procs[2] > 1 ? 0 : 1;
  Fbfm::bfm_args[1].local_comm[3] = procs[3] > 1 ? 0 : 1;

  for(int i=0;i<4;i++) Fbfm::bfm_args[1].ncoor[i] = ncoor[i];

  if(GJP.Gparity()){
    Fbfm::bfm_args[1].gparity = 1;
    printf("G-parity directions: ");
    for(int d=0;d<3;d++)
      if(GJP.Bc(d) == BND_CND_GPARITY){ Fbfm::bfm_args[1].gparity_dir[d] = 1; printf("%d ",d); }
      else Fbfm::bfm_args[1].gparity_dir[d] = 0;
    for(int d=0;d<4;d++){
      Fbfm::bfm_args[1].nodes[d] = procs[d];
    }
    printf("\n");
  }
}

//A string to int that ensures the string contains no starting or trailing characters
int str_to_int(const std::string &str, bool &fail){
  int idx;
  std::stringstream ss; ss << str; ss >> std::noskipws >> idx;
  fail = (ss.bad() || ss.fail() || !ss.eof());
  return idx;
}

bool copy_file(const std::string &to, const std::string &from){
  std::ifstream  src(from.c_str(), std::ios::binary);
  std::ofstream  dst(to.c_str(), std::ios::binary);

  dst << src.rdbuf();
  return !(src.fail() || src.bad() || dst.fail() || dst.bad());
}

//Check a directory exists and open it on node 0
DIR* check_open_dir(const std::string &dir, const std::string &descr = ""){
  DIR* ret;
  Float fail = 0.;
  if(!UniqueID()){
    ret = opendir(dir.c_str());
    fail = (ret == NULL) ? 1. : 0.;
  }
  glb_sum_five(&fail);
  if(fail > 0.)
    ERR.General("","check_open_dir","Could not open %s directory '%s'\n", descr.c_str(), dir.c_str());
  return ret;
}

//Returns non-zero float on all nodes if the condition 'fail_node0' is not satisfied on node 0
Float node0_fail_check(const bool fail_node0){
  Float fail = 0.;
  if(!UniqueID() && fail_node0)
    fail = 1.;
  glb_sum_five(&fail);
  return fail;
}


//On KEKSC we cannot execute a bash script to sort out the vml files at the start, hence we have to do it internally. This comprises the following stages:
//1) Identify most recently-generated configuration
//2) Copy the associated vml files to the scripts directory
//It assumes the standard setup with a root directory containing 'scripts', 'configurations', and 'work' subdirectories. argv[1] should be the 'scripts' directory
void setup_vml(int argc, char **argv){
  std::string scripts_dir = argv[1];
  std::string root_dir = scripts_dir + "/..";
  std::string conf_dir = root_dir + "/configurations";
  std::string work_dir = root_dir + "/work";
  
  bool fail;
  
  //Find the most recent config
  Float conf_dir_fail = 0.;
  
  DIR *conf_dir_p = check_open_dir(conf_dir, "configurations"); //opens on node 0
  
  int largest_idx = -1;
  if(!UniqueID()){
    printf("setup_vml: Contents of configurations directory:\n");
    
    struct dirent *dptr;
    while(NULL != (dptr = readdir(conf_dir_p)) ){
      std::string file = dptr->d_name;
      printf("File %s\n",dptr->d_name);
      
      size_t s = file.find("ckpoint_lat.");
      if(s == std::string::npos) continue;

      std::string idx_str = file.substr(s+12);
      int idx = str_to_int(idx_str, fail);
      if(fail){
	printf("Warning setup_vml : Could not convert index string '%s' to int\n",idx_str.c_str());
	continue;
      }
      
      printf("setup_vml : Found config '%s' with index %d\n",file.c_str(),idx);
      
      if(idx > largest_idx) largest_idx = idx;    
    }
  }

  if(node0_fail_check(largest_idx == -1) >0.)
    ERR.General("","setup_vml", "Could not find any configurations in directory %s\n",conf_dir.c_str());

  if(!UniqueID()){
    printf("setup_vml : Found most recent config with index %d\n",largest_idx);
    closedir(conf_dir_p);
  }
    
  //Find all files ending with .<idx> in the work directory and copy them to the scripts dir with ending .vml
  DIR *work_dir_p = check_open_dir(work_dir,"work");

  bool copy_fail = false;
  if(!UniqueID()){
    printf("setup_vml: Copying files ending '.%d' from work directory %s\n",largest_idx,work_dir.c_str());
    std::string find_end;  { std::ostringstream os; os << '.' << largest_idx; find_end = os.str(); }
    struct dirent *dptr;
    while(NULL != (dptr = readdir(work_dir_p)) ){
      std::string file = dptr->d_name;
      printf("File %s\n",dptr->d_name);
      
      size_t end_pos;
      if( (end_pos = file.find(find_end)) == std::string::npos)
	continue;

      printf("This file '%s' has the correct ending\n",dptr->d_name);
      
      std::string file_from = work_dir + '/' + file;
      std::string file_to = scripts_dir + '/' + file.substr(0,end_pos) + ".vml";
      printf("Copying '%s' to '%s'\n",file_from.c_str(),file_to.c_str());
    
      if(!copy_file(file_to,file_from)){
	copy_fail = true;
	break;
      }
    }
    closedir(work_dir_p);
  }
  
  if(node0_fail_check(copy_fail) > 0.) //this should also sync the nodes
    ERR.General("","setup_vml","File copy failed!\n");  
}


bool fexists(const std::string &filename){
  std::ifstream ifile(filename.c_str());
  return ifile;
}

//If no lock file exists in work directory, add a lock file. If it exists, exit the job.
//Assumes argv[1] is the scripts directory
//Returns a pair where the first entry indicates if the lock file was created and the second its filename
std::pair<bool, std::string> check_lock(int argc, char **argv){
  std::pair<bool, std::string> ret(false,"");
  bool do_lock = false;
  for(int i=1;i<argc;i++)
    if(std::string(argv[i]) == "-lock")
      do_lock = true;

  if(!do_lock) return ret;

  std::string scripts_dir = argv[1];
  std::string root_dir = scripts_dir + "/..";
  std::string work_dir = root_dir + "/work"; 

  std::string &lock_file = ret.second;
  lock_file = work_dir + "/lock";

  //Check if directory exists  
  bool dir_fail = false;
  if(!UniqueID()){
    DIR* work_dir_p = opendir(work_dir.c_str());
    if(work_dir_p == NULL)
      dir_fail = true;
    else
      closedir(work_dir_p);
  }

  if(node0_fail_check(dir_fail) > 0.)
    ERR.General("","check_lock","Could not open work directory '%s'\n", work_dir.c_str());

  bool lock_exists = false;  
  if(!UniqueID()){
    //Check if the lock file exists; if so, exit    
    if(fexists(lock_file))
      lock_exists = true;
  }
  if(node0_fail_check(lock_exists) > 0.)
    ERR.General("","check_lock","Lock file exists; some job is presently running evolution. Exiting\n");

  //Create the lock file
  bool lock_write_fail = false;
  if(!UniqueID()){
    std::ofstream of(lock_file.c_str());
    if(!of)
      lock_write_fail = true;
    else{    
      of << "Locked!\n";
      of.close();
      printf("LOCK: Applied lock '%s'\n",lock_file.c_str());
    }
  }  
  if(node0_fail_check(lock_write_fail) > 0.)
    ERR.General("","check_lock","Failed to create new lock file\n");
  
  ret.first = true;
}

void remove_lock(const std::pair<bool, std::string> &lock_info){
  if(lock_info.first){
    bool lock_remove_fail = false;
    if(!UniqueID()){
      if(remove(lock_info.second.c_str()) != 0)
	lock_remove_fail = true;
    }
    if(node0_fail_check(lock_remove_fail) > 0.)
      ERR.General("","remove_lock","Failed to remove lock file\n");
    
    if(!UniqueID()) printf("LOCK: Removed lock '%s'\n",lock_info.second.c_str());
  }
}
  

void setup(int argc, char **argv)
{
    const char *fname = "setup()";
    
    if(argc < 2) {
        ERR.General(cname, fname, "Must provide VML directory.\n");
    }
    
    //Check for -setup_vml option
    for(int i=1;i<argc;i++)
      if(std::string(argv[i]) == "-setup_vml")
	setup_vml(argc,argv);
    
    if(chdir( argv[1]) != 0) {
      ERR.General(cname, fname, "Changing directory to %s failed.\n", (*argv)[1]);
    }

    decode_vml_all();

    if(chdir(evo_arg.work_directory) != 0) {
        ERR.General(cname, fname, "Changing directory to %s failed.\n", evo_arg.work_directory);
    }
    VRB.Result(cname, fname, "Reading VML files successfully.\n");

#if TARGET == BGQ
  LRG.setSerial();
#endif

    GJP.Initialize(do_arg);
    LRG.Initialize();

    int threads = 32;
    if ( getenv ("BFM_NUM_THREADS") ) threads = atoi(getenv("BFM_NUM_THREADS")); 
    GJP.SetNthreads(threads);

    VRB.Result(cname, fname, "VRB.Level(%d)\n", do_arg.verbose_level);
    VRB.Level(do_arg.verbose_level);

    init_bfm(&argc, &argv);
    VRB.Result(cname, fname, "Mobius scale (2c+1) = %.10f\n", bfmarg::mobius_scale);

    //Use Fbfm to do the Twisted Mass Wilson part also
    if(rat_quo_tm_arg.bi_arg.fermion == F_CLASS_BFM_TYPE2) init_bfm_wilsontm();
    else if(rat_quo_tm_arg.bi_arg.fermion != F_CLASS_WILSON_TM) ERR.General(cname,fname,"rat_quo_tm_arg fermion type can only be either F_CLASS_WILSON_TM or F_CLASS_BFM_TYPE2");
}

int main(int argc, char *argv[])
{
    const char *fname = "main()";
    Start(&argc, &argv);
    
    std::pair<bool, std::string> lock_info = check_lock(argc, argv);
    
    setup(argc, argv);

    //CK: Option to run through 'gauge_unload_period' trajectories then compare the resulting gauge field
    //    with some other loaded field. Job will continue after comparison returns true
    char* repro_config;
    bool do_config_repro = false;
    Float repro_tolerance = 1e-09;

    int i=1;
    while(i<argc){
      char* cmd = argv[i];  
      if( strncmp(cmd,"-repro_config",15) == 0){
	if(i==argc-1){ printf("-repro_config requires an argument\n"); exit(-1); }
	repro_config = argv[i+1];
	if(!UniqueID()) printf("Running reproduction test against config %s at end of first unload period\n",repro_config);
	do_config_repro = true;
	i+=2;
      }else if( strncmp(cmd,"-repro_tolerance",20) == 0){
	if(i==argc-1){ printf("-repro_tolerance requires an argument\n"); exit(-1); }
	repro_tolerance = atof(argv[i+1]);
	if(repro_tolerance == 0.0) ERR.General(cname,fname,"Either you set repro_tolerance to 0.0 or the str->double conversion failed, either way that's an invalid input!");
	if(!UniqueID()) printf("Set reproduction tolerance to %e\n",repro_tolerance);
	i+=2;
      }else if( strncmp(cmd,"-bfm_cps_eval_norm",30) == 0 && rat_quo_tm_arg.bi_arg.fermion == F_CLASS_BFM_TYPE2){
	//Modify the eigenvalues in rat_quo_tm_arg to take into account the different field normalizations. With this option active, the results using Fbfm twisted mass wilson fermions 
	//will agree with the native CPS results. If not active they will disagree slightly as the precision of the rational approximations will differ
	if(!UniqueID()) printf("Normalizing eigenvalues such that Fbfm and CPS twisted wilson fermions agree\n");
       
	int ndet = rat_quo_tm_arg.bosons.bosons_len;
	for(int j=0;j<ndet;j++){
	  //Bosons
	  Float mass = rat_quo_tm_arg.bsn_mass.bsn_mass_val[j];
	  Float epsilon = rat_quo_tm_arg.bsn_mass_epsilon.bsn_mass_epsilon_val[j];

	  Float kappa = 1.0/2.0/sqrt( (mass + 4.0)*(mass + 4.0) + epsilon * epsilon );
	  Float fbfm_factor = 4*kappa*kappa*(4+mass);
	  
	  ApproxDescr *approx = &rat_quo_tm_arg.bosons.bosons_val[j].md_approx;
	  approx->lambda_low /= fbfm_factor;
	  approx->lambda_high /= fbfm_factor;
	  
	  approx = &rat_quo_tm_arg.bosons.bosons_val[j].mc_approx;
	  approx->lambda_low /= fbfm_factor;
	  approx->lambda_high /= fbfm_factor;

	  //Fermions
	  mass = rat_quo_tm_arg.frm_mass.frm_mass_val[j];
	  epsilon = rat_quo_tm_arg.frm_mass_epsilon.frm_mass_epsilon_val[j];

	  kappa = 1.0/2.0/sqrt( (mass + 4.0)*(mass + 4.0) + epsilon * epsilon );
	  fbfm_factor = 4*kappa*kappa*(4+mass);

	  approx = &rat_quo_tm_arg.fermions.fermions_val[j].md_approx;
	  approx->lambda_low /= fbfm_factor;
	  approx->lambda_high /= fbfm_factor;

	  approx = &rat_quo_tm_arg.fermions.fermions_val[j].mc_approx;
	  approx->lambda_low /= fbfm_factor;
	  approx->lambda_high /= fbfm_factor;
	}
	i++;
      }else{
	i++;
      }
    }

    //VRB.ElapsedTime("CPS", fname);

    // OpenMP test
    VRB.Result( "CPS", "main", "omp_get_num_threads[1] -> %d", omp_get_num_threads() );
    
    #pragma omp parallel
    {
      if ( UniqueID() == 0 && omp_get_thread_num() == 0 ) {
        VRB.Result( "CPS", "main", "omp_get_num_threads[2] -> %d", omp_get_num_threads() );
      }
    }

    VRB.Result( "CPS", "main", "omp_get_num_threads[3] -> %d", omp_get_num_threads() );
    
    //////////////////////////////////////////////////////////////////////
    // creating actions and integrators
    AlgMomentum mom;
    AlgActionGauge gauge(mom, gauge_arg);
    AlgActionRationalQuotient rat_quo_tm(mom, rat_quo_tm_arg);
    AlgActionRationalQuotient rat_quo_l_hsb1(mom, rat_quo_l_hsb1_arg); // light->intermediate Hasenbusch mass
    AlgActionRationalQuotient rat_quo_l_hsb2(mom, rat_quo_l_hsb2_arg); // intermediate -> PV Hasenbusch mass, on a finer timestep
    AlgActionRationalQuotient rat_quo_h(mom, rat_quo_h_arg);

    IntABArg sum_arg;
    sum_arg.A_steps = 1;
    sum_arg.B_steps = 1;
    sum_arg.level = EMBEDDED_INTEGRATOR;
    AlgIntSum sum(rat_quo_l_hsb2, rat_quo_h, sum_arg);

    //!< Construct numerical integrators
    AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge,       ab1_arg);
    AlgIntAB &ab2 = AlgIntAB::Create(ab1, rat_quo_tm,  ab2_arg);
    AlgIntAB &ab3 = AlgIntAB::Create(ab2, sum, ab3_arg);
    AlgIntAB &ab4 = AlgIntAB::Create(ab3, rat_quo_l_hsb1, ab4_arg);

    //////////////////////////////////////////////////////////////////////

    int traj = evo_arg.traj_start;
    for(int conf = 0; conf< evo_arg.gauge_configurations; ++conf) {
        CommonArg common_arg_plaq;
        CommonArg common_arg_pbp;
        CommonArg common_arg_hmc;

        truncate_it(&common_arg_plaq , evo_arg.plaquette_stem      , traj);
        truncate_it(&common_arg_pbp  , evo_arg.pbp_stem            , traj);
        truncate_it(&common_arg_hmc  , evo_arg.evo_stem            , traj);

        // Inner trajectory loop
        for(int i = 0; i < evo_arg.gauge_unload_period; ++i, ++traj) {
            measure_plaq(common_arg_plaq);
            measure_pbp(common_arg_pbp, traj);

            //VRB.ElapsedTime("CPS", "main[2]");
            VRB.Result( "CPS", "main", "omp_get_num_threads[4] -> %d", omp_get_num_threads() );
            run_hmc(common_arg_hmc, traj, ab4);
        }//End of inter-cfg sweep

	if(do_config_repro){
	  config_repro(repro_config, repro_tolerance);
	  do_config_repro = false;
	}
        checkpoint(traj);

    } //End config loop

    AlgIntAB::Destroy(ab4);
    AlgIntAB::Destroy(ab3);
    AlgIntAB::Destroy(ab2);
    AlgIntAB::Destroy(ab1);

    End();

    remove_lock(lock_info);
    VRB.Result(cname, fname, "Program ended normally.\n");
    
    return 0;
}

#undef encode_vml
#define encode_vml(arg_name, traj) do{                                  \
        char vml_file[256];                                             \
        sprintf(vml_file, #arg_name".%d", traj);                        \
        if( !arg_name.Encode(vml_file, #arg_name) ){                    \
            ERR.General(cname, fname, #arg_name " encoding failed.\n"); \
        }                                                               \
    }while(0)

void config_repro(const char* repro_config, const Float &tolerance){
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  LatticeContainer cfg_store;
  cfg_store.Get(lat);
  
  {
    ReadLatticeParallel readLat;
    if(!UniqueID()) printf("Reading: configuration %s for reproduction test\n",repro_config);
    readLat.read(lat,repro_config);
  }
  
  Float* compare_to = (Float*)lat.GaugeField();
  Float* generated = (Float*)cfg_store.GaugeField();
  
  Float fail(0);

  long gauge_size = 18*4*GJP.VolNodeSites();

  for(int i=0;i<gauge_size;i++){
    int rem = i;
    int midx = rem % 18; rem/=18;
    int mu = rem % 4; rem/=4;
    int x[4];
    for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

    if(fabs(generated[i] - compare_to[i]) > tolerance){
      printf("Reproduce comparison fail: midx %d mu %d (%d %d %d %d). Generated: %.10f Loaded: %.10f\n",midx,mu,x[0],x[1],x[2],x[3],generated[i],compare_to[i]);
      fail=1;
    }
  }
  glb_sum(&fail);

  if(fail > 0 && !UniqueID()){ printf("Failed reproduction test\n"); exit(-1); }
  else if(!UniqueID()) printf("Passed reproduction test\n");

  cfg_store.Set(lat);
  LatticeFactory::Destroy();
}



void checkpoint(int traj)
{
    const char *fname="checkpoint()";

    char lat_file[256];
    char rng_file[256];

    Float time = -dclock();

    // Save this config to disk
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

    sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
    QioArg wt_arg(lat_file,0.001);

    wt_arg.ConcurIONumber=evo_arg.io_concurrency;
    WriteLatticeParallel wl;
#if TARGET == BGQ
  wl.setSerial();
#endif
    wl.setHeader(evo_arg.ensemble_id,evo_arg.ensemble_label,traj);
    wl.write(lat,wt_arg);

    if(!wl.good())
        ERR.General(cname,fname,"Failed write lattice %s",lat_file);

    LatticeFactory::Destroy();

    // Save the RNG's
    sprintf(rng_file,"%s.%d",evo_arg.rng_file_stem,traj);
    if ( !LRG.Write(rng_file) )
        ERR.General(cname,fname,"Failed write RNG file %s",rng_file);

    // Update the parameter files for restart
    do_arg.start_seed_filename = rng_file;
    do_arg.start_seed_kind = START_SEED_FILE;
    do_arg.start_conf_filename = lat_file;
    do_arg.start_conf_kind = START_CONF_FILE;
    evo_arg.traj_start     = traj;

    //encode_vml(b_arg, traj);
    encode_vml(hmc_arg, traj);
    encode_vml(gauge_arg, traj);
    encode_vml(rat_quo_tm_arg, traj);
    encode_vml(rat_quo_l_hsb1_arg, traj);
    encode_vml(rat_quo_l_hsb2_arg, traj);
    encode_vml(rat_quo_h_arg, traj);
    encode_vml(ab1_arg, traj);
    encode_vml(ab2_arg, traj);
    encode_vml(ab3_arg, traj);
    encode_vml(ab4_arg, traj);
    encode_vml(pbp_arg, traj);
    encode_vml(do_arg, traj);
    encode_vml(evo_arg, traj);

    time += dclock();
    print_flops("","checkpoint()",0,time);
}

void truncate_it(CommonArg *common_arg, const char stem[], int traj)
{
    char fnbuf[1024];
    sprintf(fnbuf, "%s.%d", stem, traj);
    FILE *truncate_it = Fopen(fnbuf, "w");
    Fclose(truncate_it);
    common_arg->set_filename(fnbuf);
}

void measure_plaq(CommonArg &common_arg)
{
    const char *fname = "measure_plaq()";

    Float dtime = -dclock();

    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
    AlgPlaq plaq(lat, &common_arg, &no_arg);
    plaq.run();
    LatticeFactory::Destroy();

    dtime += dclock();
    print_flops("AlgPlaq", "run()", 0, dtime);	
}

void measure_pbp(CommonArg &common_arg, int traj)
{
    const char *fname = "measure_pbp()";

    // fix pbp_arg
    pbp_arg.src_u_s = 0;
    pbp_arg.src_l_s = GJP.SnodeSites() * GJP.Snodes() - 1;
    pbp_arg.snk_u_s = GJP.SnodeSites() * GJP.Snodes() - 1;
    pbp_arg.snk_l_s = 0;

    const int g_int = evo_arg.gauge_unload_period;
    if (traj % g_int == 0 && evo_arg.measure_pbp) {
        Float dtime = -dclock();

        LRGState rng_state;
        rng_state.GetStates();

        Lattice &lat = LatticeFactory::Create(F_CLASS_BFM, G_CLASS_NONE);
        VRB.Result( "cps", "measure_pbp", "LatticeFactory::Create(F_CLASS_BFM, G_CLASS_NONE)" );

        AlgPbp pbp(lat, &common_arg, &pbp_arg);

        for(int pbp_counter = 0; pbp_counter < evo_arg.measure_pbp; pbp_counter++) {
            pbp.run();
        }
        LatticeFactory::Destroy();
        rng_state.SetStates();

        dtime += dclock();
        print_flops("AlgPbp", "run()", 0, dtime);	
    }
}

void run_hmc(CommonArg &common_arg, int traj, AlgIntAB &int_ab)
{
    const char *fname = "run_hmc()";

    Float dtime = -dclock();

    if ( (evo_arg.reproduce_interval > 0) &&
         (traj % evo_arg.reproduce_interval) == 0 ) {
        VRB.Result(cname,fname,"Running traj %d with reproduction\n",traj);
        hmc_arg.reproduce = REPRODUCE_YES;
    } else {
        VRB.Result(cname,fname,"Running traj %d without reproduction\n",traj);
        hmc_arg.reproduce = REPRODUCE_NO;	
    }
    
    AlgHmc hmc(int_ab, common_arg, hmc_arg);
    hmc.run();

    dtime += dclock();
    print_flops("AlgHmc", "run()", 0, dtime);
}



// void measure_wline(CommonArg &common_arg)
// {
//     const char *fname = "measure_wline()";

//     Float dtime = -dclock();

//     Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
//     AlgWline wline(lat, &common_arg, &no_arg);
//     wline.run();
//     LatticeFactory::Destroy();

//     dtime += dclock();
//     print_flops("AlgWline", "run()", 0, dtime);	
// }

// void measure_tc(CommonArg &common_arg, int cycle)
// {
//     const char *fname = "measure_tc()";

//     Float dtime = -dclock();

//     // calculate the topological charge. Need to copy the lattice since
//     // we need to smear it first. Use Chulwoo's "lattice container"
//     Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
//     LatticeContainer lat_cont;
//     // copy lattice
//     lat_cont.Get(lat);

//     //----- mess up lattice -------------------------
//     AlgApeSmear ape(lat, &common_arg, &ape_arg);
//     AlgTcharge  tcharge(lat, &common_arg);
//     for (int i = 0; i < cycle; ++i) {
//         VRB.Result(cname,fname,"%i\n",i);
//         VRB.Result(cname,fname,"   running tcharge\n"); tcharge.run();
//         VRB.Result(cname,fname,"   running ape\n"); ape.run();
//         VRB.Result(cname,fname,"   running ape\n"); ape.run();
//         VRB.Result(cname,fname,"   running ape\n"); ape.run();
//     }
//     tcharge.run();
//     // restore the lattice
//     lat_cont.Set(lat);
//     LatticeFactory::Destroy();

//     dtime += dclock();
//     print_flops("AlgTcharge", "run()", 0, dtime);
// }
