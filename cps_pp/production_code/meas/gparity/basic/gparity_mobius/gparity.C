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

#include <chroma.h>
#include <omp.h>
#include <pthread.h>

#include <alg/bfm_arg.h>

USING_NAMESPACE_CPS

//NOTE: Peter's BFM Shamir domain wall propagators are a factor of (5-m5) smaller than the CPS propagators


void init_bfm(int *argc, char **argv[], const BfmArg &args)
{
  // /*! IMPORTANT: BfmSolverType is not the same as the BfmSolver in the bfm package. BfmSolverType is defined in enum.x. Basically it adds a BFM_ prefix to the corresponding names in the BfmSolver enum. */

  Chroma::initialize(argc, argv);
  multi1d<int> nrow(Nd);
  
  for(int i = 0; i< Nd; ++i)
    nrow[i] = GJP.Sites(i);
  
  Layout::setLattSize(nrow);
  Layout::create();

  if(args.solver == BFM_HmCayleyTanh){
    Fbfm::bfm_args[0].solver = HmCayleyTanh;
  }else if(args.solver == BFM_DWF){
    Fbfm::bfm_args[0].solver = DWF;
  }else ERR.General("","init_bfm","CPS solver enum correspondance to bfm enum has not been added for input solver type\n");

  Fbfm::bfm_args[0].precon_5d = args.precon_5d;

  // mixed-precision CG *based on environment variable*, *true by default*
  char* use_mixed_solver_env = getenv( "use_mixed_solver" );
  Fbfm::use_mixed_solver = true;
  if ( use_mixed_solver_env && strcmp( use_mixed_solver_env, "false" ) == 0 )
    Fbfm::use_mixed_solver = false;
  VRB.Result( "cps", "init_bfm", "Fbfm::use_mixed_solver: %d\n", Fbfm::use_mixed_solver );
  
  Fbfm::bfm_args[0].Ls = GJP.SnodeSites();
  Fbfm::bfm_args[0].M5 = GJP.DwfHeight();
  Fbfm::bfm_args[0].mass = args.mass;
  Fbfm::bfm_args[0].residual = args.residual;
  Fbfm::bfm_args[0].max_iter = args.max_iter;
  Fbfm::bfm_args[0].Csw = args.Csw;
  
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
  bfmarg::mobius_scale = args.mobius_scale;

  //Fbfm::current_arg_idx = 0;

  bfmarg::Threads(args.threads); 

  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(args.verbose);
}

int main(int argc,char *argv[])
{
  Start(&argc,&argv);
  
  if(argc < 5){
    ERR.General("","main()","Not enough arguments. Require DoArg, JobPropagatorArgs, BfmArg and GparityContractArg");
  }

  bool binary_write = false; //enable output in binary where implemented

  int i=5;
  while(i<argc){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-binary_write",25) == 0){
      if(!UniqueID()) printf("Enabled binary write\n");
      binary_write = true;
      i++;
    }else{
      if(!UniqueID()) printf("Unknown argument: %s\n",cmd);
      exit(-1);
    }
  }


  DoArg do_arg;
  if(!do_arg.Decode(argv[1],"do_arg")){
    ERR.General("","main()","Failed to decode %s\n",argv[1]); exit(-1);
  }
  JobPropagatorArgs prop_args;
  if(!prop_args.Decode(argv[2],"prop_arg")){
    ERR.General("","main()","Failed to decode %s\n",argv[2]); exit(-1);
  }
  BfmArg bfm_args;
  if(!bfm_args.Decode(argv[3],"bfm_arg")){
    ERR.General("","main()","Failed to decode %s\n",argv[3]); exit(-1);
  }
  GparityContractArg contract_args;
  if(!contract_args.Decode(argv[4],"contract_arg")){
    ERR.General("","main()","Failed to decode %s\n",argv[4]); exit(-1);
  }

  if(contract_args.conf_start >= contract_args.conf_lessthan || contract_args.conf_incr == 0){
    ERR.General("","main()","Invalid configuration args");
  }

  if(!UniqueID()){
    printf("contract_args contains %d measurements:\n",contract_args.meas.meas_len);
    for(int m=0;m<contract_args.meas.meas_len;m++) contract_args.meas.meas_val[m].print();
  }

  char *c = contract_args.config_fmt;

  if(UniqueID()==0) printf("Configuration format is '%s'\n",contract_args.config_fmt);
  bool found(false);
  while(*c!='\0'){
    if(*c=='\%' && *(c+1)=='d'){ found=true; break; }
    c++;
  }
  if(!found) ERR.General("","main()","GparityContractArg config format '%s' does not contain a %%d",contract_args.config_fmt);

  GJP.Initialize(do_arg);
  GJP.SetNthreads(bfm_args.threads);

  GJP.StartConfKind(START_CONF_MEM); //we will handle the gauge field read/write thankyou!

#if TARGET == BGQ
  LRG.setSerial();
#endif

  LRG.Initialize();

  init_bfm(&argc,&argv,bfm_args);

  PropManager::setup(prop_args);
  if(UniqueID()==0){
    printf("prop_args contains %d propagators\n", prop_args.props.props_len);
    prop_args.print();
  }
  
  GnoneFbfm lattice;
  CommonArg carg("label","filename");
  char load_config_file[1000];

  for(int conf=contract_args.conf_start; conf < contract_args.conf_lessthan; conf += contract_args.conf_incr){
    PropManager::startNewTraj();

    //Read/generate the gauge configuration 
    if(do_arg.start_conf_kind == START_CONF_FILE){
    
      if(sprintf(load_config_file,contract_args.config_fmt,conf) < 0){
	ERR.General("","main()","Congfiguration filename creation problem : %s | %s",load_config_file,contract_args.config_fmt);
      }
      //load the configuration
      ReadLatticeParallel readLat;
      if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
      readLat.read(lattice,load_config_file);
      if(UniqueID()==0) printf("Config read.\n");
    }else if(do_arg.start_conf_kind == START_CONF_ORD){
      if(!UniqueID()) printf("Using unit gauge links\n");
      lattice.SetGfieldOrd();
    }else if(do_arg.start_conf_kind == START_CONF_DISORD){
      if(!UniqueID()) printf("Using random gauge links\n");
      lattice.SetGfieldDisOrd();
      printf("Gauge checksum = %d\n", lattice.CheckSum());
    }else{
      ERR.General("","main()","Invalid do_arg.start_conf_kind\n");
    }
    lattice.BondCond(); //apply BCs and import gauge field from CPS container

    //Gauge fix lattice if required
    if(contract_args.fix_gauge.fix_gauge_kind != FIX_GAUGE_NONE){
      AlgFixGauge fix_gauge(lattice,&carg,&contract_args.fix_gauge);
      fix_gauge.run();
    }

    //Perform the inversions/contractions
    AlgGparityContract contract(lattice,carg,contract_args);
    if(binary_write) contract.enable_binary_write();
    contract.run(conf);

    //Free the gauge fixing matrices and reset for next config
    if(contract_args.fix_gauge.fix_gauge_kind != FIX_GAUGE_NONE) lattice.FixGaugeFree();
  }

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}

