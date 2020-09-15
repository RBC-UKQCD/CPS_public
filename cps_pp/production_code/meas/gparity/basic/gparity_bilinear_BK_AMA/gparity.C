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

#include <chroma.h>
#include <omp.h>
#include <pthread.h>

#include <alg/bfm_arg.h>


//CK 2015
//Compute all point and wall-sink bilinears as well as B_K using wall source propagators
//Because the exact solves are only performed on certain specified timeslices, this limits the K->Kbar separations that can be used for B_K to the available separations between those
//We must therefore specify at least 2 exact solve timeslices. If there are more than 1 pair with the same time separation, the B_K result will be averaged over those.

USING_NAMESPACE_CPS

template<class T>
void decode_vml(const std::string &directory, const std::string &arg_name, T&into){
  std::string file = directory + std::string("/") + arg_name + std::string(".vml");
  if(!UniqueID()) printf("Decoding %s.vml from directory %s\n",arg_name.c_str(),directory.c_str() ); fflush(stdout); 
  if ( ! into.Decode(const_cast<char*>(file.c_str()), const_cast<char*>(arg_name.c_str()) ) )
    ERR.General("", "decode_vml", "Could not read %s\n",file.c_str());
}


BfmArg bfm_arg; //note, only the solver and mobius scale are used here
DoArg do_arg;
LancArg lanc_arg_h;
LancArg lanc_arg_l;
GparityAMAbilBKarg ama_arg;

//In the course of execution, each propagator is computed on every timeslice using a sloppy solve, and on select timeslices with an exact solve.
//NOTE: The solver and mobius scale specified in the LanczosContainer in prop_arg will be ignored; the Lanczos will use the bfm instance in fbfm

void decode_vml_all(const std::string &script_dir){
  const char* cname = "";
  const char *fname = "decode_vml_all()";
  decode_vml(script_dir,"do_arg",do_arg);
  decode_vml(script_dir,"bfm_arg",bfm_arg);
  decode_vml(script_dir,"lanc_arg_h",lanc_arg_h);
  decode_vml(script_dir,"lanc_arg_l",lanc_arg_l);
  decode_vml(script_dir,"ama_arg",ama_arg);
}

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

CGAttrArg* addCGattr(PropagatorContainer &prop){
  AttributeContainer tmp; tmp.type = CG_ATTR;
  tmp.AttributeContainer_u.cg_attr.max_num_iter = 10000;
  prop.add(tmp);
  CGAttrArg* out;
  if(!prop.getAttr(out)) ERR.General("","addCGattr","Failed");
  return out;
}



//void computePropagator(const std::string &tag, const int &tsrc, const int &mass, const double &resid, const bool &gauge_fix, Lattice &latt){


int main(int argc,char *argv[])
{
  Start(&argc,&argv);
  if(argc<2) ERR.General("","main()","Arguments: Require a directory containing the vml files");
  decode_vml_all(argv[argc-1]);
  if(ama_arg.conf_start >= ama_arg.conf_lessthan || ama_arg.conf_incr == 0) ERR.General("","main()","Invalid configuration args");

  bool lanczos_tune = false;
  bool dbl_latt_storemode = false;
  Fbfm::use_mixed_solver = false;
  {
    int i = 1;
    while(i<argc-1){
      if( std::string(argv[i]) == "-lanczos_tune" ){
	lanczos_tune = true;
	i++;
      }else if( std::string(argv[i]) == "-load_dbl_latt" ){
	if(!UniqueID()) printf("Loading double latt\n");
	dbl_latt_storemode = true;
	i++;
      }else if( std::string(argv[i]) == "-use_mixed_solver" ){
	Fbfm::use_mixed_solver = true;
	i++;
      }else{
	ERR.General("","main","Unknown argument: %s",argv[i]);
      }
    }
  }

  const char *c = ama_arg.config_fmt;
  if(UniqueID()==0) printf("Configuration format is '%s'\n",ama_arg.config_fmt);
  bool found(false);
  while(*c!='\0'){
    if(*c=='\%' && *(c+1)=='d'){ found=true; break; }
    c++;
  }
  if(!found) ERR.General("","main()","GparityAMAarg config format '%s' does not contain a %%d",ama_arg.config_fmt);

  GJP.Initialize(do_arg);
  GJP.StartConfKind(START_CONF_MEM); //we will handle the gauge field read/write thankyou!

#if TARGET == BGQ
  LRG.setSerial();
#endif

  LRG.Initialize();

  init_bfm(&argc,&argv,bfm_arg);

  if(UniqueID()==0)
    if(Fbfm::use_mixed_solver) printf("Using Fbfm mixed precision solver\n");
    else printf("Using Fbfm double precision solver\n");

  GnoneFbfm lattice;
  CommonArg carg("label","filename");
  char load_config_file[1000];

  //Setup the light and strange lanczos
  LanczosContainerArg lc_arg_l;
  lc_arg_l.tag = strdup("lanc_l");
  lc_arg_l.tbc = GJP.Bc(3);
  lc_arg_l.lanc_arg.deep_copy(lanc_arg_l);

  LanczosContainer lanc_l = PropManager::addLanczos(lc_arg_l);

  LanczosContainerArg lc_arg_h;
  lc_arg_h.tag = strdup("lanc_h");
  lc_arg_h.tbc = GJP.Bc(3);
  lc_arg_h.lanc_arg.deep_copy(lanc_arg_h);  

  LanczosContainer lanc_h = PropManager::addLanczos(lc_arg_h);

  //Attach the bfm instance from fbfm to the lanczos solver so it doesn't have to maintain separate bfm instance (also propagates mobius scale and solver type, etc)
  //(Mass is set automatically by Lanczos class to that specified in the LancArg
  lanc_l.set_bfm(&lattice.bd);
  lanc_h.set_bfm(&lattice.bd);

  //Setup the sloppy propagator containers
  for(int t=0;t<GJP.Tnodes()*GJP.TnodeSites();t++){
    std::ostringstream tag_l;
    tag_l << "sloppy_l_t" << t;
    setupPropagator(tag_l.c_str(),t, ama_arg.ml, ama_arg.sloppy_precision, ama_arg.fix_gauge.fix_gauge_kind != FIX_GAUGE_NONE);
    
    std::ostringstream tag_h;
    tag_h << "sloppy_h_t" << t;
    setupPropagator(tag_h.c_str(),t, ama_arg.mh, ama_arg.sloppy_precision, ama_arg.fix_gauge.fix_gauge_kind != FIX_GAUGE_NONE);
  }
    
  //Setup the exact propagators
  for(int tt=0;tt<ama_arg.exact_solve_timeslices.exact_solve_timeslices_len;tt++){
    int t = ama_arg.exact_solve_timeslices.exact_solve_timeslices_val[tt];

    std::ostringstream tag_l;
    tag_l << "exact_l_t" << t;
    setupPropagator(tag_l.c_str(),t, ama_arg.ml, ama_arg.exact_precision, ama_arg.fix_gauge.fix_gauge_kind != FIX_GAUGE_NONE);
    
    std::ostringstream tag_h;
    tag_h << "exact_h_t" << t;
    setupPropagator(tag_h.c_str(),t, ama_arg.mh, ama_arg.exact_precision, ama_arg.fix_gauge.fix_gauge_kind != FIX_GAUGE_NONE);
  }

  for(int conf=ama_arg.conf_start; conf < ama_arg.conf_lessthan; conf += ama_arg.conf_incr){
    PropManager::startNewTraj();

    //Read/generate the gauge configuration 
    if(do_arg.start_conf_kind == START_CONF_FILE){
    
      if(sprintf(load_config_file,ama_arg.config_fmt,conf) < 0){
	ERR.General("","main()","Congfiguration filename creation problem : %s | %s",load_config_file,ama_arg.config_fmt);
      }
      //load the configuration
      ReadLatticeParallel readLat;
      if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
      if(dbl_latt_storemode){
	if(!UniqueID()) printf("Disabling U* field reconstruction\n");
	readLat.disableGparityReconstructUstarField();
      }
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
    lattice.ImportGauge(); //make sure the bfm lattice is up-to-date 

    //Gauge fix lattice if required
    if(ama_arg.fix_gauge.fix_gauge_kind != FIX_GAUGE_NONE){
      AlgFixGauge fix_gauge(lattice,&carg,&ama_arg.fix_gauge);
      fix_gauge.run();
    }

    //Compute the eigenvectors
    LanczosContainer *lancs[] = { &lanc_l, &lanc_h };
    for(int l=0;l<2;l++){
      LanczosContainer& lanc = *lancs[l];
      lanc.calcEig(lattice);

      if(Fbfm::use_mixed_solver){
	//Convert eigenvectors to single precision
	lanc.setPrecision(1);
      }
    }

    if(lanczos_tune) return 0;

    //Compute the AMA point-sink bilinears
    if(GJP.Gparity()) do_bilinears<SpinColorFlavorMatrix>(lattice,conf);
    else do_bilinears<WilsonMatrix>(lattice,conf);

    //Free the gauge fixing matrices and reset for next config
    if(ama_arg.fix_gauge.fix_gauge_kind != FIX_GAUGE_NONE) lattice.FixGaugeFree();
  }

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}

