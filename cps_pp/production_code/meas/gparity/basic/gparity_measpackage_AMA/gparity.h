#ifndef _GPARITY_H
#define _GPARITY_H

CPS_START_NAMESPACE

template<class T>
void decode_vml(const std::string &directory, const std::string &arg_name, T&into){
  std::string file = directory + std::string("/") + arg_name + std::string(".vml");
  if(!UniqueID()) printf("Decoding %s.vml from directory %s\n",arg_name.c_str(),directory.c_str() ); fflush(stdout); 
  if ( ! into.Decode(const_cast<char*>(file.c_str()), const_cast<char*>(arg_name.c_str()) ) ){
    std::string templ = directory + std::string("/") + arg_name + std::string(".templ");
    into.Encode(const_cast<char*>(templ.c_str()),const_cast<char*>(arg_name.c_str()));
    ERR.General("", "decode_vml", "Could not read %s\n",file.c_str());
  }
}


void decode_vml_all(DoArg &do_arg, BfmArg &bfm_arg, LancArg &lanc_arg_l, LancArg &lanc_arg_h, GparityAMAarg2 &ama_arg, const std::string &script_dir){
  const char* cname = "";
  const char *fname = "decode_vml_all()";
  decode_vml(script_dir,"do_arg",do_arg);
  decode_vml(script_dir,"bfm_arg",bfm_arg);
  decode_vml(script_dir,"lanc_arg_l",lanc_arg_l);
  decode_vml(script_dir,"lanc_arg_h",lanc_arg_h);
  decode_vml(script_dir,"ama_arg",ama_arg);
}

void check_bk_tsources(GparityAMAarg2 &ama_arg){
  const int Lt = GJP.Tnodes()*GJP.TnodeSites();
  const int lens[2] = { ama_arg.exact_solve_timeslices.exact_solve_timeslices_len, ama_arg.sloppy_solve_timeslices.sloppy_solve_timeslices_len };
  const int* vals[2] = { ama_arg.exact_solve_timeslices.exact_solve_timeslices_val, ama_arg.sloppy_solve_timeslices.sloppy_solve_timeslices_val };

  for(int p=0;p<2;p++){
    std::string status = p==0 ? "exact" : "sloppy";
    if(!UniqueID()) printf("Checking %s BK sources\n",status.c_str());

    if(lens[p] == 0){
      if(!UniqueID()) printf("Skipping %s BK as no sources\n",status.c_str());
      continue;
    }

    std::vector<int> tsep_meas_count(ama_arg.bk_tseps.bk_tseps_len, 0);

    for(int tt=0;tt<lens[p];tt++){
      int t0 = vals[p][tt];
      for(int tsi=0;tsi<ama_arg.bk_tseps.bk_tseps_len;tsi++){
	int tsep = ama_arg.bk_tseps.bk_tseps_val[tsi];
	int t1 = (t0 + tsep) % Lt;
	bool found = false;
	for(int tt2=0;tt2<lens[p];tt2++) if(vals[p][tt2]==t1){ found = true; break; }
	if(found) tsep_meas_count[tsi]++;

	if(!UniqueID()) printf("BK t0=%d tsep=%d, want t1=%d, status %d, count now %d\n",t0,tsep,t1,found,tsep_meas_count[tsi]);
      }
    }

    for(int tsi=0;tsi<ama_arg.bk_tseps.bk_tseps_len;tsi++){
      int tsep = ama_arg.bk_tseps.bk_tseps_val[tsi];
      if(!UniqueID()) printf("BK with tsep %d will have %d translations\n",tsep,tsep_meas_count[tsi]);
      if(tsep_meas_count[tsi] == 0){
	if(!UniqueID()){
	  printf("Error: BK time sep %d and %s props, no measurements will be performed due to lack of source timeslices\n",tsep,status.c_str());
	  std::cout.flush();
	}
	exit(-1);
      }
    }
  }
}

void init_fbfm(int *argc, char **argv[], const BfmArg &args)
{
  if(!UniqueID()) printf("Initializing Fbfm\n");
  // /*! IMPORTANT: BfmSolverType is not the same as the BfmSolver in the bfm package. BfmSolverType is defined in enum.x. Basically it adds a BFM_ prefix to the corresponding names in the BfmSolver enum. */
  cps_qdp_init(argc,argv);
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
    if(!UniqueID()) printf("G-parity directions: ");
    for(int d=0;d<3;d++)
      if(GJP.Bc(d) == BND_CND_GPARITY){ Fbfm::bfm_args[0].gparity_dir[d] = 1; printf("%d ",d); }
      else Fbfm::bfm_args[0].gparity_dir[d] = 0;
    for(int d=0;d<4;d++){
      Fbfm::bfm_args[0].nodes[d] = procs[d];
    }
    printf("\n");
  }else if(!UniqueID()) printf("Standard boundary conditions\n");

  // mobius_scale = b + c in Andrew's notation
  bfmarg::mobius_scale = args.mobius_scale;

  //Fbfm::current_arg_idx = 0;

  bfmarg::Threads(args.threads); 
  omp_set_num_threads(args.threads);

  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(args.verbose);
}


inline int toInt(const char* a){
  std::stringstream ss; ss << a; int o; ss >> o;
  return o;
}

inline bool contains_pctd(const std::string &c){
  return c.find("%d") != std::string::npos;
}


CPS_END_NAMESPACE


#endif
