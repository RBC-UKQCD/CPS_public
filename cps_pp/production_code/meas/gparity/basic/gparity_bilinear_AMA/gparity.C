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

USING_NAMESPACE_CPS

// #define decode_vml(arg_name) \
  // if(!UniqueID()) printf("Decoding %s.vml\n",#arg_name); fflush(stdout); \
  // do{									\
  //       if ( ! arg_name.Decode(#arg_name".vml", #arg_name) )            \
  //           ERR.General(cname, fname, "Bad " #arg_name ".vml.\n");      \
  //   } while(0)


template<class T>
void decode_vml(const std::string &directory, const std::string &arg_name, T&into){
  std::string file = directory + std::string("/") + arg_name + std::string(".vml");
  if(!UniqueID()) printf("Decoding %s.vml from directory %s\n",arg_name.c_str(),directory.c_str() ); fflush(stdout); 
  if ( ! into.Decode(const_cast<char*>(file.c_str()), const_cast<char*>(arg_name.c_str()) ) )
    ERR.General("", "decode_vml", "Could not read %s\n",file.c_str());
}


BfmArg bfm_arg; //note, only the solver and mobius scale are used here
DoArg do_arg;
JobPropagatorArgs prop_arg;
GparityAMAarg ama_arg;

//In the course of execution, each propagator is computed on every timeslice using a sloppy solve, and on select timeslices with an exact solve.
//In the arguments, you should not specify all of these times separately; just define one propagator, and the time coordinate specified will be over-written
//by the code

//NOTE: Each propagator should have a DeflatedCGAttrArg defining which Lanczos instance to use. If you do not specify this, the job will be very slow!
//NOTE2: Any residuals specified will be ignored/over-written by the AMA code
//NOTE3: The solver and mobius scale specified in the LanczosContainer in prop_arg will be ignored; the Lanczos will use the bfm instance in fbfm

void decode_vml_all(const std::string &script_dir){
  const char* cname = "";
  const char *fname = "decode_vml_all()";
  decode_vml(script_dir,"do_arg",do_arg);
  decode_vml(script_dir,"bfm_arg",bfm_arg);
  decode_vml(script_dir,"prop_arg",prop_arg);
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


template<typename MatrixType>
void do_bilinears(Lattice &latt, const int &conf_idx){
  int global_T = GJP.Tnodes()*GJP.TnodeSites();

  //Typically we are going to be dealing with a large number of propagators, plus the eigenvectors, so memory is going to be tight. We therefore don't want to 
  //store all the propagators in memory. Thus we only compute the ones we immediately need from a given source timeslice then throw them away
  
  //Create storage for sums of sloppy and 'rest' (exact-sloppy) bilinear sets for each propagator pair
  int nbils = ama_arg.bilinear_args.bilinear_args_len;
  std::vector< ContractedBilinearSimple<MatrixType> > bil_sum(nbils);
  std::vector< ContractedBilinearSimple<MatrixType> > rest_sum(nbils);

  for(int tsrc = 0; tsrc < global_T; tsrc++){ //timeslice of source
    //Do all sloppy calculations
    std::vector< ContractedBilinearSimple<MatrixType> > bils_t(nbils); 
    
    for(int bil_meas = 0; bil_meas < nbils; bil_meas++){
      if(!UniqueID()){ printf("OUTPUT: Calculating sloppy propagators %d on timeslice %d\n",bil_meas,tsrc); fflush(stdout); }

      ContractionTypeAllBilinears &bil_arg = ama_arg.bilinear_args.bilinear_args_val[bil_meas];
      PropagatorContainer &prop_A = PropManager::getProp(bil_arg.prop_1);
      PropagatorContainer &prop_B = PropManager::getProp(bil_arg.prop_2);
    
      //Get the CGAttr (add if not already present) to set sloppy solve residuals
      CGAttrArg* cg_attr_A;
      CGAttrArg* cg_attr_B;
      if(!prop_A.getAttr(cg_attr_A)) cg_attr_A = addCGattr(prop_A);
      if(!prop_B.getAttr(cg_attr_B)) cg_attr_B = addCGattr(prop_B);

      //Set the sloppy residuals
      cg_attr_A->stop_rsd = ama_arg.sloppy_precision;
      cg_attr_A->true_rsd = ama_arg.sloppy_precision;
    
      cg_attr_B->stop_rsd = ama_arg.sloppy_precision;
      cg_attr_B->true_rsd = ama_arg.sloppy_precision;
    
      //Set source timeslices
      dynamic_cast<QPropWcontainer&>(prop_A).setSourceTimeslice(tsrc);
      dynamic_cast<QPropWcontainer&>(prop_B).setSourceTimeslice(tsrc);

      //Calculate props
      prop_A.calcProp(latt);	
      prop_B.calcProp(latt);

      if(!UniqueID()){ printf("OUTPUT: Calculating sloppy bilinears %d on timeslice %d\n",bil_meas,tsrc); fflush(stdout); }
      
      //Calculate sloppy bilinear
      ContractedBilinearSimple<MatrixType> &bil = bils_t[bil_meas];
      _multimom_helper<ContractedBilinearSimple<MatrixType> >::add_momenta(bil, bil_arg.momenta.momenta_val, bil_arg.momenta.momenta_len);
      bil.calculateBilinears(latt, bil_arg.prop_1, PropDFT::Dagger, bil_arg.prop_2, PropDFT::None); 

      //Temporal shift so prop A always lives at 0 (take advantage of time translation symmetry)
      if(tsrc>0) bil.Tshift(-tsrc);

      bil_sum[bil_meas] += bil;

      //Don't discard the propagators yet as they may need to be re-used for other bilinear measurements
    }
    //Discard all sloppy propagators
    for(int l=0;l<prop_arg.props.props_len;l++)
      PropManager::getProp(prop_arg.props.props_val[l].generics.tag).deleteProp();

    //Now check if we are on a timeslice for an exact solve. If so repeat the above with exact props and compute the 
    bool do_exact = false;
    for(int i=0; i< ama_arg.exact_solve_timeslices.exact_solve_timeslices_len; i++)
      if(tsrc == ama_arg.exact_solve_timeslices.exact_solve_timeslices_val[i]){ do_exact = true; break; }

    if(do_exact){
      for(int bil_meas = 0; bil_meas < nbils; bil_meas++){
	if(!UniqueID()) printf("OUTPUT: Calculating exact propagators %d on timeslice %d\n",bil_meas,tsrc);

	ContractionTypeAllBilinears &bil_arg = ama_arg.bilinear_args.bilinear_args_val[bil_meas];
	PropagatorContainer &prop_A = PropManager::getProp(bil_arg.prop_1);
	PropagatorContainer &prop_B = PropManager::getProp(bil_arg.prop_2);
	
	CGAttrArg* cg_attr_A; prop_A.getAttr(cg_attr_A);
	CGAttrArg* cg_attr_B; prop_B.getAttr(cg_attr_B);

	cg_attr_A->stop_rsd = ama_arg.exact_precision;
	cg_attr_A->true_rsd = ama_arg.exact_precision;
    
	cg_attr_B->stop_rsd = ama_arg.exact_precision;
	cg_attr_B->true_rsd = ama_arg.exact_precision;
    
	prop_A.calcProp(latt);	
	prop_B.calcProp(latt);

	if(!UniqueID()) printf("OUTPUT: Calculating exact bilinears %d on timeslice %d\n",bil_meas,tsrc);

	//Calculate exact bilinear
	ContractedBilinearSimple<MatrixType> ebil;
	_multimom_helper<ContractedBilinearSimple<MatrixType> >::add_momenta(ebil, bil_arg.momenta.momenta_val, bil_arg.momenta.momenta_len);
	ebil.calculateBilinears(latt, bil_arg.prop_1, PropDFT::Dagger, bil_arg.prop_2, PropDFT::None); 

	if(tsrc>0) ebil.Tshift(-tsrc);

	ebil -= bils_t[bil_meas]; //just the 'rest' part

	rest_sum[bil_meas] += ebil;
      }
      for(int l=0;l<prop_arg.props.props_len;l++)
	PropManager::getProp(prop_arg.props.props_val[l].generics.tag).deleteProp();
    }
  }

  //Take timeslice average of the sums
  for(int bil_meas = 0; bil_meas < nbils; bil_meas++){
    bil_sum[bil_meas]/=Float(global_T); //timeslice average
    rest_sum[bil_meas]/=Float(ama_arg.exact_solve_timeslices.exact_solve_timeslices_len);
  }
  //Write out results
  for(int bil_meas = 0; bil_meas < nbils; bil_meas++){
    std::string avg_file;
    { 
      std::ostringstream filestr; 
      filestr << ama_arg.bilinear_args.bilinear_args_val[bil_meas].file << "." << conf_idx; 
      avg_file = filestr.str();
    }
    std::string rest_file;
    { 
      std::ostringstream filestr; 
      filestr << ama_arg.bilinear_args.bilinear_args_val[bil_meas].file << "_rest." << conf_idx; 
      rest_file = filestr.str();
    }   
    bil_sum[bil_meas].write(avg_file);
    rest_sum[bil_meas].write(rest_file);
  }
}
    
    



int main(int argc,char *argv[])
{
  Start(&argc,&argv);
  if(argc<2) ERR.General("","main()","Arguments: Require a directory containing the vml files");
  decode_vml_all(argv[argc-1]);
  if(ama_arg.conf_start >= ama_arg.conf_lessthan || ama_arg.conf_incr == 0) ERR.General("","main()","Invalid configuration args");

  //Props must be QPropW type
  for(int l=0;l<prop_arg.props.props_len;l++){
    if(prop_arg.props.props_val[l].generics.type != QPROPW_TYPE) ERR.General("","main()","Propagators must be QPropW type");
  }

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

  PropManager::setup(prop_arg);
  
  GnoneFbfm lattice;
  CommonArg carg("label","filename");
  char load_config_file[1000];

  //Attach the bfm instance from fbfm to the lanczos solvers so they don't have to maintain separate bfm instances
  for(int l=0;l<prop_arg.lanczos.lanczos_len;l++){
    LanczosContainer& lanc = PropManager::getLanczos(prop_arg.lanczos.lanczos_val[l].tag);
    lanc.set_bfm(&lattice.bd);
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
    for(int l=0;l<prop_arg.lanczos.lanczos_len;l++){
      LanczosContainer& lanc = PropManager::getLanczos(prop_arg.lanczos.lanczos_val[l].tag);
      lanc.calcEig(lattice);

      if(Fbfm::use_mixed_solver){
	//Convert eigenvectors to single precision
	lanc.setPrecision(1);
      }
    }

    if(lanczos_tune) return 0;

    //Compute the AMA bilinears
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

