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
#include <util/dirac_op.h>
#include <util/time_cps.h>
#include <alg/do_arg.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>
#include <alg/alg_plaq.h>
#include <alg/alg_rnd_gauge.h>
#include <alg/threept_arg.h>
#include <alg/threept_prop_arg.h>
#include <alg/alg_threept.h>
#include <util/smalloc.h>

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <util/command_line.h>
#include <sstream>

#include<unistd.h>
#include<config.h>

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
#include <inttypes.h>
//something defines these elsewhere, so the bfm header gets screwed up
#undef ND
#undef SPINOR_SIZE
#undef HALF_SPINOR_SIZE
#undef GAUGE_SIZE
#undef Nmu
#undef Ncb
#undef NMinusPlus
#undef Minus
#undef Plus
#undef DaggerYes
#undef DaggerNo
#undef SingleToDouble
#undef DoubleToSingle
#undef Odd
#undef Even


#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_eigcg.h>
#include <util/lattice/fbfm.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/enum_func.h>
#include <util/sproj_tr.h>
#include <util/time_cps.h>
#include <util/lattice/fforce_wilson_type.h>

#include<sstream>

#include<omp.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

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
void setup_double_5d_vector(Vector *double_vect, Vector* orig_vect, bool gparity_X, bool gparity_Y){
  if(!UniqueID()){ printf("Setting up 1f vector field.\n"); fflush(stdout); }
  SingleToDouble5dVectorField doubler(gparity_X, gparity_Y, orig_vect, double_vect, CANONICAL);
  doubler.Run();
  if(!UniqueID()){ printf("Finished setting up 1f vector field\n"); fflush(stdout); }
}
  
void GaugeTransformU(Matrix *gtrans, Lattice &lat);

void convert_ferm_cpsord_sord(Float *cps, Float* &sord, bfm_evo<Float> &bfm){
  Fermion_t handle[2] = { bfm.allocFermion(), bfm.allocFermion() };
  bfm.cps_impexFermion(cps,handle,1);
  
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  if(GJP.Gparity()) f_size*=2;
  sord = (Float *)pmalloc(sizeof(Float) * f_size);
  bfm.cps_impexFermion_s(sord,handle,0);

  bfm.freeFermion(handle[0]);
  bfm.freeFermion(handle[1]);
}
void convert_ferm_sord_cpsord(Float *sord, Float* &cps, bfm_evo<Float> &bfm){
  Fermion_t handle[2] = { bfm.allocFermion(), bfm.allocFermion() };
  bfm.cps_impexFermion_s(sord,handle,1);
  
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  if(GJP.Gparity()) f_size*=2;
  cps = (Float *)pmalloc(sizeof(Float) * f_size);
  bfm.cps_impexFermion(cps,handle,0);

  bfm.freeFermion(handle[0]);
  bfm.freeFermion(handle[1]);
}


void setup_bfmargs(bfmarg &dwfa, const BfmSolver &solver){
  printf("Setting up bfmargs\n");

   int nthreads = 1; 
#if TARGET == BGQ
   nthreads = 64;
#endif
   omp_set_num_threads(nthreads);

  dwfa.node_latt[0]  = GJP.XnodeSites();
  dwfa.node_latt[1]  = GJP.YnodeSites();
  dwfa.node_latt[2]  = GJP.ZnodeSites();
  dwfa.node_latt[3]  = GJP.TnodeSites();
  
  multi1d<int> ncoor(4);
  multi1d<int> procs(4);
  for(int i=0;i<4;i++){ ncoor[i] = GJP.NodeCoor(i); procs[i] = GJP.Nodes(i); }

  if(GJP.Gparity()){
    dwfa.gparity = 1;
    printf("G-parity directions: ");
    for(int d=0;d<3;d++)
      if(GJP.Bc(d) == BND_CND_GPARITY){ dwfa.gparity_dir[d] = 1; printf("%d ",d); }
      else dwfa.gparity_dir[d] = 0;
    for(int d=0;d<4;d++){
      dwfa.nodes[d] = procs[d];
      dwfa.ncoor[d] = ncoor[d];
    }
    printf("\n");
  }

  dwfa.verbose=1;
  dwfa.reproduce=0;
  bfmarg::Threads(nthreads);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);

  for(int mu=0;mu<4;mu++){
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
      printf("Non-local comms in direction %d\n",mu);
    } else { 
      dwfa.local_comm[mu] = 1;
      printf("Local comms in direction %d\n",mu);
    }
  }

  dwfa.precon_5d = 1;
  if(solver == HmCayleyTanh){
    dwfa.precon_5d = 0; //mobius uses 4d preconditioning
    dwfa.mobius_scale = 2.0; //b = 0.5(scale+1) c=0.5(scale-1), hence this corresponds to b=1.5 and c=0.5, the params used for the 48^3
  }
  dwfa.Ls   = GJP.SnodeSites();
  dwfa.solver = solver;
  dwfa.M5   = toDouble(GJP.DwfHeight());
  dwfa.mass = toDouble(0.001);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 5000;
  dwfa.residual = 1e-08;
  printf("Finished setting up bfmargs\n");
}

Float* rand_5d_canonical_fermion(Lattice &lat){
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  if(GJP.Gparity()) f_size*=2;
  Float *v1 = (Float *)pmalloc(sizeof(Float) * f_size);
  printf("Making random gaussian 5d vector\n");
  lat.RandGaussVector((Vector*)v1, 0.5, 2, CANONICAL, FIVE_D);
  printf("Finished making random gaussian vector\n");
  return v1;
}

#include <bfm_vmx.h>

template<typename PREC>
static void cost_breakdown(Lattice* lattice, const BfmSolver &solver){
  //Time the various operations involved in the CG
  
  bfmarg dwfa;
  setup_bfmargs(dwfa,solver);
  
  bfm_evo<PREC> bfm;
  bfm.init(dwfa);
  bfm.verbose = 1;
    
  lattice->BondCond();
  Float* gauge = (Float*) lattice->GaugeField();
  bfm.cps_importGauge(gauge);  

  LatRanGen LRGbak(LRG);
  Float* v1 = rand_5d_canonical_fermion(*lattice);
  printf("Restoring RNG\n"); fflush(stdout);
  LRG = LRGbak;

  printf("Allocating fermions\n"); fflush(stdout);
  Fermion_t src[2] = {bfm.allocFermion(), bfm.allocFermion()}; //odd/even

  Fermion_t tmp1 = bfm.allocFermion();
  Fermion_t tmp2 = bfm.allocFermion();

  bfm.cps_impexFermion(v1,src,1);

  //typedef long long unsigned int lint;
#pragma omp parallel
  {
    struct timeval start, stop, diff;
    
    double sum = 0.0, sumsq = 0.0;
    int niter = 1000;

    //Do Mprecs
    for(int i=0;i<niter;i++){
      gettimeofday(&start,NULL); 
      bfm.Mprec(src[0],tmp1,tmp2,0,1);
      gettimeofday(&stop,NULL);
      timersub(&stop,&start,&diff);
      
      double time = (double)diff.tv_sec * 1000000 + (double)diff.tv_usec;
      sum += time;
      sumsq += time*time;
    }

    double avg = sum / double(niter);
    double stddev = sqrt(sumsq/ double(niter) - avg*avg);
    double stderr = stddev / sqrt(double(niter));

    int me = omp_get_thread_num();

    if(!me && !UniqueID()){
      printf("Mprec time avg %.3f us, std.err %.3f us\n",avg,stderr);
    }

    
    //Do axpy_norms
    sum=0.0; sumsq=0.0;
    for(int i=0;i<niter;i++){
      gettimeofday(&start,NULL); 
      bfm.axpy_norm(src[0],src[1],src[0],1e-07);
      
      gettimeofday(&stop,NULL);
      timersub(&stop,&start,&diff);
      
      double time = (double)diff.tv_sec * 1000000 + (double)diff.tv_usec;
      sum += time;
      sumsq += time*time;
    }

    avg = sum / double(niter);
    stddev = sqrt(sumsq/ double(niter) - avg*avg);
    stderr = stddev / sqrt(double(niter));

    if(!me && !UniqueID()){
      printf("axpy_norm time avg %.3f us, std.err %.3f us\n",avg,stderr);
    }

    //Do axpy
    sum=0.0; sumsq=0.0;
    for(int i=0;i<niter;i++){
      gettimeofday(&start,NULL); 
      bfm.axpy(src[0],src[1],src[0],1e-07);
      
      gettimeofday(&stop,NULL);
      timersub(&stop,&start,&diff);
      
      double time = (double)diff.tv_sec * 1000000 + (double)diff.tv_usec;
      sum += time;
      sumsq += time*time;
    }

    avg = sum / double(niter);
    stddev = sqrt(sumsq/ double(niter) - avg*avg);
    stderr = stddev / sqrt(double(niter));

    if(!me && !UniqueID()){
      printf("axpy time avg %.3f us, std.err %.3f us\n",avg,stderr);
    }

    //Do axpby
    sum=0.0; sumsq=0.0;
    for(int i=0;i<niter;i++){
      gettimeofday(&start,NULL); 
      bfm.axpby(src[0],src[1],src[0],1e-07,1.0000001);
      
      gettimeofday(&stop,NULL);
      timersub(&stop,&start,&diff);
      
      double time = (double)diff.tv_sec * 1000000 + (double)diff.tv_usec;
      sum += time;
      sumsq += time*time;
    }

    avg = sum / double(niter);
    stddev = sqrt(sumsq/ double(niter) - avg*avg);
    stderr = stddev / sqrt(double(niter));

    if(!me && !UniqueID()){
      printf("axpby time avg %.3f us, std.err %.3f us\n",avg,stderr);
    }

    //Do norm
    sum=0.0; sumsq=0.0;
    for(int i=0;i<niter;i++){
      gettimeofday(&start,NULL); 
      bfm.norm(src[0]);
      
      gettimeofday(&stop,NULL);
      timersub(&stop,&start,&diff);
      
      double time = (double)diff.tv_sec * 1000000 + (double)diff.tv_usec;
      sum += time;
      sumsq += time*time;
    }

    avg = sum / double(niter);
    stddev = sqrt(sumsq/ double(niter) - avg*avg);
    stderr = stddev / sqrt(double(niter));

    if(!me && !UniqueID()){
      printf("norm time avg %.3f us, std.err %.3f us\n",avg,stderr);
    }


  }

  bfm.freeFermion(tmp1);
  bfm.freeFermion(tmp2);
  bfm.freeFermion(src[0]);
  bfm.freeFermion(src[1]);
  pfree(v1);
}


#if 1
int main(int argc,char *argv[])
{
  Start(&argc,&argv); //initialises QMP

#ifdef HAVE_BFM
  Chroma::initialize(&argc,&argv);
#endif

  CommandLine::is(argc,argv);

  bool gparity_X(false);
  bool gparity_Y(false);

  int arg0 = CommandLine::arg_as_int(0);
  printf("Arg0 is %d\n",arg0);
  if(arg0==0){
    gparity_X=true;
    printf("Doing G-parity HMC test in X direction\n");
  }else if(arg0==1){
    printf("Doing G-parity HMC test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
  }else{
    printf("Doing No G-parity test\n");
  }

  bool dbl_latt_storemode(false);
  bool save_config(false);
  bool load_config(false);
  bool load_lrg(false);
  bool save_lrg(false);
  char *load_config_file;
  char *save_config_file;
  char *save_lrg_file;
  char *load_lrg_file;
  bool gauge_fix(false);
  bool verbose(false);
  bool skip_gparity_inversion(false);
  bool unit_gauge(false);

  BfmSolver solver = DWF;

  int size[] = {2,2,2,2,2};

  double min_fp_resid = 1e-05;

  bool single = false;

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
      size[0] = CommandLine::arg_as_int(i); //CommandLine ignores zeroth input arg (i.e. executable name)
      size[1] = CommandLine::arg_as_int(i+1);
      size[2] = CommandLine::arg_as_int(i+2);
      size[3] = CommandLine::arg_as_int(i+3);
      size[4] = CommandLine::arg_as_int(i+4);
      i+=6;
    }else if( strncmp(cmd,"-save_double_latt",20) == 0){
      dbl_latt_storemode = true;
      i++;
    }else if( strncmp(cmd,"-load_lrg",15) == 0){
      if(i==argc-1){ printf("-load_lrg requires an argument\n"); exit(-1); }
      load_lrg=true;
      load_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-min_fp_resid",15) == 0){
      std::stringstream ss;
      ss << argv[i+1];
      ss >> min_fp_resid;
      if(UniqueID()) printf("Set minimum floating point residual for mixed-prec multi-mass shift to %e\n",min_fp_resid);
      i+=2;
    }else if( strncmp(cmd,"-save_lrg",15) == 0){
      if(i==argc-1){ printf("-save_lrg requires an argument\n"); exit(-1); }
      save_lrg=true;
      save_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-gauge_fix",15) == 0){
      gauge_fix=true;
      i++;   
    }else if( strncmp(cmd,"-verbose",15) == 0){
      verbose=true;
      i++;
    }else if( strncmp(cmd,"-skip_gparity_inversion",30) == 0){
      skip_gparity_inversion=true;
      i++;
    }else if( strncmp(cmd,"-unit_gauge",15) == 0){
      unit_gauge=true;
      i++;
    }else if( strncmp(cmd,"-mobius",15) == 0){
      solver= HmCayleyTanh;
      i++;
    }else if( strncmp(cmd,"-single",15) == 0){
      single = true;
      i++;
      if(UniqueID()) printf("Doing single precision\n");
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }
  

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
  do_arg.start_conf_kind = START_CONF_ORD;
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
  do_arg.verbose_level = -1202; //VERBOSE_DEBUG_LEVEL; //-1202;
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

  if(verbose) do_arg.verbose_level = VERBOSE_DEBUG_LEVEL;

  if(gparity_X) do_arg.x_bc = BND_CND_GPARITY;
  if(gparity_Y) do_arg.y_bc = BND_CND_GPARITY;

  //if(!gparity_X && !gparity_Y) ERR.General("","","Must have G-parity in at least one direction!\n");

  GJP.Initialize(do_arg);

  SerialIO::dbl_latt_storemode = dbl_latt_storemode;
  
  LRG.Initialize(); //usually initialised when lattice generated, but I pre-init here so I can load the state from file

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }
  
  GwilsonFdwf* lattice = new GwilsonFdwf;
					       
  if(!load_config){
    printf("Creating gauge field\n");
    if(!unit_gauge) lattice->SetGfieldDisOrd();
    else lattice->SetGfieldOrd();
  }else{
    ReadLatticeParallel readLat;
    if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
    readLat.read(*lattice,load_config_file);
    if(UniqueID()==0) printf("Config read.\n");
  }

  if(save_config){
    if(UniqueID()==0) printf("Saving config to %s\n",save_config_file);

    QioArg wt_arg(save_config_file,0.001);
    
    wt_arg.ConcurIONumber=32;
    WriteLatticeParallel wl;
    wl.setHeader("disord_id","disord_label",0);
    wl.write(*lattice,wt_arg);
    
    if(!wl.good()) ERR.General("main","()","Failed write lattice %s",save_config_file);

    if(UniqueID()==0) printf("Config written.\n");
  }

  if(gauge_fix){
    lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    lattice->FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }
  cps_qdp_init(&argc,&argv);

  if(single){
    cost_breakdown<float>(lattice,solver);
  }else{
    cost_breakdown<double>(lattice,solver);
  }


#ifdef HAVE_BFM
  Chroma::finalize();
#endif

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}




void GaugeTransformU(Matrix *gtrans, Lattice &lat){
  Matrix recv_buf;
  Matrix tmp;
  //apply the gauge transformation to U
  int nflav = 1;
  if(GJP.Gparity()) nflav = 2;

  for(int flav=0;flav<nflav;flav++){
    for(int t=0;t<GJP.TnodeSites();t++){
      for(int z=0;z<GJP.ZnodeSites();z++){
	for(int y=0;y<GJP.YnodeSites();y++){
	  for(int x=0;x<GJP.XnodeSites();x++){
	    int pos[4] = {x,y,z,t};
	    int v_x_off = x + GJP.XnodeSites()*(y+GJP.YnodeSites()*(z+GJP.ZnodeSites()*t)) + flav*GJP.VolNodeSites();
	    Matrix &v_x = *(gtrans + v_x_off);

	    for(int mu=0;mu<4;mu++){
	      int u_x_off = lat.GsiteOffset(pos) + mu + flav*4*GJP.VolNodeSites();
	      Matrix &u_x = *(lat.GaugeField() + u_x_off);

	      //get V_x+mu
	      int posp[4] = {x,y,z,t};
	      posp[mu] = (posp[mu]+1)%GJP.NodeSites(mu);

	      Matrix *v_xpmu_ptr = gtrans + posp[0] + GJP.XnodeSites()*(posp[1]+GJP.YnodeSites()*(posp[2]+GJP.ZnodeSites()*posp[3])) + flav*GJP.VolNodeSites();
	      if(pos[mu] == GJP.NodeSites(mu)-1){
		//if node is on the left wall, send the opposite flavour 
		if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu) == 0){
		  if(flav == 1)
		    v_xpmu_ptr-= GJP.VolNodeSites();
		  else
		    v_xpmu_ptr+= GJP.VolNodeSites();		  
		}

		//doesnt need to be fast!
		getPlusData((double *)&recv_buf, (double *)v_xpmu_ptr, 18, mu);
		v_xpmu_ptr = &recv_buf; 
	      }

	      //dagger/transpose it
	      Matrix vdag_xpmu;
	      vdag_xpmu.Dagger(*v_xpmu_ptr);

	      //gauge transform link
	      tmp.DotMEqual(v_x,u_x);
	      u_x.DotMEqual(tmp,vdag_xpmu);
	    }
	  }
	}
      }
    }

  }

}

#endif
