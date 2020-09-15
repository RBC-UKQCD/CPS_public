//Test the Lanczos-5d code with G-parity BCs

#include <alg/lanc_arg.h>
#include <alg/eigen/Krylov_5d.h>



#include<chroma.h>
//bfm headers
#include<actions/ferm/invert/syssolver_linop_cg_array.h>
#include<bfm.h>
#include<bfm_qdp.h>
#include<bfm_cg.h>
#include<bfm_mprec.h>


//cps headers
#include<alg/fermion_vector.h>
#include<alg/do_arg.h>
#include<alg/meas_arg.h>
#include<util/qioarg.h>
#include<util/ReadLatticePar.h>

//c++ classes
#include<sys/stat.h>
#include<util/qcdio.h>
//#include<fftw3.h>

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

#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>

#include <util/gparity_singletodouble.h>

//some piece of **** defines these elsewhere, so the bfm header gets screwed up
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

#include<omp.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

int nget = 4;
int nuse = 11;

void setup_double_latt(Lattice &double_latt, cps::Matrix* orig_gfield, bool gparity_X, bool gparity_Y){
  //orig latt ( U_0 U_1 ) ( U_2 U_3 ) ( U_4 U_5 ) ( U_6 U_7 )
  //double tatt ( U_0 U_1 U_2 U_3 ) ( U_4 U_5 U_6 U_7 ) ( U_0* U_1* U_2* U_3* ) ( U_4* U_5* U_6* U_7* )

  cps::Matrix *dbl_gfield = double_latt.GaugeField();

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
void setup_double_matrixfield(cps::Matrix* double_mat, cps::Matrix* orig_mat, int nmat_per_site, bool gparity_X, bool gparity_Y){
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
  
void GaugeTransformU(cps::Matrix *gtrans, Lattice &lat);

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


void setup_bfmargs(bfmarg &dwfa, const BfmSolver &solver = DWF){
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
  dwfa.mass = toDouble(0.5);
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

void lanczos_arg(LancArg &into, const bool &precon){
  into.mass = 0.01;
  into.stop_rsd = 1e-06;
  into.qr_rsd = 1e-14; ///convergence of intermediate QR solves, defaults to 1e-14
  into.EigenOper = DDAGD;
  into.precon = precon; //also try this with true
  into.N_get = nget;///Want K converged vectors
  into.N_use = nuse;///Dimension M Krylov space
  into.N_true_get = nget;//Actually number of eigen vectors you will get
  into.ch_ord = 50;///Order of Chebyshev polynomial
  into.ch_alpha = 7;///Spectral radius
  into.ch_beta = 2;///Spectral offset (ie. find eigenvalues of magnitude less than this)
  into.ch_sh = false;///Shifting or not
  into.ch_mu = 0;///Shift the peak
  into.lock = false;///Use locking transofrmation or not
  into.maxits =10000;///maxiterations
  into.fname = "Lanczos";
}

//Returns the number of eigenvectors calculated
int lanczos_test_2f(GwilsonFdwf* lattice, Float** &eigenvectors, std::vector<Float> &evals, const bool &precon){ //assumes GJP and LRG have already been set up appropriately
  bfm_evo<double> dwf;
  bfmarg dwfa;
  setup_bfmargs(dwfa);

  dwf.init(dwfa);
  assert(dwf.gparity == 1);

  lattice->BondCond(); //Don't forget to apply the boundary conditions!
  Float* gauge = (Float*) lattice->GaugeField();
  dwf.cps_importGauge(gauge); 
  lattice->BondCond(); //Don't forget to un-apply the boundary conditions!

  //Setup and run the lancsoz algorithm
  LancArg lanc_arg; lanczos_arg(lanc_arg,precon);
  BFM_Krylov::Lanczos_5d<double> eig(dwf,lanc_arg);
  eig.Run();

  multi1d<bfm_fermion> &eigenvecs = eig.bq;
  
  long f_size = (long)2*24 * GJP.VolNodeSites() * GJP.SnodeSites();

  evals.resize(nget);
  eigenvectors = (Float**)pmalloc(sizeof(Float*)*nget);
  for(int i=0;i<nget;i++){
  //bfm_fermion is an array of 2 Fermion_t
  //For preconditioned solve the memory for the first checkerboard is not allocated
  //which causes the export to unpreconditioned format to SEGV
  //Hence we allocate zeroes to the first checkerboard in this case

    if(precon){ eigenvecs[i][0] = dwf.allocCompactFermion();  dwf.set_zero(eigenvecs[i][0]); }

    eigenvectors[i] = (Float *)pmalloc(sizeof(Float) * f_size);
    dwf.cps_impexFermion(eigenvectors[i],eigenvecs[i],0);
    evals[i] = eig.evals[i];
  }
  return nget;
}
int lanczos_test_1f(GwilsonFdwf* lattice, Float** &eigenvectors,std::vector<Float> &evals,const bool &precon ){ //assumes GJP and LRG have already been set up appropriately
  bfm_evo<double> dwf;
  bfmarg dwfa;
  setup_bfmargs(dwfa);

  dwf.init(dwfa);

  lattice->BondCond(); //Don't forget to apply the boundary conditions!
  Float* gauge = (Float*) lattice->GaugeField();
  dwf.cps_importGauge(gauge); 

  //Setup and run the lancsoz algorithm
  LancArg lanc_arg; lanczos_arg(lanc_arg,precon);
  BFM_Krylov::Lanczos_5d<double> eig(dwf,lanc_arg);
  eig.Run();
  lattice->BondCond(); //Don't forget to un-apply the boundary conditions!

  multi1d<bfm_fermion> &eigenvecs = eig.bq;

  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();
  evals.resize(nget);

  eigenvectors = (Float**)pmalloc(sizeof(Float*)*nget);
  for(int i=0;i<nget;i++){
    if(precon){ eigenvecs[i][0] = dwf.allocCompactFermion();  dwf.set_zero(eigenvecs[i][0]); } //cf above for explanation

    eigenvectors[i] = (Float *)pmalloc(sizeof(Float) * f_size);
    dwf.cps_impexFermion(eigenvectors[i],eigenvecs[i],0);
    evals[i] = eig.evals[i];
  }
  return nget;
}  
void lanczos_test_compare(Float** &_2f_eigenvectors, Float** &_1f_eigenvectors, const bool &gparity_X, const bool &gparity_Y, const int &n_evecs){
  //Convert the 2f eigenvectors onto the 1f lattice format
  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites();

  for(int i=0;i<n_evecs;i++){
    Float *v_dbl = (Float *)pmalloc(sizeof(Float) * f_size);
    setup_double_5d_vector((Vector*)v_dbl,(Vector*)_2f_eigenvectors[i], gparity_X,gparity_Y);
    pfree(_2f_eigenvectors[i]);
    _2f_eigenvectors[i] = v_dbl;
  }
  //Compare them site-by-site
  bool fail(false);
  for(int i=0;i<n_evecs;i++){
    for(int s=0;s<f_size;s++){
      if( fabs(_2f_eigenvectors[i][s] - _1f_eigenvectors[i][s]) > 1e-12 ){
	printf("Lancsoz test fail: evec %d, off %d,  2f: %.12e 1f: %.12e\n",i,s, _2f_eigenvectors[i][s],_1f_eigenvectors[i][s]);
	fail = true;
      }
    }
  }
  if(fail){
    printf("Lanczos test failed\n");
    exit(-1);
  }else{
    printf("Lanczos test passed\n");
  }
}

void test_eigenvectors_unprec(Float** &eigenvectors, const int &n_evecs, GwilsonFdwf* lattice){
  printf("Testing eigenvectors with Gparity = %d\n",GJP.Gparity());
  bfm_evo<double> dwf;
  bfmarg dwfa;
  setup_bfmargs(dwfa);

  dwf.init(dwfa);

  lattice->BondCond(); //Don't forget to apply the boundary conditions!
  Float* gauge = (Float*) lattice->GaugeField();
  dwf.cps_importGauge(gauge); 
  lattice->BondCond(); //Don't forget to un-apply the boundary conditions!

  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites() * (GJP.Gparity()?2:1);
  Float* cps_tmp = (Float *)pmalloc(sizeof(Float) * f_size);

  for(int i=0;i<n_evecs;i++){
    Fermion_t bfm_ev[2] = { dwf.allocFermion(), dwf.allocFermion() };
    dwf.cps_impexFermion(eigenvectors[i],bfm_ev,1);
    Fermion_t Ddag_ev[2] = { dwf.allocFermion(), dwf.allocFermion() };
    Fermion_t DDdag_ev[2] = { dwf.allocFermion(), dwf.allocFermion() };
    Fermion_t tmp = dwf.allocFermion();

    dwf.Munprec(bfm_ev,Ddag_ev,tmp,1);
    dwf.Munprec(Ddag_ev,DDdag_ev,tmp,0);

    dwf.cps_impexFermion(cps_tmp,DDdag_ev,0);

    for(int s=0;s<f_size;s++){
      Float ratio = cps_tmp[s] / eigenvectors[i][s];
      
      printf("Eigenvector %d, offset %d: ratio %.12le\n",i,s,ratio);
    }

    dwf.freeFermion(bfm_ev[0]);    dwf.freeFermion(bfm_ev[1]);
    dwf.freeFermion(Ddag_ev[0]);    dwf.freeFermion(Ddag_ev[1]);
    dwf.freeFermion(DDdag_ev[0]);    dwf.freeFermion(DDdag_ev[1]);
    dwf.freeFermion(tmp);
  }
}

void test_eigenvectors_prec(Float** &eigenvectors, const std::vector<Float> &evals, const int &n_evecs, GwilsonFdwf* lattice){
  if(!UniqueID()){ printf("Testing eigenvectors with Gparity = %d\n",GJP.Gparity()); fflush(stdout); }
  bfm_evo<double> dwf;
  bfmarg dwfa;
  setup_bfmargs(dwfa);

  LancArg lanc_arg; lanczos_arg(lanc_arg,1); //get mass from lanc_arg!
  dwfa.mass = lanc_arg.mass;

  dwf.init(dwfa);
  assert(dwf.gparity == GJP.Gparity());

  lattice->BondCond(); //Don't forget to apply the boundary conditions!
  Float* gauge = (Float*) lattice->GaugeField();
  dwf.cps_importGauge(gauge); 
  lattice->BondCond(); //Don't forget to un-apply the boundary conditions!

  long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites() * (GJP.Gparity()?2:1);
  Float* cps_tmp = (Float *)pmalloc(sizeof(Float) * f_size);

  for(int i=0;i<n_evecs;i++){
    Fermion_t bfm_ev[2] = { dwf.allocFermion(), dwf.allocFermion() };
    dwf.cps_impexFermion(eigenvectors[i],bfm_ev,1);
    Fermion_t Ddag_ev[2] = { dwf.allocFermion(), dwf.allocFermion() };
    Fermion_t DDdag_ev[2] = { dwf.allocFermion(), dwf.allocFermion() };
    Fermion_t tmp = dwf.allocFermion();

    dwf.Mprec(bfm_ev[1],Ddag_ev[1],tmp,0);
    dwf.Mprec(Ddag_ev[1],DDdag_ev[1],tmp,1);

    dwf.copy(DDdag_ev[0],bfm_ev[0]);

    dwf.cps_impexFermion(cps_tmp,DDdag_ev,0);

    double fail = 0;

    for(int s=0;s<f_size;s++){
      if(cps_tmp[s] == 0 && eigenvectors[i][s] == 0){ printf("Eigenvector %d, offset %d: both zero\n",i,s); continue; }
      else if(cps_tmp[s] == 0){ printf("!!Eigenvector %d, offset %d: result is zero, eigenvector is not: %.12le\n",i,s,eigenvectors[i][s]); fail=1.0; continue; }
      else if(eigenvectors[i][s] == 0){ printf("!!Eigenvector %d, offset %d: eigenvector component zero, result is not: %.12le\n",i,s,cps_tmp[s]); fail=1.0; continue; }

      Float ratio = cps_tmp[s] / eigenvectors[i][s];
      
      printf("Eigenvector %d, offset %d: ratio %.12le, eigenvalue is %.12le\n",i,s,ratio,evals[i]);
      if(fabs(ratio - evals[i]) > 1e-05) fail = 1.;
    }
    glb_sum(&fail);
    if(fail > 0.0){ 
      printf("Eigenvector %d does not appear to be an eigenvector!\n",i); fflush(stdout); 
      //exit(-1);
    }

    dwf.freeFermion(bfm_ev[0]);    dwf.freeFermion(bfm_ev[1]);
    dwf.freeFermion(Ddag_ev[0]);    dwf.freeFermion(Ddag_ev[1]);
    dwf.freeFermion(DDdag_ev[0]);    dwf.freeFermion(DDdag_ev[1]);
    dwf.freeFermion(tmp);
  }
}

using namespace Chroma;
using namespace cps;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

int cout_time(char *);
void ReadGaugeField(const MeasArg &meas_arg);
void bfm_init(bfm_evo<double> &dwf,double mq);

int main (int argc,char **argv )
{
  Start(&argc, &argv);
  CommandLine::is(argc,argv);

  bool gparity_X(false);
  bool gparity_Y(false);

  int arg0 = CommandLine::arg_as_int(0);
  printf("Arg0 is %d\n",arg0);
  if(arg0==0){
    gparity_X=true;
    printf("Doing G-parity HMC test in X direction\n");
  }else{
    printf("Doing G-parity HMC test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
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

  int size[] = {2,2,2,2,2};

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

  GJP.Initialize(do_arg);

  LRG.Initialize(); //usually initialised when lattice generated, but I pre-init here so I can load the state from file
  
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

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }


  if(gauge_fix){
    lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    lattice->FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }
  cps_qdp_init(&argc,&argv);
#ifdef HAVE_BFM
  Chroma::initialize(&argc,&argv);
#endif  

  omp_set_num_threads(1);
  //#define UNPREC_TEST
#ifdef UNPREC_TEST
  Float** _2f_evecs_unprec;
  int n_evec_unprec = lanczos_test_2f(lattice, _2f_evecs_unprec,false);
#endif

  #define PREC_TEST
#ifdef PREC_TEST
  std::vector<Float> _2f_evals_prec;
  Float** _2f_evecs_prec;
  int n_evec_prec = lanczos_test_2f(lattice, _2f_evecs_prec,_2f_evals_prec,true); //cf comment at test stage
  test_eigenvectors_prec(_2f_evecs_prec, _2f_evals_prec, n_evec_prec, lattice);
#endif

  if(UniqueID()==0){ printf("Starting double lattice section\n"); fflush(stdout); }
   
  int array_size = 2*lattice->GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
  cps::Matrix *orig_lattice = (cps::Matrix *) pmalloc(array_size);
  memcpy((void*)orig_lattice, (void*)lattice->GaugeField(), array_size);

  lattice->FreeGauge(); //free memory and reset
  delete lattice; //lattice objects are singleton (scope_lock)
  
  //setup 1f model. Upon calling GJP.Initialize the lattice size will be doubled in the appropriate directions
  //and the boundary condition set to APRD
  if(gparity_X) do_arg.gparity_1f_X = 1;
  if(gparity_Y) do_arg.gparity_1f_Y = 1;

  GJP.Initialize(do_arg);

  if(GJP.Gparity()){ printf("Que?\n"); exit(-1); }
  if(UniqueID()==0) printf("Doubled lattice : %d %d %d %d\n", GJP.XnodeSites()*GJP.Xnodes(),GJP.YnodeSites()*GJP.Ynodes(),
			   GJP.ZnodeSites()*GJP.Znodes(),GJP.TnodeSites()*GJP.Tnodes());
  
#ifdef HAVE_BFM
  {
    QDP::multi1d<int> nrow(Nd);  
    for(int i = 0;i<Nd;i++) nrow[i] = GJP.Sites(i);
    //  multi1d<LatticeFermion> test(Nd);  
    //  nrow=size;
    QDP::Layout::setLattSize(nrow);
    QDP::Layout::create();
  }
#endif
  //cps_qdp_init(argc,argv);
  lattice = new GwilsonFdwf;
  setup_double_latt(*lattice,orig_lattice,gparity_X,gparity_Y);
  setup_double_rng(gparity_X,gparity_Y);
 
  if(gauge_fix){
    lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    lattice->FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }
#ifdef UNPREC_TEST
  Float** _1f_evecs_unprec;
  lanczos_test_1f(lattice, _1f_evecs_unprec,false);
#endif
#ifdef PREC_TEST
  Float** _1f_evecs_prec; std::vector<Float> _1f_evals_prec;
  lanczos_test_1f(lattice, _1f_evecs_prec, _1f_evals_prec, true);
  test_eigenvectors_prec(_1f_evecs_prec, _1f_evals_prec, n_evec_prec, lattice);
#endif 

#ifdef UNPREC_TEST
  QDPIO::cout << "Comparing un-preconditioned Lanczos\n";
  lanczos_test_compare(_2f_evecs_unprec,_1f_evecs_unprec,gparity_X,gparity_Y, 4);
#endif  

#ifdef PREC_TEST
  QDPIO::cout << "Comparing preconditioned Lanczos\n";
  lanczos_test_compare(_2f_evecs_prec,_1f_evecs_prec,gparity_X,gparity_Y, 4);
#endif

#ifdef HAVE_BFM
  Chroma::finalize();
#endif

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}






#if 0
	char *cname=argv[0];
	const int TrajStart = atoi(argv[2]);
	const int LessThanLimit = atoi(argv[3]);
	char *fname="main(int,char**)";
	chdir(argv[1]);
	CommonArg common_arg("","");
	DoArg do_arg;
	MeasArg meas_arg;
	FixGaugeArg fix_gauge_arg;
	A2AArg a2a_arg, a2a_arg_s;
	LancArg lanc_arg, lanc_arg_s;

	if(!do_arg.Decode("do_arg.vml","do_arg")){VRB.Result(cname,fname,"Can't open do_arg.vml!\n");exit(1);}
	if(!meas_arg.Decode("meas_arg.vml","meas_arg")){std::cout<<"Can't open meas_arg!"<<std::endl;exit(1);}
	if(!a2a_arg.Decode("a2a_arg.vml","a2a_arg")){VRB.Result(cname,fname,"Can't open a2a_arg.vml!\n");exit(1);}
	if(!a2a_arg_s.Decode("a2a_arg_s.vml","a2a_arg_s")){VRB.Result(cname,fname,"Can't open a2a_arg_s.vml!\n");exit(1);}
	if(!lanc_arg.Decode("lanc_arg.vml","lanc_arg")){VRB.Result(cname,fname,"Can't open lanc_arg.vml!\n");exit(1);}
	if(!lanc_arg_s.Decode("lanc_arg_s.vml","lanc_arg_s")){VRB.Result(cname,fname,"Can't open lanc_arg_s.vml!\n");exit(1);}
	if(!fix_gauge_arg.Decode("fix_gauge_arg.vml","fix_gauge_arg")){VRB.Result(cname,fname,"Can't open fix_gauge_arg.vml!\n");exit(1);}

	Chroma::initialize(&argc,&argv);
	GJP.Initialize(do_arg);
	LRG.Initialize();

	if(lanc_arg.N_true_get!=a2a_arg.nl){
		VRB.Result(cname,fname,"low modes number doesn't match!");
		exit(1);
	}
	if(lanc_arg_s.N_true_get!=a2a_arg_s.nl){
		VRB.Result(cname,fname,"low modes number doesn't match!");
		exit(1);
	}

  /********************************************************
   * Setup QDP
   ********************************************************/
  multi1d<int> nrow(Nd);
  nrow[0] = do_arg.x_sites;
  nrow[1] = do_arg.y_sites;
  nrow[2] = do_arg.z_sites;
  nrow[3] = do_arg.t_sites;
  int Ls = do_arg.s_sites;

  Layout::setLattSize(nrow);
  Layout::create();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  int threads = 64;
  bfmarg::Threads(threads); //This is just telling bfm that we have this many threads on the machine, but not setting real threads to this number.
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);

  /********************************************************
   * Setup DWF operator
   ********************************************************/
  //bfm_qdp<double> dwf;
  bfm_evo<double> dwf; // bfm_evo is another derivative of bfm_base, and it's suitable for gauge and fermion transferring between cps and bfm.

	double M5 = 1.8;

	//Physics parameters
	bfmarg dwfa;
	dwfa.solver       = DWF; //DWFrb4d = 4d preconditioning. DWF = 5d preconditioning
	dwfa.Ls           = GJP.SnodeSites()*GJP.Snodes();
	dwfa.mass         = 0.01;
	dwfa.M5           = M5;
	dwfa.Csw 	    = 0.0;
	dwfa.precon_5d    = 1;
	dwfa.max_iter     = 10000;
	dwfa.residual     = 1e-8;
	//Geometry
	dwfa.node_latt[0] = QDP::Layout::subgridLattSize()[0];
	dwfa.node_latt[1] = QDP::Layout::subgridLattSize()[1];
	dwfa.node_latt[2] = QDP::Layout::subgridLattSize()[2];
	dwfa.node_latt[3] = QDP::Layout::subgridLattSize()[3];

	multi1d<int> procs = QDP::Layout::logicalSize();
	QDPIO::cout << procs.size() << " dim machine\n\t" << endl;
	for(int mu=0;mu<4;mu++){
		QDPIO::cout << procs[mu] << " ";
		if ( procs[mu]>1 ) {
			dwfa.local_comm[mu] = 0;
		} else { 
			dwfa.local_comm[mu] = 1;
		}
		dwfa.ncoor[mu] = 0;
	}
	QDPIO::cout << "\nLocal comm = ";
	for(int mu=0;mu<4;mu++){
		QDPIO::cout << dwfa.local_comm[mu] << " ";
	}
	QDPIO::cout << endl; 

	multi1d<int> ncoor = QDP::Layout::nodeCoord();

	dwf.init(dwfa);


	//bfm_init(dwf,0.032);

 /********************************************************
	* Main Loop Begin!
	********************************************************/
	char dir[200];
	sprintf(dir,"L%d_%d%d_S%d_%d%d_Mu_%1.3f_Ms_%1.3f",a2a_arg.nl,a2a_arg.nhits,a2a_arg.src_width,a2a_arg_s.nl,a2a_arg_s.nhits,a2a_arg_s.src_width,lanc_arg.mass,lanc_arg_s.mass);
	mkdir(dir,0775);
	
	//for(int conf=meas_arg.TrajStart; conf<meas_arg.TrajLessThanLimit; conf+=meas_arg.TrajIncrement)
	for(int conf=TrajStart; conf<LessThanLimit; conf+=meas_arg.TrajIncrement)
	{
		if(!UniqueID()) cout<<"conf="<<conf<<endl;

		meas_arg.TrajCur = conf;

		char comm_arg_filename[200];
		sprintf(comm_arg_filename,"%s/traj_%d",dir,conf);
		common_arg.set_filename(comm_arg_filename);

		
		/********************************************************
		 * Read gauge field
		 ********************************************************/
		ReadGaugeField(meas_arg);
		GimprRectFdwf lat;
		dwf.cps_importGauge((cps::Float *)lat.GaugeField()); // This is a single threads function


		/********************************************************
		 * Light quark v and w
		 ********************************************************/
		Lanczos_5d<double> eig(dwf,lanc_arg);
		cout_time("Light quark low modes begins!");
		eig.Run();
		cout_time("Light quark low modes ends!");
		A2APropbfm a2aprop(lat, a2a_arg, common_arg, &eig);
		a2aprop.allocate_vw();
		a2aprop.compute_vw_low(dwf);
		cout_time("Light quark high modes begins!");
		a2aprop.compute_vw_high(dwf);
		cout_time("Light quark high modes ends!");


		/********************************************************
		 * Strange quark v and w
		 ********************************************************/
		Lanczos_5d<double> eig_s(dwf,lanc_arg_s);
		cout_time("Strange quark low modes begins!");
		eig_s.Run();
		cout_time("Strange quark low modes ends!");
		A2APropbfm a2aprop_s(lat, a2a_arg_s, common_arg, &eig_s);
		a2aprop_s.allocate_vw();
		a2aprop_s.compute_vw_low(dwf);
		cout_time("Strange quark high modes begins!");
		a2aprop_s.compute_vw_high(dwf);
		cout_time("Strange quark high modes ends!");

		/**************************************************
		 * Calculating Meson field.
		 **************************************************/
		AlgFixGauge fix_gauge(lat,&common_arg,&fix_gauge_arg);
		MesonField mf(lat, &a2aprop, &a2aprop_s, &fix_gauge, &common_arg);
		//MesonField mf(lat, &a2aprop, &fix_gauge, &common_arg);

		mf.allocate_vw_fftw();
		cout_time("prepare_vw() begins!");
		mf.prepare_vw();// Gauge fix and fftw each v and w 
		cout_time("prepare_vw() ends!");

		cout_time("Begin cal_mf_ext!");
		mf.cal_mf_ll(2, 1);
		double kaon_rad;
		for(kaon_rad = 2.0; kaon_rad <= 2.1; kaon_rad += 0.2) { 
			mf.cal_mf_sl(kaon_rad, 1);
			mf.cal_mf_ls(kaon_rad, 1);
			mf.run_kaon(kaon_rad);
		}
		mf.cal_mf_ww(2,1,LIGHT,LIGHT); 
		mf.cal_mf_ww(2,1,LIGHT,STRANGE); 
		mf.cal_mf_ww(2,1,STRANGE,LIGHT); 
		mf.free_vw_fftw();
		cout_time("End cal_mf_ext!");
		for(int sep = 1; sep < 5 ; sep++)
			mf.run_pipi(sep);
		//mf.run_type1(12,2);
		mf.run_type1_three_sink_approach(12,2);
		//mf.run_type2(12,2);
		mf.run_type2_three_sink_approach(12,2);
		//mf.run_type3(12,2);
		mf.run_type3_three_sink_approach(12,2);
		//mf.run_type4(12,2);
		mf.run_type4_three_sink_approach(12,2);
	}
	QDPIO::cout<<"Done!"<<endl;
	exit(0);
}
#endif


int cout_time(char *info)
{
	time_t t=time( 0 );
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y/%m/%d %X %A %z",localtime(&t) );
	QDPIO::cout<<tmp<<"\t"<<info<<endl;
	return 0;
}

void ReadGaugeField(const MeasArg &meas_arg)
{
	char *cname = "main";
	char *fname = "ReadGaugeField";

	GnoneFnone lat;
	//  std::stringstream lat_file;
	//  lat_file<<meas_arg.GaugeStem<<'.'<<meas_arg.TrajCur;
	//  QioArg rd_arg(lat_file.str().c_str(),0.001);
	//why do we have to check precision? what is for?
	char lat_file[100];
	sprintf(lat_file,"%s.%d",meas_arg.GaugeStem,meas_arg.TrajCur);
	QioArg rd_arg(lat_file,0.001);

	rd_arg.ConcurIONumber=meas_arg.IOconcurrency;

	ReadLatticeParallel rl;
	rl.read(lat,rd_arg);
	if(!rl.good())ERR.General(cname,fname,"Failed read lattice %s",lat_file);
}
void bfm_init(bfm_evo<double> &dwf,double mq)
{
  int threads = 64;
  bfmarg::Threads(threads); //This is just telling bfm that we have this many threads on the machine, but not setting real threads to this number.
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);

	double M5 = 1.8;

	//Physics parameters
	bfmarg dwfa;
	dwfa.solver       = DWFrb4d; //DWFrb4d = 4d preconditioning. DWF = 5d preconditioning
	dwfa.Ls           = GJP.SnodeSites()*GJP.Snodes();
	dwfa.mass         = mq;
	dwfa.M5           = M5;
	dwfa.Csw 	    = 0.0;
	dwfa.precon_5d    = 0;
	dwfa.max_iter     = 10000;
	dwfa.residual     = 1e-8;
	//Geometry
	dwfa.node_latt[0] = QDP::Layout::subgridLattSize()[0];
	dwfa.node_latt[1] = QDP::Layout::subgridLattSize()[1];
	dwfa.node_latt[2] = QDP::Layout::subgridLattSize()[2];
	dwfa.node_latt[3] = QDP::Layout::subgridLattSize()[3];

	multi1d<int> procs = QDP::Layout::logicalSize();
	QDPIO::cout << procs.size() << " dim machine\n\t" << endl;
	for(int mu=0;mu<4;mu++){
		QDPIO::cout << procs[mu] << " ";
		if ( procs[mu]>1 ) {
			dwfa.local_comm[mu] = 0;
		} else { 
			dwfa.local_comm[mu] = 1;
		}
		dwfa.ncoor[mu] = 0;
	}
	QDPIO::cout << "\nLocal comm = ";
	for(int mu=0;mu<4;mu++){
		QDPIO::cout << dwfa.local_comm[mu] << " ";
	}
	QDPIO::cout << endl; 

	multi1d<int> ncoor = QDP::Layout::nodeCoord();

	dwf.init(dwfa);
}
