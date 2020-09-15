//Test deflated propagator inversion using Lanczos eigenvectors

#include <alg/a2a/lanc_arg.h>
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
  
void GaugeTransformU(cps::Matrix *gtrans, Lattice &lat);


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
  dwfa.mass = toDouble(0.01);
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

void lanczos_arg(LancArg &into, const bool &precon = true){
  into.mass = 0.01;
  into.stop_rsd = 1e-10;
  into.qr_rsd = 1e-14; ///convergence of intermediate QR solves, defaults to 1e-14
  into.EigenOper = DDAGD;
  into.precon = precon; //also try this with true

#if 0
  //2^5 lattice DWF
  //Tuning history (following Rudy's prescription from section 146 of his thesis)
  //1) ch_alpha = 2, ch_beta = 0,  ch_ord = 40,  M=10 N=3
// 0:27.0557002268276
// 1:27.0557002268276
// 2:26.8990312652718
// Done : Total time = 0.304079055786133sec, Dlash time = 0.26964545249939sec, Shift time = 0.00147461891174316sec.

  //Now we assume we want 80 eigenvectors. We set beta slightly larger than the largest eigenvalue in this set
  //2) ch_alpha = 28, M=100 N=80, guess beta = 5
  // 76:6.8029698498132
  // 77:6.80296984981321
  // 78:6.95649173285482
  // 79:6.95649173285482
  // total number of vector inner products 17959
  // total number of matrix vector products 19360
  // Done : Total time = 2.27254414558411sec, Dlash time = 0.480455875396729sec, Shift time = 0.18367338180542sec.

  //3) Try the same but with beta = 7
  //Done : Total time = 3.21666383743286sec, Dlash time = 0.700658559799194sec, Shift time = 0.338150024414062sec.

  //4) Try the same but with beta = 7.6
  //Done : Total time = 3.27393388748169sec, Dlash time = 0.624141693115234sec, Shift time = 0.349877595901489sec.

  //5) Weird. Maybe there is some normalization difference? Try beta = 3  (1.34 is the smallest eval)
  //Done : Total time = 1.56976580619812sec, Dlash time = 0.354743719100952sec, Shift time = 0.0855429172515869sec.
  
  //6) Go even lower. Lets do one at the same size as the smallest eval: beta = 1.34
  //Got eigevalues all over the place, time is
  //Done : Total time = 5.13522791862488sec, Dlash time = 1.78019237518311sec, Shift time = 1.03217148780823sec.

  //7) Go to beta = 2
  //didn't converge in reasonable time
  
  //8) Do beta = 2.5.
  //Evals still all over
  //Done : Total time = 7.2279589176178sec, Dlash time = 3.66947960853577sec, Shift time = 1.09042930603027sec.
  
  //9) Beta = 2.8
  //Everything looks good
  //Done : Total time = 1.40910601615906sec, Dlash time = 0.319405555725098sec, Shift time = 0.0684051513671875sec.

  //10) Beta = 2.7
  //Done : Total time = 1.39550805091858sec, Dlash time = 0.317589521408081sec, Shift time = 0.0686178207397461sec.

  //11) Wonder if we gain anything by reducing alpha by, say a factor of 2 to 14
  //Done : Total time = 1.31217694282532sec, Dlash time = 0.316191911697388sec, Shift time = 0.0570271015167236sec.
  //Interesting, definitely seems to be a norm thing. Maybe 5-M5?

  into.N_get = 80;///Want K converged vectors
  into.N_use = 100;///Dimension M Krylov space
  into.N_true_get = 80;//Actually number of eigen vectors you will get
  into.ch_ord = 40;///Order of Chebyshev polynomial
  into.ch_alpha = 14;///Spectral radius
  into.ch_beta = 2.7;///Spectral offset (ie. find eigenvalues of magnitude less than this)
  into.ch_sh = false;///Shifting or not
  into.ch_mu = 0;///Shift the peak
  into.lock = true;///Use locking transofrmation or not
#endif


#if 1
  //4^4 x 2 lattice, 80 low modes  DWF
  
  into.N_get = 80;///Want K converged vectors
  into.N_use = 100;///Dimension M Krylov space
  into.N_true_get = 80;//Actually number of eigen vectors you will get
  into.ch_ord = 40;///Order of Chebyshev polynomial
  into.ch_alpha = 15;///Spectral radius
  into.ch_beta = 1.8;///Spectral offset (ie. find eigenvalues of magnitude less than this)
  into.ch_sh = false;///Shifting or not
  into.ch_mu = 0;///Shift the peak
  into.lock = true;///Use locking transofrmation or not

  //1) ord=40 M=10 K=3 alpha = 2 beta = 0 
// 0:28.5107155128075
// 1:28.5107155128074
// 2:28.4863189089512
// total number of vector inner products 3872
// total number of matrix vector products 38886
// Done : Total time = 15.1352989673615sec, Dlash time = 14.2081019878387sec, Shift time = 0.00699234008789062sec. 

  //2) ord=40 M=100 K=80 alpha=29 beta=5
//   78:1.74565286840916
// 79:1.74565286840916
// Done : Total time = 74.7630920410156sec, Dlash time = 51.0969936847687sec, Shift time = 1.44675445556641sec. 

  //3) ord=40 M=100 K=80 alpha=29 beta=1.8
  //Done : Total time = 26.5963759422302sec, Dlash time = 15.2578508853912sec, Shift time = 0.424580335617065sec.

  //4) ord=40 M=100 K=80 alpha=15 beta=1.8
  //Done : Total time = 16.7637159824371sec, Dlash time = 8.2927839756012sec, Shift time = 0.2242431640625sec.

  //4) ord=40 M=100 K=80 alpha=10 beta=1.8   long time

  //5) ord=40 M=100 K=80 alpha=15 beta=1.0  long time

  //6) ord=40 M=120 K=80 alpha=15 beta=1.8
  //Done : Total time = 20.193972826004sec, Dlash time = 8.3267023563385sec, Shift time = 0.370113611221313sec. 

  //7) ord=50 M=100 K=80 alpha=15 beta=1.8
  //Done : Total time = 19.774062871933sec, Dlash time = 11.1009893417358sec, Shift time = 0.213002681732178sec. 

  //8) ord=30 M=100 K=80 alpha=15 beta=1.8
  //Total time = 18.5850579738617sec, Dlash time = 8.86577439308167sec, Shift time = 0.311829805374146sec.

  //9) ord=40 M=90 K=80 alpha=15 beta=1.8
  //Done : Total time = 20.3546621799469sec, Dlash time = 11.4840474128723sec, Shift time = 0.178162336349487sec.

  //Final is ord=40 M=100 K=80 alpha=15 beta=1.8
#endif


#if 0
  //4^4 x 2 lattice, 300 low modes  DWF
  
  into.N_get = 300;///Want K converged vectors
  into.N_use = 400;///Dimension M Krylov space
  into.N_true_get = 300;//Actually number of eigen vectors you will get
  into.ch_ord = 40;///Order of Chebyshev polynomial
  into.ch_alpha = 15;///Spectral radius
  into.ch_beta = 1.8;///Spectral offset (ie. find eigenvalues of magnitude less than this)
  into.ch_sh = false;///Shifting or not
  into.ch_mu = 0;///Shift the peak
  into.lock = true;///Use locking transofrmation or not

  //1) ord=40 M=10 K=3 alpha = 2 beta = 0 
  //0:28.632083490553
  //1:28.585306005737
  //2:28.5569380488056

  //2) Try best from prev but modified  ord=40 M=400 K=300 alpha=15 beta=1.8
  //Done : Total time = 1002.81136989594sec, Dlash time = 17.8665146827698sec, Shift time = 9.68568181991577sec.
#endif



#if 1
  //4^4 x 2 lattice, 80 low modes  Mobius
  
  into.N_get = 80;///Want K converged vectors
  into.N_use = 100;///Dimension M Krylov space
  into.N_true_get = 80;//Actually number of eigen vectors you will get
  into.ch_ord = 40;///Order of Chebyshev polynomial
  into.ch_alpha = 36;///Spectral radius
  into.ch_beta = 1.8;///Spectral offset (ie. find eigenvalues of magnitude less than this)
  into.ch_sh = false;///Shifting or not
  into.ch_mu = 0;///Shift the peak
  into.lock = true;///Use locking transofrmation or not

  //K=3 M=10, ord=40, alpha=1, beta=0
  //0:85.8714691087006
  //1:85.6246490472095
  //2:85.5030706593131

  //K=80 M=100, ord=40, alpha=86, beta=2
//   77:2.54614238517744
// 78:2.55921702522042
// 79:2.56574359608259
//Done : Total time = 30.604474067688sec, Dlash time = 18.4262278079987sec, Shift time = 0.606805324554443sec.
  
  //K=80 M=100, ord=40, alpha=86, beta=2.6
  //Done : Total time = 32.3841209411621sec, Dlash time = 19.7235660552979sec, Shift time = 0.668405532836914sec.

  //K=80 M=100, ord=40, alpha=86, beta=1.8
  //Done : Total time = 29.5677559375763sec, Dlash time = 17.7810125350952sec, Shift time = 0.590176820755005sec.

  //K=80 M=100, ord=40, alpha=86, beta=1.6
  //too long

  //K=80 M=100, ord=40, alpha=86, beta=1.9
  //Done : Total time = 29.6688230037689sec, Dlash time = 17.7525186538696sec, Shift time = 0.595612525939941sec.

  //K=80 M=100, ord=40, alpha=86, beta=2.2
  //Done : Total time = 30.6043748855591sec, Dlash time = 18.4223539829254sec, Shift time = 0.622102737426758sec.

  //K=80 M=100, ord=40, alpha=43, beta=1.8
  //Done : Total time = 18.5669610500336sec, Dlash time = 9.92488169670105sec, Shift time = 0.248605012893677sec.

  //K=80 M=100, ord=40, alpha=36, beta=1.8
  //Done : Total time = 17.5638499259949sec, Dlash time = 9.31308150291443sec, Shift time = 0.209044456481934sec. 
#endif










  into.maxits =10000;///maxiterations
  into.fname = "Lanczos";






}

void lanczos_container_arg(LanczosContainerArg &arg, const bfmarg &dwfa){
  arg.tag = "lanc";
  lanczos_arg(arg.lanc_arg);
  arg.tbc = GJP.Bc(3);
  //used internally to recreate dwfargs
  arg.cg_max_iter = dwfa.max_iter;
  arg.cg_residual = dwfa.residual;
  arg.cg_precon_5d = dwfa.precon_5d;
  if(dwfa.solver == HmCayleyTanh) arg.solver = BFM_HmCayleyTanh;
  else arg.solver = BFM_DWF;

  arg.mobius_scale = dwfa.mobius_scale;
}

void propagator_arg(PropagatorArg &prop, const char *tag, const bool &deflated = false, const std::string &lanczos_name = ""){
  prop.generics.type = QPROPW_TYPE;
  prop.generics.tag = strdup(tag);
  prop.generics.mass = 0.01;
  for(int d=0;d<4;d++) prop.generics.bc[d] = GJP.Bc(d);

  //Point source
  prop.attributes.attributes_len = deflated ? 2 : 1;
  prop.attributes.attributes_val = new AttributeContainer[prop.attributes.attributes_len];

  AttributeContainer &attr = prop.attributes.attributes_val[0];
  attr.type = POINT_SOURCE_ATTR;
  for(int d =0; d<4; d++) attr.AttributeContainer_u.point_source_attr.pos[d] = 0;

  if(deflated){
    AttributeContainer &attr1 = prop.attributes.attributes_val[1];
    attr1.type = DEFLATED_CG_ATTR;
    attr1.AttributeContainer_u.deflated_cg_attr.lanczos_tag = strdup(lanczos_name.c_str());
  }

}
void test_eigenvectors(Lanczos_5d<double> &lanczos, GnoneFbfm &lattice){
  printf("Testing eigenvectors with Gparity = %d\n",GJP.Gparity());
  bfm_evo<double> &dwf =lanczos.dop; //lattice.bd;

  lattice.BondCond(); //Don't forget to apply the boundary conditions!
  Float* gauge = (Float*) lattice.GaugeField();
  dwf.cps_importGauge(gauge); 

  //long f_size = (long)24 * GJP.VolNodeSites() * GJP.SnodeSites() * (GJP.Gparity()?2:1);
  //Float* cps_tmp = (Float *)pmalloc(sizeof(Float) * f_size);
  Fermion_t Ddag_ev =  dwf.allocFermion() ;
  Fermion_t DDdag_ev = dwf.allocFermion();
  Fermion_t tmp = dwf.allocFermion();


  for(int i=0;i<lanczos.get;i++){
    dwf.Mprec(lanczos.bq[i][1],Ddag_ev,tmp,0);
    dwf.Mprec(Ddag_ev,DDdag_ev,tmp,1);
    //Interesting, they are correct for D D^dag v   but not  D^dag D v

    double normev = sqrt(dwf.norm(lanczos.bq[i][1]));
    double normout = sqrt(dwf.norm(DDdag_ev));

    double gotev = normout/normev;
    double expectev = lanczos.evals[i];

    double ratio = gotev/expectev;

    //Is the solution in the same direction at least?
    dwf.copy(tmp,DDdag_ev);
    dwf.scale(tmp, 1.0/expectev);
    
    double normdiff = dwf.axpy_norm(tmp, lanczos.bq[i][1], tmp, -1.0);

    if(!UniqueID()) printf("Evec %d:  norm in %f, norm out %f, eval %f, expected eval %f (from bl %f). Ratio got/expect = %f. |v - Mv/|Mv| |^2 = %f\n",i,normev,normout,gotev,expectev,lanczos.bl[i],ratio,normdiff);
  }
    
  //lattice.BondCond();

  dwf.freeFermion(Ddag_ev);
  dwf.freeFermion(DDdag_ev);
  dwf.freeFermion(tmp);

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

#ifdef HAVE_BFM

#endif

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

  bool lanczos_tune = false;

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
    }else if( strncmp(cmd,"-lanczos_tune",20) == 0){
      lanczos_tune = true;
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

  Chroma::initialize(&argc,&argv);
  multi1d<int> nrow(Nd);
  
  for(int i = 0; i< Nd; ++i)
    nrow[i] = GJP.Sites(i);
  
  Layout::setLattSize(nrow);
  Layout::create();

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }
  
  //Setup bfm args and pass to Fbfm (static instance)
  bfmarg &dwfa = Fbfm::bfm_args[0];
  BfmSolver solver = HmCayleyTanh; //DWF

  setup_bfmargs(dwfa, solver);

  if(!UniqueID()) printf("Instantiating lattice\n");
  GnoneFbfm lattice;
  if(!UniqueID()) printf("Finished instantiating lattice\n");
					       
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
  //Fbfm is an evil hack that messes with the usual way that things are handled with CPS lattice objects
  //This is because Fbfm applies the BondCond operation when it is constructed, which changes the CPS field as well as the bfm copy.
  //NO OTHER LATTICE CLASS DOES THIS AND IT IS VERY ANNOYING

  //Here I generate the field after constructing the lattice class, so at the moment the Bfm copy and the CPS copy are out of sync
  //I import the CPS field into the bfm field manually now. IT DOES NOT HAVE BCS APPLIED
  lattice.ImportGauge(); //import the new gauge field into FBFM internal bfm instance (also applies BCs to CPS field)

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

  if(gauge_fix){
    lattice.FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    lattice.FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }
  cps_qdp_init(&argc,&argv);
  
  omp_set_num_threads(1);
  
  //Setup the Lanczos solver
  LanczosContainerArg lanc_arg;
  lanczos_container_arg(lanc_arg,dwfa);

  //Give the arguments to the PropManager and generate the eigenvectors
  PropManager::addLanczos(lanc_arg);

  if(!UniqueID()) printf("Computing eigenvectors\n");
  Float time = -dclock();
  PropManager::calcProps(lattice);
  time += dclock();
  print_flops("","Lanczos",0,time);

  Lanczos_5d<double> &lanczos = PropManager::getLanczos(lanc_arg.tag).getEig(lattice); //doesn't recompute

  if(lanczos_tune){
#ifdef HAVE_BFM
    Chroma::finalize();
#endif
    return 0;
  }
    
  //Test the norm of those darned eigenvectors
  test_eigenvectors(lanczos, lattice);


  //Setup a propagator without deflation to measure the time
  PropagatorArg prop_undef;
  propagator_arg(prop_undef,"prop_undef");

  PropManager::addProp(prop_undef);
  PropagatorContainer & prop_undef_pc = PropManager::getProp("prop_undef");

  time = -dclock();
  prop_undef_pc.calcProp(lattice);
  time += dclock();
  print_flops("","Undeflated prop solve",0,time);

  //EDIT: THE BELOW IS NOW DONE INSIDE PROPAGATORCONTAINER
  // //Pass the eigenvectors and eigenvalues to Fbfm
  // int N_use = lanc_arg.lanc_arg.N_true_get;
  // lattice.set_deflation(&lanczos.bq,&lanczos.bl,N_use);
  
  //Now we should be able to use MatInv as usual via the QPropW interface
  PropagatorArg prop_def;
  propagator_arg(prop_def,"prop_def",true,"lanc");

  PropManager::addProp(prop_def);
  PropagatorContainer & prop_def_pc = PropManager::getProp("prop_def");

  time = -dclock();
  prop_def_pc.calcProp(lattice);
  time += dclock();
  print_flops("","Deflated prop solve",0,time);

  Chroma::finalize();
  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}



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
