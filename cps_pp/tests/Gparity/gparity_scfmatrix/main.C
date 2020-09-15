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

#include <alg/lanc_arg.h>
#include <util/spincolorflavormatrix.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>


#include<sstream>
#include<cassert>
#include<omp.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

static bool equal(const FlavorSpinColorMatrix &a, const SpinColorFlavorMatrix &b){
  for(int f1=0;f1<2;++f1) 
    for(int s1=0;s1<4;++s1) 
      for(int c1=0;c1<3;++c1)
	for(int f2=0;f2<2;++f2) 
	  for(int s2=0;s2<4;++s2) 
	    for(int c2=0;c2<3;++c2){
	      const Rcomplex &aa = a(f1,s1,c1,f2,s2,c2);
	      const Rcomplex &bb = b(s1,c1,f1,s2,c2,f2);

	      if( fabs(aa.real()-bb.real())>1e-12 || fabs(aa.imag()-bb.imag())>1e-12 ) return false;
	    }
  return true;
}
static bool equal(const FlavorSpinColorMatrix &mat1, const FlavorSpinColorMatrix &mat2){
  for(int f1=0;f1<2;++f1) 
    for(int s1=0;s1<4;++s1) 
      for(int c1=0;c1<3;++c1)
	for(int f2=0;f2<2;++f2) 
	  for(int s2=0;s2<4;++s2) 
	    for(int c2=0;c2<3;++c2)
	      if( mat1(f1,s1,c1,f2,s2,c2) != mat2(f1,s1,c1,f2,s2,c2) ){
		printf("Equal fail  %d %d %d %d %d %d: Expect  %f,%f   got   %f,%f\n",f1,s1,c1,f2,s2,c2,mat1(f1,s1,c1,f2,s2,c2).real(), mat1(f1,s1,c1,f2,s2,c2).imag(), mat2(f1,s1,c1,f2,s2,c2).real(), mat2(f1,s1,c1,f2,s2,c2).imag());
		return false;
	      }
  return true;
}
static bool equal(const WilsonMatrix &a, const WilsonMatrix &b){
  for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
  for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
	if( a(s1,c1,s2,c2) != b(s1,c1,s2,c2) ) return false;
  return true;
}



static void random_matrix(FlavorSpinColorMatrix &a, SpinColorFlavorMatrix &b){
  for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
  for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
    Rcomplex c( LRG.Urand(), LRG.Urand() );
    a(f1,s1,c1,f2,s2,c2) = c;
    b(s1,c1,f1,s2,c2,f2) = c;
  }
}


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

  //Test the new FlavorSpinColorMatrix against the old SpinColorFlavorMatrix and time the methods
  assert( sizeof(FlavorSpinColorMatrix) == 1152 * sizeof(Float) );
  
  int ntest = 10000;
  FlavorSpinColorMatrix *fsc = new FlavorSpinColorMatrix[ntest];
  SpinColorFlavorMatrix *scf = new SpinColorFlavorMatrix[ntest];
  Float *rfloat = new Float[2*ntest];

  for(int n=0;n<ntest;++n){
    random_matrix(fsc[n],scf[n]);
    rfloat[2*n] = LRG.Urand();
    rfloat[2*n+1] = LRG.Urand();
  }
  struct timeval start,stop,diff1,diff2;
  
  //Time to loop over all elements
  if(0){
    gettimeofday(&start,NULL);
    Rcomplex f;

    for(int n=0;n<ntest;++n){
      for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
							for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
													  f = fsc[n](f1,s1,c1,f2,s2,c2);
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff1);  
  
    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      for(int s1=0;s1<4;++s1) for(int f1=0;f1<2;++f1) for(int c1=0;c1<3;++c1)
							for(int s2=0;s2<4;++s2) for(int f2=0;f2<2;++f2) for(int c2=0;c2<3;++c2)
													  f = scf[n](s1,c1,f1,s2,c2,f2);
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff2); 
  
    if(!UniqueID()) printf("Loop over element: FSC %d.%6.6d s, SCF %d.%6.6d s\n",diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);
  }


  //Test copy constructor
  if(0){  
    {
      FlavorSpinColorMatrix cp(fsc[0]);
      assert( equal(cp,fsc[0]) );
    }

    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      FlavorSpinColorMatrix cp(fsc[n]); 
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff1);  

    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      SpinColorFlavorMatrix cp(scf[n]); 
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff2);  
  
    if(!UniqueID()) printf("Copy constructor times: FSC %d.%6.6d s, SCF %d.%6.6d s\n",diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);
  }

  //Test construct to float time
  if(0){
    {
      FlavorSpinColorMatrix cp(rfloat[0]);
      for(int s1=0;s1<4;++s1) for(int f1=0;f1<2;++f1) for(int c1=0;c1<3;++c1)
							for(int s2=0;s2<4;++s2) for(int f2=0;f2<2;++f2) for(int c2=0;c2<3;++c2)
													  assert( cp(f1,s1,c1,f2,s2,c2).real() == rfloat[0] && cp(f1,s1,c1,f2,s2,c2).imag() == rfloat[0] );
    }

    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      FlavorSpinColorMatrix cp(rfloat[n]); 
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff1);  

    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      SpinColorFlavorMatrix cp(rfloat[n]); 
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff2);  

    if(!UniqueID()) printf("Float constructor times: FSC %d.%6.6d s, SCF %d.%6.6d s\n",diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);
  }

  //Test operator=
  if(0){
    {
      FlavorSpinColorMatrix cp;
      cp = fsc[0];
      assert( equal(cp,fsc[0]) );
    }

    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      FlavorSpinColorMatrix cp;
      cp = fsc[n]; 
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff1);  

    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      SpinColorFlavorMatrix cp(0);
      cp = scf[n]; 
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff2);  
  
    if(!UniqueID()) printf("operator= times: FSC %d.%6.6d s, SCF %d.%6.6d s\n",diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);
  }

  FlavorMatrixType fmats[] = {F0, F1, Fud, sigma0, sigma1, sigma2, sigma3};
  const char* fmat_nm[] = {"F0","F1","Fud","sigma0","sigma1","sigma2","sigma3"};

  //Test flavor matrix pl
  if(0){

    for(int i=0;i<7;i++){
      FlavorSpinColorMatrix m1(fsc[0]);
      m1.pl(fmats[i]);
    
      SpinColorFlavorMatrix m2(scf[0]);
      m2.pl(fmats[i]);
    
      if(! equal(m1,m2) ){
	if(!UniqueID()) printf("pl test fail for mat s\n",fmat_nm[i]); 
	exit(-1);
      }
    }

    for(int i=0;i<7;i++){
      FlavorSpinColorMatrix *fsc_cp = new FlavorSpinColorMatrix[ntest];
      SpinColorFlavorMatrix *scf_cp = new SpinColorFlavorMatrix[ntest];
      for(int n=0;n<ntest;++n){
	fsc_cp[n] = fsc[n];
	scf_cp[n] = scf[n];
      }

      gettimeofday(&start,NULL);
      for(int n=0;n<ntest;++n){
	fsc_cp[n].pl(fmats[i]);
      }
      gettimeofday(&stop,NULL);
      timersub(&stop,&start,&diff1);  
    
      gettimeofday(&start,NULL);
      for(int n=0;n<ntest;++n){
	scf_cp[n].pl(fmats[i]);
      }
      gettimeofday(&stop,NULL);
      timersub(&stop,&start,&diff2);  
    
      if(!UniqueID()) printf("Left multiply by flavor matrix type %s: FSC %d.%6.6d s, SCF %d.%6.6d s\n",fmat_nm[i],diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);

      delete[] fsc_cp;
      delete[] scf_cp;
    }
  }

  //Test flavor matrix pr
  if(0){

    for(int i=0;i<7;i++){
      FlavorSpinColorMatrix m1(fsc[0]);
      m1.pr(fmats[i]);
    
      SpinColorFlavorMatrix m2(scf[0]);
      m2.pr(fmats[i]);
    
      if(! equal(m1,m2) ){
	if(!UniqueID()) printf("pr test fail for mat %s\n",fmat_nm[i]); 
	exit(-1);
      }
    }

    for(int i=0;i<7;i++){
      FlavorSpinColorMatrix *fsc_cp = new FlavorSpinColorMatrix[ntest];
      SpinColorFlavorMatrix *scf_cp = new SpinColorFlavorMatrix[ntest];
      for(int n=0;n<ntest;++n){
	fsc_cp[n] = fsc[n];
	scf_cp[n] = scf[n];
      }

      gettimeofday(&start,NULL);
      for(int n=0;n<ntest;++n){
	fsc_cp[n].pr(fmats[i]);
      }
      gettimeofday(&stop,NULL);
      timersub(&stop,&start,&diff1);  
    
      gettimeofday(&start,NULL);
      for(int n=0;n<ntest;++n){
	scf_cp[n].pr(fmats[i]);
      }
      gettimeofday(&stop,NULL);
      timersub(&stop,&start,&diff2);  
    
      if(!UniqueID()) printf("Right multiply by flavor matrix type %s: FSC %d.%6.6d s, SCF %d.%6.6d s\n",fmat_nm[i],diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);

      delete[] fsc_cp;
      delete[] scf_cp;
    }
  }

  //Test flavour trace
  if(0){
    {
      WilsonMatrix wm1 = fsc[0].FlavourTrace();
      WilsonMatrix wm2 = scf[0].FlavourTrace();
    
      assert( equal(wm1,wm2) );
    }
 
    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      WilsonMatrix m = fsc[n].FlavourTrace();
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff1);  
  
    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      WilsonMatrix m = scf[n].FlavourTrace();
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff2);  

    if(!UniqueID()) printf("Flavour trace: FSC %d.%6.6d s, SCF %d.%6.6d s\n",diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);
  }

  //Test full trace
  if(0){

    {
      Rcomplex wm1 = fsc[0].Trace();
      Rcomplex wm2 = scf[0].Trace();
      assert( fabs(wm1.real()-wm2.real())<1e-12 && fabs(wm1.imag()-wm2.imag())<1e-12 );
    }
 
    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      Rcomplex m = fsc[n].Trace();
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff1);  
  
    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      Rcomplex m = scf[n].Trace();
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff2);  

    if(!UniqueID()) printf("Full trace: FSC %d.%6.6d s, SCF %d.%6.6d s\n",diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);
  }

  //Test operator* matrix
  if(0){
    {
      FlavorSpinColorMatrix cp = fsc[0]*fsc[1];
      SpinColorFlavorMatrix cp2 = scf[0]*scf[1];
    
      assert( equal(cp,cp2) );
    }
    
    gettimeofday(&start,NULL);
    for(int n=0;n<ntest-1;++n){
      FlavorSpinColorMatrix m = fsc[n] * fsc[n+1];
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff1);  
  
    gettimeofday(&start,NULL);
    for(int n=0;n<ntest-1;++n){
      SpinColorFlavorMatrix m = scf[n] * scf[n+1];
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff2);  

    if(!UniqueID()) printf("operator* matrix: FSC %d.%6.6d s, SCF %d.%6.6d s\n",diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);
  }
  
  //Test operator*= matrix
  if(0){
    {
      FlavorSpinColorMatrix cp = fsc[0];
      cp*=fsc[1];
      SpinColorFlavorMatrix cp2 = scf[0];
      cp2*=scf[1];
    
      assert( equal(cp,cp2) );
    }
    FlavorSpinColorMatrix *fsc_cp = new FlavorSpinColorMatrix[ntest];
    SpinColorFlavorMatrix *scf_cp = new SpinColorFlavorMatrix[ntest];
    for(int n=0;n<ntest;++n){
      fsc_cp[n] = fsc[n];
      scf_cp[n] = scf[n];
    }

    gettimeofday(&start,NULL);
    for(int n=0;n<ntest-1;++n){
      fsc_cp[n]*= fsc[n+1];
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff1);  
  
    gettimeofday(&start,NULL);
    for(int n=0;n<ntest-1;++n){
      scf_cp[n]*= scf[n+1];
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff2);  

    if(!UniqueID()) printf("operator*= matrix: FSC %d.%6.6d s, SCF %d.%6.6d s\n",diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);

    delete[] fsc_cp;
    delete[] scf_cp;
  }

  //Test operator*= float
  if(1){
    {
      FlavorSpinColorMatrix cp = fsc[0];
      cp*=rfloat[0];
      SpinColorFlavorMatrix cp2 = scf[0];
      cp2*=rfloat[0];
    
      assert( equal(cp,cp2) );
    }
    FlavorSpinColorMatrix *fsc_cp = new FlavorSpinColorMatrix[ntest];
    SpinColorFlavorMatrix *scf_cp = new SpinColorFlavorMatrix[ntest];
    for(int n=0;n<ntest;++n){
      fsc_cp[n] = fsc[n];
      scf_cp[n] = scf[n];
    }

    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      fsc_cp[n]*= rfloat[n];
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff1);  
  
    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      scf_cp[n]*= rfloat[n];
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff2);  

    if(!UniqueID()) printf("operator*= float: FSC %d.%6.6d s, SCF %d.%6.6d s\n",diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);

    delete[] fsc_cp;
    delete[] scf_cp;
  }

  //Test operator*= complex
  if(1){
    {
      FlavorSpinColorMatrix cp = fsc[0];
      cp*= Rcomplex(rfloat[0],rfloat[1]);
      SpinColorFlavorMatrix cp2 = scf[0];
      cp2*=Rcomplex(rfloat[0],rfloat[1]);
    
      assert( equal(cp,cp2) );
    }
    FlavorSpinColorMatrix *fsc_cp = new FlavorSpinColorMatrix[ntest];
    SpinColorFlavorMatrix *scf_cp = new SpinColorFlavorMatrix[ntest];
    for(int n=0;n<ntest;++n){
      fsc_cp[n] = fsc[n];
      scf_cp[n] = scf[n];
    }

    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      fsc_cp[n]*= Rcomplex(rfloat[2*n],rfloat[2*n+1]);
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff1);  
  
    gettimeofday(&start,NULL);
    for(int n=0;n<ntest;++n){
      scf_cp[n]*= Rcomplex(rfloat[2*n],rfloat[2*n+1]);
    }
    gettimeofday(&stop,NULL);
    timersub(&stop,&start,&diff2);  

    if(!UniqueID()) printf("operator*= complex: FSC %d.%6.6d s, SCF %d.%6.6d s\n",diff1.tv_sec,diff1.tv_usec, diff2.tv_sec,diff2.tv_usec);

    delete[] fsc_cp;
    delete[] scf_cp;
  }


  delete[] fsc;
  delete[] scf;
  delete[] rfloat;











#ifdef HAVE_BFM
  Chroma::finalize();
#endif

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}


