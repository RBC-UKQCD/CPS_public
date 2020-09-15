//CK: In this test we check that the 1-flavour doubled-lattice approach gives the same HMC measurements (psibar-psi, plaquette) as the 2-flavour single-lattice approach

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
#if(0==1)
 #include <ReadLattice.h>
 #include <WriteLattice.h>
#endif

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
#include <alg/alg_pbp.h>

#include <util/data_shift.h>

#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>

#include <util/gparity_singletodouble.h>
#include <alg/alg_meas.h>

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
void GaugeTransformU(Matrix *gtrans, Lattice &lat);


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
  
  //GwilsonFdwf* lattice = new GwilsonFdwf;
  {
    Lattice &lattice = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_IMPR_RECT); 
    				       
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
    LatticeFactory::Destroy();
  }

#define ALG_PLAQ_TEST
#ifdef ALG_PLAQ_TEST
  
  CommonArg common_arg_plaq;
  NoArg no_arg;

  Float _2f_plaq[6];
  {
    LatRanGen LRGbak(LRG);

    common_arg_plaq.set_filename("plaq.2f.dat");

    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON); //*lattice;
    AlgPlaq plaq(lat, &common_arg_plaq, &no_arg);
    plaq.run(_2f_plaq);
    LatticeFactory::Destroy();

    LRG = LRGbak;
  }
#endif

#define ALG_PBP_TEST
#ifdef ALG_PBP_TEST
  CommonArg common_arg_pbp;
  PbpArg pbp_arg;

  pbp_arg.pattern_kind = ARRAY;
  pbp_arg.n_masses = 2;
  pbp_arg.mass[0] = 0.4;
  pbp_arg.mass[1] = 0.5;
  pbp_arg.max_num_iter = 5000;
  pbp_arg.stop_rsd = 1e-06;
  pbp_arg.src_u_s = 0;
  pbp_arg.src_l_s = GJP.SnodeSites()-1;
  pbp_arg.snk_u_s = 0;
  pbp_arg.snk_l_s  = GJP.SnodeSites()-1;
  pbp_arg.snk_loop = 0;

  int pbpsz = 2*pbp_arg.n_masses;

  Float _2f_pbp[pbpsz];
  {
    LatRanGen LRGbak(LRG);

    common_arg_pbp.set_filename("pbp.2f.dat");
    //Lattice &lat = *lattice;
    Lattice &lat = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_WILSON);
    AlgPbp pbp(lat,&common_arg_pbp,&pbp_arg);
    pbp.run(_2f_pbp);
    LatticeFactory::Destroy();
    LRG = LRGbak;
  }
#endif

  if(UniqueID()==0) printf("Starting double lattice inversion\n");
  Matrix *orig_lattice;
  int array_size;

  {
    Lattice &lattice = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_IMPR_RECT); 
    if(gauge_fix) lattice.FixGaugeFree();

    array_size = 2*lattice.GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
    orig_lattice = (Matrix *) pmalloc(array_size);
    memcpy((void*)orig_lattice, (void*)lattice.GaugeField(), array_size);
    
    lattice.FreeGauge(); //free memory and reset
    LatticeFactory::Destroy();
  }

  //delete lattice; //lattice objects are singleton (scope_lock)
  
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

  {
    Lattice &doubled_lattice = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_IMPR_RECT); 
    setup_double_latt(doubled_lattice,orig_lattice,gparity_X,gparity_Y);
    setup_double_rng(gparity_X,gparity_Y);
    
    if(gauge_fix){
      doubled_lattice.FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
      doubled_lattice.FixGauge(1e-06,2000);
      if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
    }
    LatticeFactory::Destroy();
  }
  
#ifdef ALG_PLAQ_TEST
  {
    LatRanGen LRGbak(LRG);
    Float _1f_plaq[6];

    common_arg_plaq.set_filename("plaq.1f.dat");
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);  //doubled_lattice;
    AlgPlaq plaq(lat, &common_arg_plaq, &no_arg);
    plaq.run(_1f_plaq);
    LatticeFactory::Destroy();
    Float err(0.0);
    for(int i=0;i<6;i++){
      if(fabs(_1f_plaq[i] - _2f_plaq[i])>1e-08){
	if(!UniqueID()) printf("Error: plaq %d difference: 2f:1f %f:%f\n",i,_2f_plaq[i],_1f_plaq[i]);
	err = 1.0;
      }
    }
    fflush(stdout);

    glb_sum(&err);

    LRG = LRGbak;
    if(err>0.0){
      printf("Plaq check failed\n");
      exit(-1);
    }else printf("Plaq check passed\n");
  }
#endif

#ifdef ALG_PBP_TEST
  {
    LatRanGen LRGbak(LRG);
    Float _1f_pbp[pbpsz];

    common_arg_pbp.set_filename("pbp.1f.dat");
    //Lattice &lat = doubled_lattice;
    Lattice &lat = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_WILSON);
    AlgPbp pbp(lat,&common_arg_pbp,&pbp_arg);
    pbp.run(_1f_pbp);
    LatticeFactory::Destroy();

    Float err(0.0);
    for(int i=0;i<pbpsz;i++){
      if(fabs(_1f_pbp[i] - _2f_pbp[i])>1e-08){
	if(!UniqueID()) printf("Error: pbp %d difference: 2f:1f %f:%f\n",i,_2f_pbp[i],_1f_pbp[i]);
	err = 1.0;
      }else if(!UniqueID()) printf("Success: pbp %d match: 2f:1f %f:%f\n",i,_2f_pbp[i],_1f_pbp[i]);
    }

    glb_sum(&err);

    LRG = LRGbak;
    if(err>0.0){
      printf("Pbp check failed\n");
      exit(-1);
    }else printf("Pbp check passed\n");
  }
#endif

  {
    Lattice &doubled_lattice = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_IMPR_RECT); 
    if(gauge_fix) doubled_lattice.FixGaugeFree();
    LatticeFactory::Destroy();
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



















// void print(const WilsonMatrix &w){
//   for(int i=0;i<4;i++){
//     for(int j=0;j<4;j++){
//       Complex c = w(i,0,j,0);
//       printf("(%.2f %.2f) ",c.real(),c.imag());
//     }
//     printf("\n");
//   }
//   printf("\n");
// }

// bool test_equals(const WilsonMatrix &a, const WilsonMatrix &b, const double &eps){
//   for(int i=0;i<4;i++){
//     for(int j=0;j<4;j++){
//       for(int aa=0;aa<3;aa++){
// 	for(int bb=0;bb<3;bb++){
// 	  Complex ca = a(i,aa,j,bb);
// 	  Complex cb = b(i,aa,j,bb);
// 	  if( fabs(ca.real()-cb.real()) > eps || fabs(ca.imag()-cb.imag()) > eps ) return false;
// 	}
//       }
//     }
//   }
//   return true;
// }





//   if(gparity_X && !gparity_Y){
//     if(GJP.Xnodes()>1){
//       

//       SingleToDoubleLattice lattdoubler(gparity_X,gparity_Y,orig_gfield,double_latt);
//       lattdoubler.Run();

//       if(!UniqueID()){ printf("Finished setting up doubled lattice\n"); fflush(stdout); }
//     }else{
//       //only one node in X-direction
//       //copy data from orig_latt stored on this node
//       int pos[4];
//       for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]++){
// 	for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]++){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]++){
// 	    for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]++){
// 	      for(int mu=0;mu<4;mu++){
// 		int dbl_site_idx = 4*(pos[0]+GJP.XnodeSites() *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		if(pos[0] < GJP.XnodeSites()/2){
// 		  //site is stored on-node in orig_gfield
// 		  int orig_site_idx = 4*(pos[0]+GJP.XnodeSites()/2 *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		  dbl_gfield[dbl_site_idx] = orig_gfield[orig_site_idx];
// 		}else{
// 		  //site is stored on-node in orig_gfield but needs complex conjugating
// 		  int orig_site_idx = 4*(pos[0]-GJP.XnodeSites()/2 +GJP.XnodeSites()/2 *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		  dbl_gfield[dbl_site_idx].Conj(orig_gfield[orig_site_idx]);
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }
// #if 0
//       printf("Converted single to double lattice. Scan across X-direction of original:\n");
//       for(int mu=0;mu<4;mu++){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]++){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]++){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]++){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites()/2;pos[0]++){
// 		int orig_site_idx = 4*(pos[0]+GJP.XnodeSites()/2 *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		printf("U(%d,%d,%d,%d)_%d = (%f,%f) ",pos[0],pos[1],pos[2],pos[3],mu,*((IFloat*)(orig_gfield+orig_site_idx)),*((IFloat*)(orig_gfield+orig_site_idx)+1));
// 	      }
// 	      printf("\n");
// 	    }
// 	  }
// 	}
//       }
//       printf("\nDoubled lattice:\n");
//       for(int mu=0;mu<4;mu++){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]++){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]++){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]++){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]++){
// 		int dbl_site_idx = 4*(pos[0]+GJP.XnodeSites() *(pos[1]+GJP.YnodeSites()*(pos[2]+GJP.ZnodeSites()*pos[3])))+mu;
// 		printf("U(%d,%d,%d,%d)_%d = (%f,%f) ",pos[0],pos[1],pos[2],pos[3],mu,*((IFloat*)(dbl_gfield+dbl_site_idx)),*((IFloat*)(dbl_gfield+dbl_site_idx)+1));
// 	      }
// 	      printf("\n");
// 	    }
// 	  }
// 	}
// 	printf("\n");
//       }
// #endif
//     }


//   }else if(gparity_X && gparity_Y){
//     if(!UniqueID()){ printf("Setting up quad lattice. sizeof(Float) %d sizeof(IFloat) %d\n",sizeof(Float), sizeof(IFloat)); fflush(stdout); }

//       SingleToDoubleLattice lattdoubler(gparity_X,gparity_Y,orig_gfield,double_latt);
//       lattdoubler.Run();

//       if(!UniqueID()){ printf("Finished setting up quad lattice\n"); fflush(stdout); }
//   }

// }












  // if(gparity_X && !gparity_Y){
  //   if(!UniqueID()) printf("Setting up RNG from original stacked version\n");


  // //orig 4D rng 2 stacked 4D volumes
  // //orig   ([R_0 R_1][R'_0 R'_1]) ([R_2 R_3][R'_2 R'_3]) ([R_4 R_5][R'_4 R'_5]) ([R_6 R_7][R'_6 R'_7]) ([R_8 R_9][R'_8 R'_9]) ([R_10 R_11][R'_10 R'_11]) ([R_12 R_13][R'_12 R'_13]) ([R_14 R_15][R'_14 R'_15])
  // //double (R_0 R_1 R_2 R_3)      (R_4 R_5 R_6 R_7)      (R_8 R_9 R_10 R_11)    (R_12 R_13 R_13 R_15)  (R'_0 R'_1 R'_2 R'_3)  (R'_4 R'_5 R'_6 R'_7)      (R'_8 R'_9 R'_10 R'_11)    (R'_12 R'_13 R'_14 R'_15)

  //   if(GJP.Xnodes()>1){
  //     SingleToDouble4dRNG fourDsetup;
  //     SingleToDouble5dRNG fiveDsetup;

  //     LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form
      
  //     if(!UniqueID()){ printf("Setting up 4D RNG\n"); fflush(stdout); }
  //     fourDsetup.Run(gparity_X,gparity_Y);      
  //     if(!UniqueID()){ printf("Setting up 5D RNG\n"); fflush(stdout); }
  //     fiveDsetup.Run(gparity_X,gparity_Y);    
  //   }else{    
  //     int n_rgen_4d = GJP.VolNodeSites()/16; //applies both to original and doubled latt
  //     int n_rgen = n_rgen_4d;
  //     if (GJP.SnodeSites()>=2)
  // 	n_rgen = GJP.VolNodeSites()*GJP.SnodeSites() / 32;

  //     int stk_index_4d_off = n_rgen_4d/2; //offset for R' on 4D orig latt
  //     int blocks_per_s_layer = n_rgen /( GJP.SnodeSites() / 2 ); //also same for original and doubled latt
  //     int stk_index_5d_off = blocks_per_s_layer/2; //offset for R' on 5D orig latt

  //     //copy the originals
  //     UGrandomGenerator *ugran_4d_orig = new UGrandomGenerator[n_rgen_4d];
  //     for(int i=0;i<n_rgen_4d;i++) ugran_4d_orig[i] = LRG.UGrandGen4D(i);

  //     UGrandomGenerator *ugran_orig = new UGrandomGenerator[n_rgen];
  //     for(int i=0;i<n_rgen;i++) ugran_orig[i] = LRG.UGrandGen(i);

  //     LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form

  
  //     int pos[5];

  //     for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
  // 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
  // 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
  // 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
  // 	      for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
  // 		//do the 4D RNG
  // 		if(pos[4]==0){
  // 		  if(pos[0]>=GJP.XnodeSites()/2){
  // 		    int orig_idx = (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + stk_index_4d_off;
  // 		    int new_idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
  // 		    LRG.UGrandGen4D(new_idx) = ugran_4d_orig[orig_idx];
  // 		  }else{
  // 		    int orig_idx = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
  // 		    int new_idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
  // 		    LRG.UGrandGen4D(new_idx) = ugran_4d_orig[orig_idx];
  // 		  }
  // 		}
  // 		//do the 5D RNG
  // 		if(pos[0]>=GJP.XnodeSites()/2){
  // 		  int orig_idx = pos[4]/2*blocks_per_s_layer + (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + stk_index_5d_off;
  // 		  int new_idx = pos[4]/2*blocks_per_s_layer + pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
  // 		  LRG.UGrandGen(new_idx) = ugran_orig[orig_idx];
  // 		}else{
  // 		  int orig_idx = pos[4]/2*blocks_per_s_layer + pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
  // 		  int new_idx = pos[4]/2*blocks_per_s_layer + pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
  // 		  LRG.UGrandGen(new_idx) = ugran_orig[orig_idx];
  // 		}

  // 	      }
  // 	    }
  // 	  }
  // 	}
  //     }
  //     delete[] ugran_4d_orig;
  //     delete[] ugran_orig;
  //   }//single node

  // }//gpx and gpy
  // else if(gparity_X && gparity_Y){
  //   SingleToDouble4dRNG fourDsetup;
  //   SingleToDouble5dRNG fiveDsetup;
    
  //   LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form
    
  //   if(!UniqueID()){ printf("Setting up 4D RNG\n"); fflush(stdout); }
  //   fourDsetup.Run(gparity_X,gparity_Y);      
  //   if(!UniqueID()){ printf("Setting up 5D RNG\n"); fflush(stdout); }
  //   fiveDsetup.Run(gparity_X,gparity_Y);  
  // }















// #ifdef DOUBLE_RNG_TEST
//   if(gparity_X && !gparity_Y){
//     if(!UniqueID()) printf("Testing RNG\n");

//     //generate rands again and compare to originals
//     if(GJP.Xnodes()==1){
//       int pos[5];

//       for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
// 		LRG.AssignGenerator(pos);
// 		if(pos[4]==0){
// 		  int origflav = 0;
// 		  int origidx = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]+GJP.ZnodeSites()/2*pos[3]));
// 		  if(pos[0]>=GJP.XnodeSites()/2){
// 		    origflav = 1;
// 		    origidx = (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*pos[3]/2));
// 		  }
// 		  IFloat &origval = orig_4d_rands[origflav][origidx];
// 		  IFloat newval = LRG.Urand(FOUR_D);//do the 4D RNG
// 		  if(origval != newval){ 
// 		    printf("4D RNG disparity: (%d %d %d %d): orig %f new %f\n",pos[0],pos[1],pos[2],pos[3],origval,newval); exit(-1);
// 		  }
// 		}
// 		int origflav = 0;
// 		int origidx = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*(pos[3]/2+GJP.SnodeSites()/2*pos[4]/2)));
// 		if(pos[0]>=GJP.XnodeSites()/2){
// 		  origflav = 1;
// 		  origidx = (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*(pos[3]/2+GJP.SnodeSites()/2*pos[4]/2)));
// 		}
// 		IFloat &origval = orig_5d_rands[origflav][origidx];
// 		IFloat newval = LRG.Urand(FIVE_D);//do the 5D RNG
// 		if(origval != newval){ 
// 		  printf("5D RNG disparity: (%d %d %d %d): orig %f new %f\n",pos[0],pos[1],pos[2],pos[3],origval,newval); exit(-1);
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }

//     }else{
//       bool printnode = false;
//       if(GJP.YnodeCoor()==0 && GJP.ZnodeCoor()==0 && GJP.TnodeCoor()==0) printnode=true;

//       { //4D RNG check
// 	int n_rgen_4d = GJP.VolNodeSites()/16;// both flavours (remember volume doubled now)      
// 	int buf_size = n_rgen_4d * sizeof(Float);
// 	Float *recv_buf = (Float *) pmalloc(buf_size);
// 	Float *send_buf = (Float *) pmalloc(buf_size);
    
// 	for(int f=0;f<2;f++){
// 	  int foff = n_rgen_4d/2;
// 	  for(int site =0; site < n_rgen_4d/2; site++){
// 	    send_buf[site+f*foff] = orig_4d_rands[f][site];
// 	    //if(printnode) printf("Node %d: f %d site %d rand %f\n",GJP.XnodeCoor(),f,site,orig_4d_rands[f][site]);
// 	  }
// 	}

// 	// for(int i=0;i<n_rgen_4d;i++){
// 	// 	if(printnode) printf("Node %d: send_buf site %d = %f\n",GJP.XnodeCoor(),i,send_buf[i]);
// 	// }
      
// 	int data_nodecoor_hf1;  //what xnode coor is this nodes data for the first halves currently stored on? (second half is always on the next node)
// 	int data_nodecoor_hf2;
// 	int data_flav = 0; 
// 	int x_origin = GJP.XnodeCoor()*GJP.XnodeSites(); //x position of start of first half
// 	if(GJP.XnodeCoor()>=GJP.Xnodes()/2){  
// 	  x_origin = (GJP.XnodeCoor()-GJP.Xnodes()/2)*GJP.XnodeSites(); 
// 	  data_flav = 1;
// 	}
// 	data_nodecoor_hf1 = (x_origin/(GJP.XnodeSites()/2) ) % GJP.Xnodes();
// 	data_nodecoor_hf2 = (data_nodecoor_hf1+1) % GJP.Xnodes();

// 	Float nodes_unhappy = 1.0;
// 	Float *cur_data_buf = send_buf;
// 	Float *send_buf_p = send_buf;
// 	Float *recv_buf_p = recv_buf;
// 	int xnode_coor_of_buf_data = GJP.XnodeCoor(); 
// 	int got_hf1 = 0;
// 	int got_hf2 = 0;

// 	int pos[5]; pos[4] = 0;
// 	while(nodes_unhappy != 0.0){
// 	  if(xnode_coor_of_buf_data == data_nodecoor_hf1 || xnode_coor_of_buf_data == data_nodecoor_hf2 ){
// 	    if(xnode_coor_of_buf_data == data_nodecoor_hf1 && printnode) printf("Node %d: Buffer contains data from node %d, testing 1st half (flav %d)\n",GJP.XnodeCoor(),xnode_coor_of_buf_data,data_flav);
// 	    else if(xnode_coor_of_buf_data == data_nodecoor_hf2 && printnode) printf("Node %d: Buffer contains data from node %d, testing 2nd half (flav %d)\n",GJP.XnodeCoor(),xnode_coor_of_buf_data,data_flav);

// 	    for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 	      for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 		for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 		  for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
// 		    LRG.AssignGenerator(pos);

// 		    if(pos[0]>=GJP.XnodeSites()/2 && xnode_coor_of_buf_data == data_nodecoor_hf2){
// 		      int orig_idx = (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + data_flav * n_rgen_4d/2;

// 		      IFloat &origval = cur_data_buf[orig_idx];
// 		      IFloat newval = LRG.Urand(FOUR_D);//do the 4D RNG
// 		      if(origval != newval){ 
// 			printf("Node %d: 4D RNG disparity: (%d %d %d %d): orig %f new %f (idx %d)\n",GJP.XnodeCoor(),pos[0],pos[1],pos[2],pos[3],origval,newval,orig_idx); exit(-1);
// 		      }
// 		      got_hf2 = 1;
// 		    }else if(pos[0]<GJP.XnodeSites()/2 && xnode_coor_of_buf_data == data_nodecoor_hf1){
// 		      int orig_idx = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + data_flav * n_rgen_4d/2;

// 		      IFloat &origval = cur_data_buf[orig_idx];
// 		      IFloat newval = LRG.Urand(FOUR_D);//do the 4D RNG
// 		      if(origval != newval){ 
// 			printf("Node %d: 4D RNG disparity: (%d %d %d %d): orig %f new %f (idx %d)\n",GJP.XnodeCoor(),pos[0],pos[1],pos[2],pos[3],origval,newval,orig_idx); exit(-1);
// 		      }
// 		      got_hf1 = 1;
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	  if(got_hf1 && got_hf2) nodes_unhappy = 0.0;
// 	  else nodes_unhappy = 1.0;
	
// 	  glb_sum(&nodes_unhappy);
// 	  if(!UniqueID()) printf("nodes_unhappy = %f\n",nodes_unhappy);

// 	  if(nodes_unhappy!=0.0){
// 	    if(!UniqueID()) printf("Passing data left\n");
// 	    cur_data_buf = recv_buf_p;
// 	    getPlusData((IFloat *)recv_buf_p, (IFloat *)send_buf_p, buf_size/sizeof(Float), 0);
// 	    xnode_coor_of_buf_data = (xnode_coor_of_buf_data+1) % GJP.Xnodes();
// 	    //swap buffers over for next send
// 	    Float *tmp = recv_buf_p;
// 	    recv_buf_p = send_buf_p;
// 	    send_buf_p = tmp;

// 	    // for(int i=0;i<n_rgen_4d;i++){
// 	    //   if(printnode) printf("Node %d: post-pass recv_buf site %d = %f\n",GJP.XnodeCoor(),i,cur_data_buf[i]);
// 	    // }
// 	  }


// 	}
// 	pfree(recv_buf);
// 	pfree(send_buf);
//       }//end of 4d RNG check

//       { //5D RNG check
// 	int n_rgen = GJP.SnodeSites()/2*GJP.VolNodeSites()/16;// both flavours (remember volume doubled now)      
// 	int buf_size = n_rgen * sizeof(Float);
// 	Float *recv_buf = (Float *) pmalloc(buf_size);
// 	Float *send_buf = (Float *) pmalloc(buf_size);
	
// 	int pos[5];
// 	for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
// 	  for(int f=0;f<2;f++){
// 	    for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 	      for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 		for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 		  for(pos[0]=0;pos[0]<GJP.XnodeSites()/2;pos[0]+=2){
// 		    int osite = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*(pos[3]/2+GJP.SnodeSites()/2*pos[4]/2)));
// 		    int bsite = pos[4]/2*GJP.VolNodeSites()/16 +  pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*pos[3]/2)) + f*GJP.VolNodeSites()/32;
// 		    send_buf[bsite] = orig_5d_rands[f][osite];
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	}
      
// 	int data_nodecoor_hf1;  //what xnode coor is this nodes data for the first halves currently stored on? (second half is always on the next node)
// 	int data_nodecoor_hf2;
// 	int data_flav = 0; 
// 	int x_origin = GJP.XnodeCoor()*GJP.XnodeSites(); //x position of start of first half
// 	if(GJP.XnodeCoor()>=GJP.Xnodes()/2){  
// 	  x_origin = (GJP.XnodeCoor()-GJP.Xnodes()/2)*GJP.XnodeSites(); 
// 	  data_flav = 1;
// 	}
// 	data_nodecoor_hf1 = (x_origin/(GJP.XnodeSites()/2) ) % GJP.Xnodes();
// 	data_nodecoor_hf2 = (data_nodecoor_hf1+1) % GJP.Xnodes();

// 	Float nodes_unhappy = 1.0;
// 	Float *cur_data_buf = send_buf;
// 	Float *send_buf_p = send_buf;
// 	Float *recv_buf_p = recv_buf;
// 	int xnode_coor_of_buf_data = GJP.XnodeCoor(); 
// 	int got_hf1 = 0;
// 	int got_hf2 = 0;
	
// 	while(nodes_unhappy != 0.0){
// 	  if(xnode_coor_of_buf_data == data_nodecoor_hf1 || xnode_coor_of_buf_data == data_nodecoor_hf2 ){
// 	    if(xnode_coor_of_buf_data == data_nodecoor_hf1 && printnode) printf("Node %d: Buffer contains data from node %d, testing 1st half (flav %d)\n",GJP.XnodeCoor(),xnode_coor_of_buf_data,data_flav);
// 	    else if(xnode_coor_of_buf_data == data_nodecoor_hf2 && printnode) printf("Node %d: Buffer contains data from node %d, testing 2nd half (flav %d)\n",GJP.XnodeCoor(),xnode_coor_of_buf_data,data_flav);

// 	    for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
// 	      for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 		for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 		  for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 		    for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
// 		      LRG.AssignGenerator(pos);

// 		      if(pos[0]>=GJP.XnodeSites()/2 && xnode_coor_of_buf_data == data_nodecoor_hf2){
// 			int orig_idx = pos[4]/2 * GJP.VolNodeSites()/16 + (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + data_flav * GJP.VolNodeSites()/32;

// 			IFloat &origval = cur_data_buf[orig_idx];
// 			IFloat newval = LRG.Urand(FIVE_D);
// 			if(origval != newval){ 
// 			  printf("Node %d: 5D RNG disparity: (%d %d %d %d %d): orig %f new %f (idx %d)\n",GJP.XnodeCoor(),pos[0],pos[1],pos[2],pos[3],pos[4],origval,newval,orig_idx); exit(-1);
// 			}
// 			got_hf2 = 1;
// 		      }else if(pos[0]<GJP.XnodeSites()/2 && xnode_coor_of_buf_data == data_nodecoor_hf1){
// 			int orig_idx = pos[4]/2 * GJP.VolNodeSites()/16 + pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + data_flav * GJP.VolNodeSites()/32;

// 			IFloat &origval = cur_data_buf[orig_idx];
// 			IFloat newval = LRG.Urand(FIVE_D);
// 			if(origval != newval){ 
// 			  printf("Node %d: 5D RNG disparity: (%d %d %d %d %d): orig %f new %f (idx %d)\n",GJP.XnodeCoor(),pos[0],pos[1],pos[2],pos[3],pos[4],origval,newval,orig_idx); exit(-1);
// 			}
// 			got_hf1 = 1;
// 		      }
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	  if(got_hf1 && got_hf2) nodes_unhappy = 0.0;
// 	  else nodes_unhappy = 1.0;
	
// 	  glb_sum(&nodes_unhappy);
// 	  if(!UniqueID()) printf("nodes_unhappy = %f\n",nodes_unhappy);

// 	  if(nodes_unhappy!=0.0){
// 	    if(!UniqueID()) printf("Passing data left\n");
// 	    cur_data_buf = recv_buf_p;
// 	    getPlusData((IFloat *)recv_buf_p, (IFloat *)send_buf_p, buf_size/sizeof(Float), 0);
// 	    xnode_coor_of_buf_data = (xnode_coor_of_buf_data+1) % GJP.Xnodes();
// 	    //swap buffers over for next send
// 	    Float *tmp = recv_buf_p;
// 	    recv_buf_p = send_buf_p;
// 	    send_buf_p = tmp;

// 	    // for(int i=0;i<n_rgen_4d;i++){
// 	    //   if(printnode) printf("Node %d: post-pass recv_buf site %d = %f\n",GJP.XnodeCoor(),i,cur_data_buf[i]);
// 	    // }
// 	  }


// 	}
// 	pfree(recv_buf);
// 	pfree(send_buf);
//       }//end of 5d RNG check
//     }
//     printf("Passed RNG test\n");
//   }

// #endif






// #define DOUBLE_RNG_TEST
// #ifdef DOUBLE_RNG_TEST
//   IFloat orig_4d_rands[2][GJP.VolNodeSites()/16]; //one per RNG
//   IFloat orig_5d_rands[2][GJP.SnodeSites()/2*GJP.VolNodeSites()/16];
//   if(gparity_X && !gparity_Y){
//     //generate fields of 4D and 5D random numbers and reset the generator. After lattice doubling compare the new and old random fields to ensure they are the same
//     LatRanGen LRGbak(LRG);
//     int pos[5];

//     for(int flav=0;flav<2;flav++){
//       for(pos[4]=0;pos[4]<GJP.SnodeSites();pos[4]+=2){
// 	for(pos[3]=0;pos[3]<GJP.TnodeSites();pos[3]+=2){
// 	  for(pos[2]=0;pos[2]<GJP.ZnodeSites();pos[2]+=2){
// 	    for(pos[1]=0;pos[1]<GJP.YnodeSites();pos[1]+=2){
// 	      for(pos[0]=0;pos[0]<GJP.XnodeSites();pos[0]+=2){
// 		LRG.AssignGenerator(pos,flav);
// 		if(pos[4]==0){
// 		  int idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*pos[3]/2));
// 		  orig_4d_rands[flav][idx] = LRG.Urand(FOUR_D);//do the 4D RNG
// 		}
// 		int idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2+GJP.ZnodeSites()/2*(pos[3]/2+GJP.SnodeSites()/2*pos[4]/2)));
// 		orig_5d_rands[flav][idx] = LRG.Urand(FIVE_D);//do the 5D RNG
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//     LRG = LRGbak;
//   }
// #endif
