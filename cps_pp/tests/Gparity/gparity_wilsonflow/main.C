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
#include <alg/alg_wilsonflow.h>
#include <alg/alg_actiondensity.h>

#include <util/data_shift.h>
#include <util/gparity_singletodouble.h>

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
    printf("Doing G-parity Wilson Flow test in X direction\n");
  }else{
    printf("Doing G-parity Wilson Flow test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
  }

  bool smartrun(true);
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
    }else if( strncmp(cmd,"-wf_no_smartrun",15) == 0){
      smartrun = false;
      i++;
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

  //2f section

  //Backup original lattice
  int array_size = 2*lattice->GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
  Matrix *orig_lattice = (Matrix *) pmalloc(array_size);
  memcpy((void*)orig_lattice, (void*)lattice->GaugeField(), array_size);
  
  //Run the Wilson Flow code
  Matrix *trans_lat_2f = (Matrix *) pmalloc(array_size);
  Float action_density_2f;
  {
    CommonArg c;
    AlgWilsonFlow awf(*lattice,&c,0.05);
    AlgActionDensity ad(*lattice,&c);

    if(smartrun){
      awf.smartrun();
      ad.smartrun(&action_density_2f);
    }else{
      awf.run();
      ad.run(&action_density_2f);
    }

    memcpy((void*)trans_lat_2f, (void*)lattice->GaugeField(), array_size);
  }

  //Start preparing for 1f section
  lattice->FreeGauge(); //free memory and reset
  delete lattice; //lattice objects are singleton (scope_lock)

  //1f section
  if(UniqueID()==0) printf("Starting double lattice inversion\n");  
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
    QDP::Layout::setLattSize(nrow);
    QDP::Layout::create();
  }
#endif

  GwilsonFdwf doubled_lattice;
  setup_double_latt(doubled_lattice,orig_lattice,gparity_X,gparity_Y);
  setup_double_rng(gparity_X,gparity_Y);
 
  //Also bring the transformed 2f lattice into 1f format
  array_size = doubled_lattice.GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
  Matrix *trans_lat_2f_1fform = (Matrix *) pmalloc(array_size);
  setup_double_matrixfield(trans_lat_2f_1fform , trans_lat_2f, 4,gparity_X, gparity_Y);
  pfree(trans_lat_2f);

  Float action_density_1f;

  //Do the Wilson Flow on the 1f lattice
  {
    CommonArg c;
    AlgWilsonFlow awf(doubled_lattice,&c,0.05);
    AlgActionDensity ad(doubled_lattice,&c);

    if(smartrun){
      awf.smartrun();
      ad.smartrun(&action_density_1f);
    }else{
      awf.run();
      ad.run(&action_density_1f);
    }

  }
  
  //Compare 1f and 2f transformed lattices
  bool fail = false;
  for(int i=0;i< 4*18*GJP.VolNodeSites(); i++){
    Float* _1f = (Float*)doubled_lattice.GaugeField() + i;
    Float* _2f = (Float*)trans_lat_2f_1fform + i;

    int rem = i;
    int screim = rem%18; rem/=18;
    int mu = rem%4;   rem/=4;
    int x[4];
    for(int d=0;d<4;d++){
      x[d] = rem % GJP.NodeSites(d) + GJP.NodeCoor(d)*GJP.NodeSites(d);  rem/=GJP.NodeSites(d);
    }

    if( fabs(*_1f-*_2f) > 1e-10 ){
      printf("!Wilson Flow comparison fail screim %d, mu %d, coor %d %d %d %d: %.12le %.12le\n",screim,mu,x[0],x[1],x[2],x[3],*_1f,*_2f); fail=true;
    }// else{
    //   printf("Wilson Flow comparison pass screim %d, mu %d, coor %d %d %d %d: %.12le %.12le\n",screim,mu,x[0],x[1],x[2],x[3],*_1f,*_2f);
    // }

  }
  if(fail){
    printf("Failed Wilson Flow test\n"); exit(-1); 
  }else{
    printf("Passed Wilson Flow test\n");
  }

  //Compare action densities

  if( fabs(action_density_1f - action_density_2f) > 1e-10){
    printf("!Action Density comparison fail: %.12le %.12le\n",action_density_1f,action_density_2f); fail=true;
  }
  if(fail){
    printf("Failed Action Density test\n"); exit(-1); 
  }else{
    printf("Passed Action Density test\n");
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
