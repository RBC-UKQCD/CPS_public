#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/random.h>
#include<util/time_cps.h>

#include<alg/alg_hmc.h>
#include<alg/common_arg.h>
#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>

#include<alg/alg_int.h>
#include<alg/int_arg.h>

#include<alg/no_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_plaq.h>
#include<alg/alg_pbp.h>
#include<alg/pbp_arg.h>
#include<alg/alg_remez.h>
#include<alg/alg_wline.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>
#include<util/command_line.h>
//#include <sys/bgl/bgl_sys_all.h>

#undef USE_SCU_CHECKSUMS
#ifdef USE_SCU_CHECKSUMS
#include <qcdocos/scu_checksum.h>
#endif
//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

HmcArg hmc_arg;
HmcArg hmc_arg_pass;

ActionGaugeArg gauge_arg;
IntABArg ab1_arg;
DoArg do_arg;


void checkpoint(int traj);

void LRGcheck(FermionFieldDimension frm_dim = FIVE_D){
  if(frm_dim == FIVE_D)
    printf("Doing 5d LRG check\n");
  else
    printf("Doing 4d LRG check\n");

  int smax = GJP.SnodeSites();
  if(frm_dim == FOUR_D) smax = 1;

  int nstacked = 1;
  if(GJP.Gparity()) nstacked =2;

  int i5,i4;

  for(int s=0;s<smax;s++){
    for(int t=0;t<GJP.TnodeSites();t++){
      for(int z=0;z<GJP.ZnodeSites();z++){
	for(int y=0;y<GJP.YnodeSites();y++){
	  for(int x=0;x<GJP.XnodeSites();x++){  
	        
	    for(int field_idx = 0;field_idx < nstacked;field_idx++){
	      LRG.AssignGenerator(x,y,z,t,s,field_idx);
	      LRG.GetIndices(i5,i4);
	      int idx= i5;
	      if(frm_dim==FOUR_D) idx = i4;
	      //int idx=0;

	      printf("%d %d %d %d %d : field %d,   %e    %e\n",x,y,z,t,s,field_idx,LRG.Urand(frm_dim),LRG.Urand(frm_dim));
	    }

	  }
	}
      }
    }
  }

  printf("Finished LRG check\n");
  //exit(0);
}











int main(int argc, char *argv[])
{ 
  Start(&argc,&argv);
  CommandLine::is(argc,argv);

  char *cname=argv[0];
  char *fname="main()";
  Float dtime;

  bool gparity_X(false);
  bool gparity_Y(false);

  int arg0 = CommandLine::arg_as_int(0);
  printf("Arg0 is %d\n",arg0);
  if(arg0==0){
    gparity_X=true;
    printf("Doing G-parity HMC test in X direction\n");
  }else if(arg0==1){
    printf("Doing standard HMC test\n");
  }else{
    printf("Doing G-parity HMC test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
  }

  bool dbl_latt_storemode(false);
  bool load_config(false);
  char *load_config_file;
  int size[] = {4,4,4,4,4};
  bool save_lrg(false);
  char *save_lrg_file;
  bool load_lrg(false);
  char *load_lrg_file;
  
  int i=2;
  while(i<argc){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-load_config",15) == 0){
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
    }else if( strncmp(cmd,"-save_lrg",15) == 0){
      if(i==argc-1){ printf("-save_lrg requires an argument\n"); exit(-1); }
      save_lrg=true;
      save_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-load_lrg",15) == 0){
      if(i==argc-1){ printf("-load_lrg requires an argument\n"); exit(-1); }
      load_lrg=true;
      load_lrg_file = argv[i+1];
      i+=2;
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }
  SerialIO::dbl_latt_storemode = dbl_latt_storemode;

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
  if(load_config){
    do_arg.start_conf_kind = START_CONF_FILE;
    do_arg.start_conf_filename = load_config_file;
  }else{
    do_arg.start_conf_kind = START_CONF_ORD;
  }
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
  do_arg.verbose_level = -1202;//VERBOSE_DEBUG_LEVEL; //-1202;
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

  if(gparity_X) do_arg.x_bc = BND_CND_GPARITY;
  if(gparity_Y) do_arg.y_bc = BND_CND_GPARITY;

  CommonArg common_arg_hmc;
 
  // do_arg.verbose_level=VERBOSE_RESULT_LEVEL;
  GJP.Initialize(do_arg);
  // VRB.Level(VERBOSE_RESULT_LEVEL);
  LRG.Initialize();
  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }

  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }


  //LRGcheck(FOUR_D);
  //LRGcheck(FIVE_D);
  //exit(0);

  // Outer config loop

  gauge_arg.gluon = G_CLASS_IMPR_RECT;
  gauge_arg.action_arg.force_measure = FORCE_MEASURE_NO;
  gauge_arg.action_arg.force_label = "Gauge";

  hmc_arg.steps_per_traj = 1;
  hmc_arg.step_size =   1.2500000000000000e-01;
  hmc_arg.metropolis = METROPOLIS_NO;//YES;
  hmc_arg.reunitarize = REUNITARIZE_YES;
  hmc_arg.reverse = REVERSE_NO;
  hmc_arg.reproduce = REPRODUCE_NO;
  hmc_arg.reproduce_attempt_limit = 1;
  hmc_arg.wfm_md_sloppy = 1;

  hmc_arg_pass = hmc_arg;

  //!< Create fictitous Hamiltonian (mom + action)
  AlgMomentum mom;
  AlgActionGauge gauge(mom, gauge_arg);
  
  //!< Construct numerical integrators
  ab1_arg.type = INT_OMELYAN;
  ab1_arg.A_steps = 1;
  ab1_arg.B_steps = 1;
  ab1_arg.level = TOP_LEVEL_INTEGRATOR;
  ab1_arg.lambda =   2.2000000000000000e-01;
    
  AlgIntOmelyan ab1(mom, gauge, ab1_arg);
  
  //FILE *fp = fopen("links.dat","w");

  for(int conf=0; conf< 1; conf ++ ) {
    for(int traj=0;traj< 10;traj++) {
      VRB.Result("","main()","Running traj %d without reproduction\n",traj);
      hmc_arg_pass.reproduce = REPRODUCE_NO;
      
      //!< Run hybrid Monte Carlo
      AlgHmc hmc(ab1, common_arg_hmc, hmc_arg_pass);
      Float time = -dclock();
      hmc.run();
      time += dclock();
      print_flops("AlgHmc","run()",0,time);

      {
	char lat_file[256];
	// Save this config to disk
	Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	
	sprintf(lat_file,"gauge.%d.traj%d",conf,traj);
	QioArg wt_arg(lat_file,0.001);
	
	wt_arg.ConcurIONumber=32;
	WriteLatticeParallel wl;
	wl.setHeader("ens_idx","ens_label",conf);
	wl.write(lat,wt_arg);
	
	if(!wl.good()) 
	  ERR.General(cname,fname,"Failed write lattice %s",lat_file);
	LatticeFactory::Destroy();
      }


      if(GJP.Gparity()){
	Lattice &lattice = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	int base_loc[4];
	for(int i=0;i<4;i++) base_loc[i] = GJP.NodeCoor(i)*GJP.NodeSites(i); //absolute pos of start point
    
	//printf("Lattice StrOrd is %s\n",StrOrdType_map[(int)lattice.StrOrd()].name);

	int pos[4];
	
	for(int t=0;t<GJP.NodeSites(3);t++){
	  pos[3] = t;
	  for(int z=0;z<GJP.NodeSites(2);z++){
	    pos[2] = z;
	    for(int y=0;y<GJP.NodeSites(1);y++){
	      pos[1] = y;
	      for(int x=0;x<GJP.NodeSites(0);x++){
		pos[0] = x;
		
		for(int mu=0;mu<4;mu++){
		  Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu));
		  Matrix *Ustar = const_cast<Matrix *>(lattice.GetLink(pos, mu,1));

		  Matrix tmp;
		  tmp.Conj((IFloat*)U);
		    
		  for(int ii=0;ii<18;ii++){
		    if(tmp.elem(ii)!=Ustar->elem(ii) ){
		      printf("Symmetry broken on traj idx %d at location %d %d %d %d, %d (elem %d)\n",traj,x,y,z,t,mu,ii); 
		      printf("%.16f != %.16f\n",tmp.elem(ii),Ustar->elem(ii));
		      printf("Matrix should be:\n");
		      for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
			  const Complex &c = tmp(i,j);
			  printf("[%.3f %.3f] ",c.real(), c.imag());
			}
			printf("\n");
		      }
		      printf("\n");
		      
		      printf("But got:\n");
		      for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
			  const Complex &c = (*Ustar)(i,j);
			  printf("[%.3f %.3f] ",c.real(), c.imag());
			}
			printf("\n");
		      }
		      printf("\n");
			
		      exit(-1);
		    }
		  }
		  
		}
		
	      }
	      
	      
	    }
	  }
	}
	
      // 	fprintf(fp,"Traj idx %d\n",traj);
      // 	fprintf(fp,"Gauge links U\n");

      // 	for(int t=0;t<GJP.NodeSites(3);t++){
      // 	  pos[3] = t;
      // 	  for(int z=0;z<GJP.NodeSites(2);z++){
      // 	    pos[2] = z;
      // 	    for(int y=0;y<GJP.NodeSites(1);y++){
      // 	      pos[1] = y;
      // 	      for(int x=0;x<GJP.NodeSites(0);x++){
      // 		pos[0] = x;
		
      // 		for(int mu=0;mu<4;mu++){
      // 		  Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu));
      // 		  fprintf(fp,"%d %d %d %d; %d:\n",x+base_loc[0],y+base_loc[1],z+base_loc[2],t+base_loc[3],mu);
      // 		  for(int i=0;i<3;i++){
      // 		    for(int j=0;j<3;j++){
      // 		      const Complex &c = (*U)(i,j);
      // 		      fprintf(fp,"[%.3f %.3f] ",c.real(), c.imag());
      // 		    }
      // 		    fprintf(fp,"\n");
      // 		  }
      // 		  fprintf(fp,"\n");
      // 		}
      // 	      }
      // 	    }
      // 	  }
      // 	}
      	
      // 	fprintf(fp,"Conjugated gauge links U*\n");

      // 	for(int t=0;t<GJP.NodeSites(3);t++){
      // 	  pos[3] = t;
      // 	  for(int z=0;z<GJP.NodeSites(2);z++){
      // 	    pos[2] = z;
      // 	    for(int y=0;y<GJP.NodeSites(1);y++){
      // 	      pos[1] = y;
      // 	      for(int x=0;x<GJP.NodeSites(0);x++){
      // 		pos[0] = x;
		
      // 		for(int mu=0;mu<4;mu++){
      // 		  Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu,1));
      // 		  fprintf(fp,"%d %d %d %d; %d:\n",x+base_loc[0],y+base_loc[1],z+base_loc[2],t+base_loc[3],mu);
      // 		  for(int i=0;i<3;i++){
      // 		    for(int j=0;j<3;j++){
      // 		      const Complex &c = (*U)(i,j);
      // 		      fprintf(fp,"[%.3f %.3f] ",c.real(), c.imag());
      // 		    }
      // 		    fprintf(fp,"\n");
      // 		  }
      // 		  fprintf(fp,"\n");
      // 		}
      // 	      }
      // 	    }
      // 	  }
      // 	}
       	LatticeFactory::Destroy();
      }
	      
    }//End of inter-cfg sweep

    char lat_file[256];
    // Save this config to disk
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    
    sprintf(lat_file,"gauge.%d",conf);
    QioArg wt_arg(lat_file,0.001);
    
    wt_arg.ConcurIONumber=32;
    WriteLatticeParallel wl;
    wl.setHeader("ens_idx","ens_label",conf);
    wl.write(lat,wt_arg);
    
    if(!wl.good()) 
      ERR.General(cname,fname,"Failed write lattice %s",lat_file);
    LatticeFactory::Destroy();


  } //End config loop
  //lose(fp);
  //AlgIntAB::Destroy(ab1);

  End();
  printf("End of program\n");
 return(0);
}

void checkpoint(int traj)
{
  char *cname="cps::";
  char *fname="checkpoint()";


}
