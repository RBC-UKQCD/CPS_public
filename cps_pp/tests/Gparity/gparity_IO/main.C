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

#include <util/data_shift.h>

using namespace std;
USING_NAMESPACE_CPS

int main(int argc,char *argv[])
{
  Start(&argc,&argv);
  CommandLine::is(argc,argv);
  
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

  bool dbl_latt_storemode(false);//writes out RNG and gauge then exits, no test performed
  int size[] = {4,4,4,4,4};

  int i=2;
  while(i<argc){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-latt",10) == 0){
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
      //if gparity in X and Y, it stores a quadruple lattice here
      dbl_latt_storemode = true;
      i++;
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
  do_arg.verbose_level =-1202;//VERBOSE_DEBUG_LEVEL; //-1202;
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

  GJP.Initialize(do_arg);

  if(UniqueID()==0) printf("Initialising RNG\n");
  LRG.Initialize();

  if(!(LRG == LRG)){
    printf("LRG == fail!\n");
    exit(-1);
  }

  if(UniqueID()==0) printf("Creating gauge field\n");
  GwilsonFdwf lattice;
  lattice.SetGfieldDisOrd(); //unit gauge

  unsigned int cksum = lattice.CheckSum();
  if(UniqueID()==0) printf("Initial lattice checksum: %u\n",cksum);

  for(int parallelread = 0; parallelread < 2; parallelread++){
    for(int parallelwrite = 0; parallelwrite < 2; parallelwrite++){
      if(UniqueID()==0){
	if(!parallelread && !parallelwrite) printf("Starting serial read and serial write\n");
	else if(!parallelread && parallelwrite) printf("Starting serial read and parallel write\n");
	else if(parallelread && !parallelwrite) printf("Starting parallel read and serial write\n");
	else printf("Starting parallel read and parallel write\n");
      }

      //Write the RNG state
      char *rng_file = "rng.dat";
      if(UniqueID()==0) printf("Writing RNG state to %s\n",rng_file);
      {
	if(parallelwrite) LRG.setParallel();
	else LRG.setSerial();
			    
	if(LRG.parIO()) printf("Using parallel IO\n");
	else printf("Using serial IO\n");
	LRG.Write(rng_file,32);
      }
  
      //Write the gauge field
      char *lat_file = "gauge.dat";
      if(UniqueID()==0) printf("Writing gauge field to %s\n",lat_file);
      {
	QioArg wt_arg("gauge.dat",0.001);
    
	wt_arg.ConcurIONumber=32;
	WriteLatticeParallel wl;
	if(parallelwrite) wl.setParallel();
	else wl.setSerial();

	wl.setHeader("disord_id","disord_label",0);
	wl.write(lattice,wt_arg);
    
	if(!wl.good()) ERR.General("main","()","Failed write lattice %s",lat_file);

	if(UniqueID()==0) printf("Config written.\n");
      }
  
      if(dbl_latt_storemode){
	if(UniqueID()==0) printf("Double latt output written, ending program.\n");
	exit(0);
      }

      //Save a copy of the current RNG
      LatRanGen LRGcopy(LRG);
  
      //Load the RNG from file
      if(UniqueID()==0) printf("Loading RNG state from %s\n",rng_file);
      {
	if(parallelread) LRG.setParallel();
	else LRG.setSerial();

	LRG.Read(rng_file,32);
      }
      if(!(LRG == LRGcopy)){
	printf("RNG state load/save fail\n"); exit(-1);
      }
      if(UniqueID()==0) printf("RNG load successful\n");

      //Save a copy of the current gauge field
      Matrix gaugecopy[4*GJP.VolNodeSites()];
      Matrix *gauge = lattice.GaugeField();
      for(int i=0;i<4*GJP.VolNodeSites();i++){
	gaugecopy[i] = gauge[i];
      }


      if(UniqueID()==0) printf("Reloading gauge field\n");
      {
	ReadLatticeParallel readLat;
	if(parallelread) readLat.setParallel();
	else readLat.setSerial();

	if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",lat_file);
	readLat.read(lattice,lat_file);
	if(UniqueID()==0) printf("Config read.\n");
      }
  
      for(int i=0;i<4*GJP.VolNodeSites();i++){
	for(int j=0;j<18;j++){
	  if( fabs(gaugecopy[i].elem(j) - gauge[i].elem(j) ) > 1e-08 ){
	    int mu = i%4; int rest = i/4;
	    int x = rest%GJP.XnodeSites(); rest/=GJP.XnodeSites();
	    int y = rest%GJP.YnodeSites(); rest/=GJP.YnodeSites();
	    int z = rest%GJP.ZnodeSites(); rest/=GJP.ZnodeSites();
	    int t = rest%GJP.TnodeSites();

	    printf("Gauge load/save fail (%d,%d,%d,%d; %d): orig %f loaded %f\n",x,y,z,t,mu,gaugecopy[i].elem(j), gauge[i].elem(j)); exit(-1);
	  }
	}
      }
      if(UniqueID()==0) printf("Gauge load successful\n");

      if(UniqueID()==0){
	if(!parallelread && !parallelwrite) printf("Serial read and serial write successful\n");
	else if(!parallelread && parallelwrite) printf("Serial read and parallel write successful\n");
	else if(parallelread && !parallelwrite) printf("Parallel read and serial write successful\n");
	else printf("Parallel read and parallel write sucessful\n");
      }
    }
  }

  if(UniqueID()==0){
    printf("IO check complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}



//  int volnodesites = GJP.VolNodeSites();


  // printf("U:\n");
  // for(int i=0;i<4*volnodesites;i++){
  //   printf("i=%d ",i);
  //   for(int j=0;j<18;j++){
  //     if(gauge[i].elem(j)>=0)
  // 	printf("%d: %.3f ",j,gauge[i].elem(j));
  //     else
  // 	printf("%d:%.3f ",j,gauge[i].elem(j));
  //   }
  //   printf("\n");
  // }
  // if(GJP.Gparity()){
  //   printf("U*:\n");
  //   for(int i=4*volnodesites;i<8*volnodesites;i++){
  //     printf("i=%d ",i-4*volnodesites);
  //     for(int j=0;j<18;j++){
	
  // 	if(j%2==0 && gauge[i].elem(j) != gauge[i-4*volnodesites].elem(j)){
  // 	  printf("Error re %d,%d : %.3f %.3f\n",i,j,gauge[i].elem(j),gauge[i-4*volnodesites].elem(j));
  // 	}else if(j%2==1 && gauge[i].elem(j) != gauge[i-4*volnodesites].elem(j)){
  // 	  printf("Error im %d,%d : %.3f %.3f\n",i,j,gauge[i].elem(j),gauge[i-4*volnodesites].elem(j));
  // 	}

  // 	if(gauge[i].elem(j)>=0)
  // 	  printf("%d: %.3f ",j,gauge[i].elem(j));
  // 	else
  // 	  printf("%d:%.3f ",j,gauge[i].elem(j));
  //     }
  //     printf("\n");
  //   }
  // }














  // //exit(0);
  
  // Matrix *gauge = lattice.GaugeField();
  // int volnodesites = GJP.VolNodeSites();
  
  // int pos[4];

  // for(int t=0;t<GJP.NodeSites(3);t++){
  //   pos[3] = t;
  //   for(int z=0;z<GJP.NodeSites(2);z++){
  //     pos[2] = z;
  //     for(int y=0;y<GJP.NodeSites(1);y++){
  // 	pos[1] = y;
  // 	for(int x=0;x<GJP.NodeSites(0);x++){
  // 	  pos[0] = x;

  // 	  for(int mu=0;mu<4;mu++){
  // 	    Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu));
  // 	    Matrix *Uconj = const_cast<Matrix *>(lattice.GetLink(pos, mu,1));

  // 	    for(int i=0;i<3;i++){
  // 	    	for(int j=0;j<3;j++){
  // 	    	  const Complex &c = (*U)(i,j);
  // 		  const Complex &cstar = (*Uconj)(i,j);
  // 		  if(c.real() != cstar.real() || c.imag() != -cstar.imag() ){
  // 		    ERR.General("main","()","Non-conj problem at %d %d %d %d; %d :  %d %d\n",x,y,z,t,mu,i,j);
  // 		  }
  // 	    	}
  // 	    }
  // 	  }
	    

  // 	}
  //     }
  //   }
  // }


  // // printf("U:\n");
  // // for(int i=0;i<4*volnodesites;i++){
  // //   printf("i=%d ",i);
  // //   for(int j=0;j<18;j++){
  // //     if(gauge[i].elem(j)>=0)
  // // 	printf("%d: %.3f ",j,gauge[i].elem(j));
  // //     else
  // // 	printf("%d:%.3f ",j,gauge[i].elem(j));
  // //   }
  // //   printf("\n");
  // // }
  // // if(GJP.Gparity()){
  // //   printf("U*:\n");
  // //   for(int i=4*volnodesites;i<8*volnodesites;i++){
  // //     printf("i=%d ",i-4*volnodesites);
  // //     for(int j=0;j<18;j++){
	
  // // 	if(j%2==0 && gauge[i].elem(j) != gauge[i-4*volnodesites].elem(j)){
  // // 	  printf("Error re %d,%d : %.3f %.3f\n",i,j,gauge[i].elem(j),gauge[i-4*volnodesites].elem(j));
  // // 	}else if(j%2==1 && gauge[i].elem(j) != gauge[i-4*volnodesites].elem(j)){
  // // 	  printf("Error im %d,%d : %.3f %.3f\n",i,j,gauge[i].elem(j),gauge[i-4*volnodesites].elem(j));
  // // 	}

  // // 	if(gauge[i].elem(j)>=0)
  // // 	  printf("%d: %.3f ",j,gauge[i].elem(j));
  // // 	else
  // // 	  printf("%d:%.3f ",j,gauge[i].elem(j));
  // //     }
  // //     printf("\n");
  // //   }
  // // }
  
  // // Matrix orig[8*volnodesites];
  // // for(int i=0;i<8*volnodesites;i++) orig[i] = gauge[i];

  // // lattice.Convert(WILSON);
  // // lattice.Convert(CANONICAL);

  // // gauge = lattice.GaugeField();
  
  // // for(int i=0;i<8*volnodesites;i++){
  // //   bool linepass =true;
  // //   for(int j=0;j<18;j++){
  // //     if(gauge[i].elem(j) != orig[i].elem(j)){
  // // 	printf("Err %d %d; %e %e\n",i,j,gauge[i].elem(j),orig[i].elem(j));
  // // 	linepass = false;
  // //     }
  // //   }
  // //   if(linepass) printf("Matrix %d pass\n",i);
  // // }
  

  // // exit(0);
  // // printf("Post convert\n");

  // // printf("cb0:\n");
  // // int cbsize = 2*volnodesites;

  // // for(int i=0;i<cbsize;i++){
  // //   printf("U  i=%d ",i);
  // //   for(int j=0;j<18;j++){
  // //     if(gauge[i].elem(j)>=0)
  // // 	printf("%d: %.3f ",j,gauge[i].elem(j));
  // //     else
  // // 	printf("%d:%.3f ",j,gauge[i].elem(j));
  // //   }
  // //   printf("\n");

  // //   if(GJP.Gparity()){
  // //     printf("U* i=%d ",i+cbsize);
  // //     for(int j=0;j<18;j++){
  // // 	if(gauge[i+cbsize].elem(j)>=0)
  // // 	  printf("%d: %.3f ",j,gauge[i+cbsize].elem(j));
  // // 	else
  // // 	  printf("%d:%.3f ",j,gauge[i+cbsize].elem(j));
  // //     }
  // //     printf("\n\n");

  // //   }


  // // }
  // // printf("cb1:\n");
  // // int offset = cbsize;
  // // if(GJP.Gparity()) offset*=2;

  // // for(int i=offset;i<offset+cbsize;i++){
  // //   printf("U  i=%d ",i);
  // //   for(int j=0;j<18;j++){
  // //     if(gauge[i].elem(j)>=0)
  // // 	printf("%d: %.3f ",j,gauge[i].elem(j));
  // //     else
  // // 	printf("%d:%.3f ",j,gauge[i].elem(j));
  // //   }
  // //   printf("\n");

  // //   if(GJP.Gparity()){
  // //     printf("U* i=%d ",i+cbsize);
  // //     for(int j=0;j<18;j++){
  // // 	if(gauge[i+cbsize].elem(j)>=0)
  // // 	  printf("%d: %.3f ",j,gauge[i+cbsize].elem(j));
  // // 	else
  // // 	  printf("%d:%.3f ",j,gauge[i+cbsize].elem(j));
  // //     }
  // //     printf("\n\n");

  // //   }


  // // }


  // // printf("Gauge field:\n\n");

  // // int base_loc[4];
  // // for(int i=0;i<4;i++) base_loc[i] = GJP.NodeCoor(i)*GJP.NodeSites(i); //absolute pos of start point
    
  // // int pos[4];

  // // for(int t=0;t<GJP.NodeSites(3);t++){
  // //   pos[3] = t;
  // //   for(int z=0;z<GJP.NodeSites(2);z++){
  // //     pos[2] = z;
  // //     for(int y=0;y<GJP.NodeSites(1);y++){
  // // 	pos[1] = y;
  // // 	for(int x=0;x<GJP.NodeSites(0);x++){
  // // 	  pos[0] = x;

  // // 	  for(int mu=0;mu<4;mu++){
  // // 	    Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu));
  // // 	    printf("%d %d %d %d; %d:\n",x+base_loc[0],y+base_loc[1],z+base_loc[2],t+base_loc[3],mu);
  // // 	    for(int i=0;i<3;i++){
  // // 	    	for(int j=0;j<3;j++){
  // // 	    	  const Complex &c = (*U)(i,j);
  // // 	    	  printf("[%.3f %.3f] ",c.real(), c.imag());
  // // 	    	}
  // // 	    	printf("\n");
  // // 	    }
  // // 	    printf("\n");
  // // 	  }
	    

  // // 	}
  // //     }
  // //   }
  // // }

  // // if(GJP.Gparity()){
  // //   printf("Conjugated Gauge field:\n\n");

  // //   for(int t=0;t<GJP.NodeSites(3);t++){
  // //     pos[3] = t;
  // //     for(int z=0;z<GJP.NodeSites(2);z++){
  // // 	pos[2] = z;
  // // 	for(int y=0;y<GJP.NodeSites(1);y++){
  // // 	  pos[1] = y;
  // // 	  for(int x=0;x<GJP.NodeSites(0);x++){
  // // 	    pos[0] = x;
	    
  // // 	    for(int mu=0;mu<4;mu++){
  // // 	      Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu, 1));
  // // 	      printf("%d %d %d %d; %d:\n",x+base_loc[0],y+base_loc[1],z+base_loc[2],t+base_loc[3],mu);
  // // 	      for(int i=0;i<3;i++){
  // // 	    	for(int j=0;j<3;j++){
  // // 	    	  const Complex &c = (*U)(i,j);
  // // 	    	  printf("[%.3f %.3f] ",c.real(), c.imag());
  // // 	    	}
  // // 	    	printf("\n");
  // // 	      }
  // // 	      printf("\n");
  // // 	    }
	    

  // // 	  }
  // // 	}
  // //     }
  // //   }

  // // }







  // // exit(0);

  // CgArg cg;
  // cg.mass =  0.5;
  // cg.max_num_iter = 5000;
  // cg.stop_rsd =   1.0000000000000000e-06;
  // cg.true_rsd =   1.0000000000000000e-06;
  // cg.RitzMatOper = NONE;
  // cg.Inverter = CG;
  // cg.bicgstab_n = 0;


  // QPropWArg qpropw_arg;
  // qpropw_arg.cg = cg;
  // qpropw_arg.x = 0;
  // qpropw_arg.y = 0;
  // qpropw_arg.z = 0;
  // qpropw_arg.t = 0;
  // if(gparity){
  //   qpropw_arg.flavor = 0; //point source on d field
  // }
  // qpropw_arg.ensemble_label = "ens";
  // qpropw_arg.ensemble_id = "ens_id";
  // qpropw_arg.StartSrcSpin = 0;
  // qpropw_arg.EndSrcSpin = 4;
  // qpropw_arg.StartSrcColor = 0;
  // qpropw_arg.EndSrcColor = 3;

  // qpropw_arg.save_prop = 1;
  // qpropw_arg.file = "prop.dat";
  
  // CommonArg common_arg("label","filename");
  
  // QPropW* src = new QPropWPointSrc(lattice, &qpropw_arg, &common_arg);
  
  // if(UniqueID()==0) printf("Inversion finished\n");

  // //int id = UniqueID();
  // //char propfile[50];
  // //sprintf(&propfile,"prop.$d.dat",id);

  // FILE * ftxt = Fopen(ADD_ID,"prop_txt","w");

  // //int pos[4];
  
  // int shift_x = GJP.XnodeCoor()*GJP.XnodeSites();
  // int shift_y = GJP.YnodeCoor()*GJP.YnodeSites();
  // int shift_z = GJP.ZnodeCoor()*GJP.ZnodeSites();
  // int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  // //Local lattice dimensions:
  // int size_x = GJP.XnodeSites();
  // int size_y = GJP.YnodeSites();
  // int size_z = GJP.ZnodeSites();
  // int size_t = GJP.TnodeSites();
  // int size_xy = size_x*size_y;
  // int vol = (GJP.VolNodeSites()/GJP.TnodeSites()); // =size_x*size_y_size_z
  // //Global lattice dimensions
  // int Size_X = GJP.Xnodes()*GJP.XnodeSites();
  // int Size_Y = GJP.Ynodes()*GJP.YnodeSites();
  // int Size_Z = GJP.Znodes()*GJP.ZnodeSites();
  // int Size_T = GJP.Tnodes()*GJP.TnodeSites();


  // for (int i=0; i<GJP.VolNodeSites(); i++) {
  //   //Global coordinates
  //   pos[3] = i/vol + shift_t;
  //   pos[2] = (i%vol)/size_xy + shift_z;
  //   pos[1] = (i%size_xy)/size_x + shift_y;
  //   pos[0] = i%size_x + shift_x;
    
  //   //printf("site (%d,%d,%d,%d)\n",pos[0],pos[1],pos[2],pos[3]);  
  //   Fprintf(ADD_ID,ftxt,"site (%d,%d,%d,%d)\n",pos[0],pos[1],pos[2],pos[3]);
    
  //   WilsonMatrix &m = (*src)[i];
  //   for(int s1=0;s1<4;s1++){
  //     for(int c1=0;c1<3;c1++){
  // 	for(int s2=0;s2<4;s2++){
  // 	  for(int c2=0;c2<3;c2++){
  // 	    Complex &val = m(s1,c1,s2,c2);
  // 	    IFloat re = val.real();
  // 	    //if(re!=0.0){
  // 	      //printf("u %d %d %d %d : %f %f\n",s1,c1,s2,c2,val.real(),val.imag());
  // 	      Fprintf(ADD_ID,ftxt,"u %d %d %d %d : %f %f\n",s1,c1,s2,c2,val.real(),val.imag());
  // 	      //}
  // 	  }
  // 	}
  //     }
  //   }
  // }
  // if(GJP.Gparity()){
  //   for (int i=0; i<GJP.VolNodeSites(); i++) {
  //     pos[3] = i/vol + shift_t;
  //     pos[2] = (i%vol)/size_xy + shift_z;
  //     pos[1] = (i%size_xy)/size_x + shift_y;
  //     pos[0] = i%size_x + shift_x;

  //     //printf("site (%d,%d,%d,%d)\n",pos[0],pos[1],pos[2],pos[3]);  
  //     Fprintf(ADD_ID,ftxt,"site (%d,%d,%d,%d)\n",pos[0],pos[1],pos[2],pos[3]);

  //     WilsonMatrix &m = (*src)[i+GJP.VolNodeSites()];
  //     for(int s1=0;s1<4;s1++){
  // 	for(int c1=0;c1<3;c1++){
  // 	  for(int s2=0;s2<4;s2++){
  // 	    for(int c2=0;c2<3;c2++){
  // 	      Complex &val = m(s1,c1,s2,c2);
  // 	      IFloat re = val.real();
  // 	      //if(re!=0.0){
  // 		//printf("d %d %d %d %d : %f %f\n",s1,c1,s2,c2,val.real(),val.imag());
  // 		Fprintf(ADD_ID,ftxt,"d %d %d %d %d : %f %f\n",s1,c1,s2,c2,val.real(),val.imag());
  // 		//}
  // 	    }
  // 	  }
  // 	}
  //     }
  //   }
  // }

  // Fclose(ADD_ID,ftxt);

