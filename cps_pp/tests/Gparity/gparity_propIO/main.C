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

#include <util/data_shift.h>

#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>

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

  int size[] = {4,4,4,4,4};

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

  if(verbose) do_arg.verbose_level = VERBOSE_DEBUG_LEVEL;

  if(gparity_X) do_arg.x_bc = BND_CND_GPARITY;
  if(gparity_Y) do_arg.y_bc = BND_CND_GPARITY;

  GJP.Initialize(do_arg);

  SerialIO::dbl_latt_storemode = dbl_latt_storemode;
  
  LRG.Initialize();

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }


  GwilsonFdwf lattice;
					       
  if(!load_config){
    printf("Creating gauge field\n");
    lattice.SetGfieldDisOrd(); //unit gauge
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

  printf("Beginning decode\n");
  
  JobPropagatorArgs prop_args;
  if(!prop_args.Decode("parg.vml","prop_arg")){
    printf("Failed to decode parg.vml\n"); exit(-1);
  }
  printf("Decode finished\n");

  if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args.props.props_len);
  
  PropManager::setup(prop_args);

  for(int i=0;i<prop_args.props.props_len;i++){
    PropagatorArg &arg = prop_args.props.props_val[i];
    PropagatorContainer &prop = PropManager::getProp(arg.generics.tag);
    PropIOAttrArg *ioarg;
    if(!prop.getAttr(ioarg)) continue;

    if(UniqueID()) printf("Propagator %s is on disk in file %s\n",prop_args.props.props_val[i].generics.tag,ioarg->qio_filename);

    //create a new prop args
    PropagatorArg new_parg;
    new_parg.deep_copy(arg);
    delete new_parg.generics.tag;
    char *newtag = new char[strlen(arg.generics.tag)+10];
    sprintf(newtag,"%s load",arg.generics.tag);
    new_parg.generics.tag = newtag;

    for(int x=0;x<new_parg.attributes.attributes_len;x++){
      if(new_parg.attributes.attributes_val[x].type == PROP_IO_ATTR){
  	new_parg.attributes.attributes_val[x].AttributeContainer_u.prop_io_attr.prop_on_disk = true;
      }
    }  
    
    PropManager::addProp(new_parg);

    PropagatorContainer &prop_loaded = PropManager::getProp(new_parg.generics.tag);

    if(!UniqueID()) printf("Inverting propagator %s\n",arg.generics.tag);
    QPropW &prop_q = prop.getProp(lattice);

    if(!UniqueID()) printf("Reloading propagator %s\n",arg.generics.tag);
    QPropW &read_prop_q = prop_loaded.getProp(lattice);
    
    if(!UniqueID()) printf("Props %p %p\n",&prop_q,&read_prop_q);

    Float errCnt(0.);
    Float sumerr(0.);

    int nstacked = 1;
    if(GJP.Gparity()) nstacked = 2;
    int stkoff = GJP.VolNodeSites();

    for(int stk=0;stk<nstacked;stk++){
      for(int index(0); index < GJP.VolNodeSites(); ++index){
       
	WilsonMatrix mat_read, mat_calc;
       
	mat_calc = prop_q[index + stk*stkoff];
	mat_read = read_prop_q[index + stk*stkoff];
       
	for(int s_src(0); s_src < 4; ++s_src)
	  for(int c_src(0); c_src < 3; ++c_src)
	    for(int s_snk(0); s_snk < 4; ++s_snk)
	      for(int c_snk(0); c_snk < 3; ++c_snk){
	       
		Complex tmp_calc = mat_calc(s_snk,c_snk,s_src,c_src);
		Complex tmp_read = mat_read(s_snk,c_snk,s_src,c_src);
	       
		Float diff;
		diff = ( fabs(tmp_calc.real()-tmp_read.real()) + fabs(tmp_calc.imag()-tmp_read.imag()) )/ sqrt((tmp_calc.real()*tmp_calc.real() + tmp_calc.imag()*tmp_calc.imag()));
	       
		if( diff > 1e-12){
		 
		  errCnt += 1.0;
		  sumerr+=diff;
		  if(GJP.Gparity()){
		    VRB.Result("main","()","mismatch propagator: stacked flavor idx %d index %i snk %i %i src %i %i\n %f: (%f,%f) <-> (%f,%f)\n",
			       stk,index, s_snk,c_snk,s_src,c_src,
			       diff, tmp_calc.real(), tmp_calc.imag(), tmp_read.real(), tmp_read.imag() );
		  }else{
		    VRB.Result("main","()","mismatch propagator: index %i snk %i %i src %i %i\n %f: (%f,%f) <-> (%f,%f)\n",
			       index, s_snk,c_snk,s_src,c_src,
			       diff, tmp_calc.real(), tmp_calc.imag(), tmp_read.real(), tmp_read.imag() );
		  }
		}
	       
	      }
      }
    }
    
    glb_sum_five(&errCnt);
    glb_sum_five(&sumerr);
    Float averr=sumerr/errCnt;
     
    if( fabs(errCnt) > 0.){
      VRB.Result("main","()"," Reload prop. with TOTAL NUMBER OF ERRORS: %f\n",errCnt);
      VRB.Result("main","()"," Average error: %e\n",averr);
      VRB.Result("main","()"," The precision is set at: %e\n",1e-10);
    } else {
      VRB.Result("main","()"," Reload prop. successfully!\n");
      VRB.Result("main","()"," The precision is set at: %e\n",1e-10);
    }
    

  }

  //GparityContractArg contract_args;
  //if(!contract_args.Decode("contractarg.vml","contract_arg")){
  //  printf("Failed to decode contractarg.vml\n"); exit(-1);
  //}

  //CommonArg carg("label","filename");
  //AlgGparityContract contract(lattice,carg,contract_args);
  //contract.run(0);

  if(gauge_fix) lattice.FixGaugeFree();

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}




  // prop_args.props.props_len = 1;
  // prop_args.props.props_val = new PropagatorArg[1];

  // prop_args.props.props_val[0].attributes.attributes_len = 2;
  // prop_args.props.props_val[0].attributes.attributes_val = new AttributeSet[2];
  // prop_args.props.props_val[0].attributes.attributes_val[0].type = POINT_SOURCE_ATTR;
  // prop_args.props.props_val[0].attributes.attributes_val[1].type = MOMENTUM_ATTR;

  // prop_args.Encode("parg.vml","prop_arg");



  
  // Matrix *gauge = lattice.GaugeField();
  // int volnodesites = GJP.VolNodeSites();
  
  //int pos[4];

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
  
  // Matrix orig[8*volnodesites];
  // for(int i=0;i<8*volnodesites;i++) orig[i] = gauge[i];

  // lattice.Convert(WILSON);
  // lattice.Convert(CANONICAL);

  // gauge = lattice.GaugeField();
  
  // for(int i=0;i<8*volnodesites;i++){
  //   bool linepass =true;
  //   for(int j=0;j<18;j++){
  //     if(gauge[i].elem(j) != orig[i].elem(j)){
  // 	printf("Err %d %d; %e %e\n",i,j,gauge[i].elem(j),orig[i].elem(j));
  // 	linepass = false;
  //     }
  //   }
  //   if(linepass) printf("Matrix %d pass\n",i);
  // }
  

  // exit(0);
  // printf("Post convert\n");

  // printf("cb0:\n");
  // int cbsize = 2*volnodesites;

  // for(int i=0;i<cbsize;i++){
  //   printf("U  i=%d ",i);
  //   for(int j=0;j<18;j++){
  //     if(gauge[i].elem(j)>=0)
  // 	printf("%d: %.3f ",j,gauge[i].elem(j));
  //     else
  // 	printf("%d:%.3f ",j,gauge[i].elem(j));
  //   }
  //   printf("\n");

  //   if(GJP.Gparity()){
  //     printf("U* i=%d ",i+cbsize);
  //     for(int j=0;j<18;j++){
  // 	if(gauge[i+cbsize].elem(j)>=0)
  // 	  printf("%d: %.3f ",j,gauge[i+cbsize].elem(j));
  // 	else
  // 	  printf("%d:%.3f ",j,gauge[i+cbsize].elem(j));
  //     }
  //     printf("\n\n");

  //   }


  // }
  // printf("cb1:\n");
  // int offset = cbsize;
  // if(GJP.Gparity()) offset*=2;

  // for(int i=offset;i<offset+cbsize;i++){
  //   printf("U  i=%d ",i);
  //   for(int j=0;j<18;j++){
  //     if(gauge[i].elem(j)>=0)
  // 	printf("%d: %.3f ",j,gauge[i].elem(j));
  //     else
  // 	printf("%d:%.3f ",j,gauge[i].elem(j));
  //   }
  //   printf("\n");

  //   if(GJP.Gparity()){
  //     printf("U* i=%d ",i+cbsize);
  //     for(int j=0;j<18;j++){
  // 	if(gauge[i+cbsize].elem(j)>=0)
  // 	  printf("%d: %.3f ",j,gauge[i+cbsize].elem(j));
  // 	else
  // 	  printf("%d:%.3f ",j,gauge[i+cbsize].elem(j));
  //     }
  //     printf("\n\n");

  //   }


  // }


  // printf("Gauge field:\n\n");

  // int base_loc[4];
  // for(int i=0;i<4;i++) base_loc[i] = GJP.NodeCoor(i)*GJP.NodeSites(i); //absolute pos of start point
    
  // //int pos[4];

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
  // 	    printf("%d %d %d %d; %d:\n",x+base_loc[0],y+base_loc[1],z+base_loc[2],t+base_loc[3],mu);
  // 	    for(int i=0;i<3;i++){
  // 	    	for(int j=0;j<3;j++){
  // 	    	  const Complex &c = (*U)(i,j);
  // 	    	  printf("[%.3f %.3f] ",c.real(), c.imag());
  // 	    	}
  // 	    	printf("\n");
  // 	    }
  // 	    printf("\n");
  // 	  }
	    

  // 	}
  //     }
  //   }
  // }

  // if(GJP.Gparity()){
  //   printf("Conjugated Gauge field:\n\n");

  //   for(int t=0;t<GJP.NodeSites(3);t++){
  //     pos[3] = t;
  //     for(int z=0;z<GJP.NodeSites(2);z++){
  // 	pos[2] = z;
  // 	for(int y=0;y<GJP.NodeSites(1);y++){
  // 	  pos[1] = y;
  // 	  for(int x=0;x<GJP.NodeSites(0);x++){
  // 	    pos[0] = x;
	    
  // 	    for(int mu=0;mu<4;mu++){
  // 	      Matrix *U = const_cast<Matrix *>(lattice.GetLink(pos, mu, 1));
  // 	      printf("%d %d %d %d; %d:\n",x+base_loc[0],y+base_loc[1],z+base_loc[2],t+base_loc[3],mu);
  // 	      for(int i=0;i<3;i++){
  // 	    	for(int j=0;j<3;j++){
  // 	    	  const Complex &c = (*U)(i,j);
  // 	    	  printf("[%.3f %.3f] ",c.real(), c.imag());
  // 	    	}
  // 	    	printf("\n");
  // 	      }
  // 	      printf("\n");
  // 	    }
	    

  // 	  }
  // 	}
  //     }
  //   }

  // }



  // // PropagatorContainer prop;
  // //  prop.setup(prop_args.props.props_val[0]);

   

  //  QPropW *src = &  PropManager::getProp(prop_args.props.props_val[0].generics.tag).getProp(lattice);  

  //  //QPropW *src = &prop.getProp(lattice);


  // // exit(0);

  // // CgArg cg;
  // // cg.mass =  0.5;
  // // cg.max_num_iter = 5000;
  // // cg.stop_rsd =   1.0000000000000000e-06;
  // // cg.true_rsd =   1.0000000000000000e-06;
  // // cg.RitzMatOper = NONE;
  // // cg.Inverter = CG;
  // // cg.bicgstab_n = 0;


  // // QPropWArg qpropw_arg;
  // // qpropw_arg.cg = cg;
  // // qpropw_arg.x = 0;
  // // qpropw_arg.y = 0;
  // // qpropw_arg.z = 0;
  // // qpropw_arg.t = 0;
  // // if(gparity){
  // //   qpropw_arg.flavor = 0; //point source on d field
  // // }
  // // qpropw_arg.ensemble_label = "ens";
  // // qpropw_arg.ensemble_id = "ens_id";
  // // qpropw_arg.StartSrcSpin = 0;
  // // qpropw_arg.EndSrcSpin = 4;
  // // qpropw_arg.StartSrcColor = 0;
  // // qpropw_arg.EndSrcColor = 3;

  // // qpropw_arg.save_prop = 1;
  // // qpropw_arg.file = "prop.dat";
  
  // // CommonArg common_arg("label","filename");
  
  // // QPropW* src = new QPropWPointSrc(lattice, &qpropw_arg, &common_arg);
  

  
  

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







// class PropagatorContainer{
// protected:
//   QPropW *prop;
//   AttributeContainer* attributes[50]; //positions are mapped by the integer values of the enum  AttrType with no duplicate entries
//   AttributeContainer* findAttr(const AttrType &type) const; //returns NULL if not present
// public:
//   void add(const AttributeContainer &p);
//   void setup(PropagatorArg &arg);
//   template<typename T>
//   bool getAttr(T* &to) const;
  
//   void readProp(); //loads prop if on disk, does nothing otherwise
//   void calcProp(Lattice &latt);
//   void deleteProp(); //deletes the prop from memory, used to save space. Will be automatically recalculated if getProp is called again

//   QPropW & getProp(Lattice &latt); //get the prop, calculate or load if necessary

//   bool tagEquals(const char* what); //check if tag is equal to input

//   int flavor() const;

//   PropagatorContainer();
//   ~PropagatorContainer();
// };

// PropagatorContainer::PropagatorContainer(): prop(NULL){ for(int i=0;i<50;i++) attributes[i]=NULL;}


// void PropagatorContainer::setup(PropagatorArg &arg){
//   for(int i=0;i<50;i++) if(attributes[i]) delete attributes[i];
  
//   //special case, generic attributes are always present and maintained outside of the array
//   AttributeContainer* p = new AttributeContainer;
//   p->type = GENERIC_PROP_ATTR;
//   p->AttributeContainer_u.generic_prop_attr = arg.generics;
//   attributes[ (int)p->type ] = p;

//   printf("Got %d additional attributes\n",arg.attributes.attributes_len);

//   for(int i=0;i<arg.attributes.attributes_len;i++){
//     //no duplicates, later entries of same type overwrite earlier
//     p = &arg.attributes.attributes_val[i];
//     int idx = (int)p->type;
//     if(attributes[idx]) delete attributes[idx];
//     attributes[idx] = new AttributeContainer(*p); //make a copy
//   }
// }

// void PropagatorContainer::add(const AttributeContainer &p){
//   int idx = (int)p.type;
//   if(attributes[idx]!=NULL) delete attributes[idx]; //no duplicates
//   attributes[idx] = new AttributeContainer(p);
// }

// AttributeContainer*  PropagatorContainer::findAttr(const AttrType &type) const{
//   return attributes[(int)type];
// }

// template<typename T>
// bool PropagatorContainer::getAttr(T* &to) const{
//   AttributeContainer *c = findAttr(T::getType());
//   if(c==NULL) return false;
//   to = reinterpret_cast<T*>(&(c->AttributeContainer_u)); //C99 standard: A pointer to a union object, suitably converted, points to each of its members
//   return true;
// }

// PropagatorContainer::~PropagatorContainer(){
//   for(int i=0;i<50;i++) if(attributes[i]) delete attributes[i];
//   if(prop!=NULL) delete prop;
// }


// void PropagatorContainer::readProp(){
//   if(prop) return; //don't load if already inverted
//   PropIOAttrArg *io;
//   if(!getAttr(io)) return;
//   if(!io->prop_on_disk) return;

//   //load prop with info from io
//   prop->ReLoad(io->qio_filename);
// }

// void PropagatorContainer::calcProp(Lattice &latt){
//   //function acts as factory for QPropW objects depending on the attribute objects
//   if(prop) return; //don't calculate twice

//   char *cname = "PropagatorContainer";
//   char *fname = "calcProp()";

//   CommonArg c_arg("label","filename");//find out what this does!

//   GenericPropAttrArg *generics;
//   if(!getAttr(generics)) ERR.General(cname,fname,"Propagator attribute list does not contain a GenericPropAttr\n");

//   printf("Calculating propagator %s\n",generics->tag);

//   CgArg cg;
//   cg.mass =  generics->mass;
//   cg.max_num_iter = 5000;
//   cg.stop_rsd =   1.0000000000000000e-06;
//   cg.true_rsd =   1.0000000000000000e-06;
//   cg.RitzMatOper = NONE;
//   cg.Inverter = CG;
//   cg.bicgstab_n = 0;

//   QPropWArg qpropw_arg;
//   qpropw_arg.cg = cg;
//   qpropw_arg.x = 0;
//   qpropw_arg.y = 0;
//   qpropw_arg.z = 0;
//   qpropw_arg.t = 0;
//   qpropw_arg.flavor = 0; //default on d field  
//   qpropw_arg.ensemble_label = "ens";
//   qpropw_arg.ensemble_id = "ens_id";
//   qpropw_arg.StartSrcSpin = 0;
//   qpropw_arg.EndSrcSpin = 4;
//   qpropw_arg.StartSrcColor = 0;
//   qpropw_arg.EndSrcColor = 3;

//   //fill out qpropw_arg arguments
//   GparityFlavorAttrArg *flav;
//   if(getAttr(flav)) qpropw_arg.flavor = flav->flavor;

//   PropIOAttrArg *io;
//   if(getAttr(io)){ 
//     qpropw_arg.save_prop = io->save_to_disk;
//     qpropw_arg.file = io->qio_filename;
//   }

//   PointSourceAttrArg *pt;
//   WallSourceAttrArg *wl;
//   if(getAttr(pt) && getAttr(wl) )
//     ERR.General(cname,fname,"Propagator %s cannot have both WallSourceAttrArg and PointSourceAttrArg\n",generics->tag);

//   if(getAttr(pt)){
//     qpropw_arg.x = pt->pos[0];
//     qpropw_arg.y = pt->pos[1];
//     qpropw_arg.z = pt->pos[2];
//     qpropw_arg.t = pt->pos[3];
    
//     //assemble the arg objects
//     prop = new QPropWPointSrc(latt,&qpropw_arg,&c_arg);

//     MomentumAttrArg *mom;
//     if(getAttr(mom)){
//       //for a point source with momentum, apply the e^-ipx factor at the source location 
//       //such that it can be treated in the same was as a momentum source
//       const Float Pi_const=3.141592654;
//       Float pdotx=0.0;
//       pdotx +=((Float) mom->p[0]*pt->pos[0])/((Float) GJP.XnodeSites()*GJP.Xnodes() )
// 	    + ((Float) mom->p[1]*pt->pos[1])/((Float) GJP.YnodeSites()*GJP.Ynodes() )
//       	    + ((Float) mom->p[2]*pt->pos[2])/((Float) GJP.ZnodeSites()*GJP.Znodes() );
//       Rcomplex mom_fac(cos(pdotx),sin(pdotx));
//       for(int i=0;i<GJP.VolNodeSites();i++) (*prop)[i] *= mom_fac; 
//     }
//   }else if(getAttr(wl)){
//     qpropw_arg.t = wl->t;
    
//     MomentumAttrArg *mom;
//     if(getAttr(mom)){
//       prop = new QPropWMomSrc(latt,&qpropw_arg,mom->p,&c_arg);
//     }else{
//       prop = new QPropWWallSrc(latt,&qpropw_arg,&c_arg);
//     }
//   }else{
//     ERR.General(cname,fname,"Propagator %s has no source type AttrArg\n",generics->tag);
//   }
  
//   if(getAttr(io)){
//     io->prop_on_disk = true; //QPropW saves the prop
//   }
// }


// int PropagatorContainer::flavor() const{
//   GparityFlavorAttrArg *flav;
//   if(getAttr(flav)) return flav->flavor;
//   return 0;
// }

// QPropW & PropagatorContainer::getProp(Lattice &latt){
//   if(!prop){
//     readProp();
//     calcProp(latt); //will calculate if prop was not read
//   }
//   return *prop;
// }

// bool PropagatorContainer::tagEquals(const char* what){
//   GenericPropAttrArg *generics;
//   if(!getAttr(generics)) ERR.General("PropagatorContainer","tagEquals(const char* what)","Propagator attribute list does not contain a GenericPropAttr\n");
//   if(strcmp(generics->tag,what)==0) return true;
//   return false;
// }

// void PropagatorContainer::deleteProp(){
//   if(prop) delete prop;
// }


// //container for multiple props with destructor that deletes all props
// class PropVector{
// public:
//   const static int MAX_SIZE = 100;
//   PropagatorContainer & operator[](const int &idx);
//   PropagatorContainer & addProp(PropagatorArg &arg);
//   const int &size() const;
//   PropVector();
//   ~PropVector();
// private:
//   PropagatorContainer *props[MAX_SIZE];
//   int sz;
// };
// PropagatorContainer & PropVector::operator[](const int &idx){ return *props[idx]; }
// PropVector::PropVector(): sz(0){ for(int i=0;i<MAX_SIZE;i++) props[i] = NULL; }
// PropVector::~PropVector(){ for(int i=0;i<MAX_SIZE;i++) if(props[i]!=NULL) delete props[i]; }
  
// const int &PropVector::size() const{ return sz; }

// PropagatorContainer & PropVector::addProp(PropagatorArg &arg){
//   if(sz==MAX_SIZE){ ERR.General("PropVector","addProp(PropagatorArg &arg)","Reached maximum number of allowed propagators: %d\n",MAX_SIZE); }
//   PropagatorContainer* p = new PropagatorContainer;
//   p->setup(arg);
//   props[sz++] = p;
//   return *p;
// }

// //static class for managing propagators
// class PropManager{
// private:
//   static PropVector props;
// public:
//   static PropagatorContainer & getProp(const char *tag);
//   static PropagatorContainer & addProp(PropagatorArg &arg); //this does not actually invert the propagator until getProp is called on the PropagatorContainer or using the above function on the vector
//   static void setup(JobPropagatorArgs &prop_args);
// };
// PropVector PropManager::props;

// PropagatorContainer & PropManager::getProp(const char *tag){
//   for(int i=0;i<props.size();i++) if(props[i].tagEquals(tag)) return props[i];
//   ERR.General("PropManager","getProp(const char *tag, Lattice &latt)","Prop '%s' does not exist!\n",tag);
// };

// PropagatorContainer & PropManager::addProp(PropagatorArg &arg){
//   return props.addProp(arg);
// }

// void PropManager::setup(JobPropagatorArgs &prop_args){
//   for(int i=0;i< prop_args.props.props_len; i++){
//     addProp(prop_args.props.props_val[i]);
//   }
// }


// class CorrelationFunction{
//   char *label;
//   int ncontract;
//   int time_size;
//   Rcomplex** wick; //the wick contractions, each a function of t

//   void sumLattice(); //sum each element of wick[contraction][t] over all nodes
// public:
//   void setNcontractions(const int &n);
//   void write(const char *file);
//   void write(FILE *fp);
//   Rcomplex & operator()(const int &contraction_idx, const int &t);
//   CorrelationFunction(const char *_label);
//   ~CorrelationFunction();
// };
// CorrelationFunction::CorrelationFunction(const char *_label): wick(NULL),ncontract(0){
//   time_size=GJP.Tnodes()*GJP.TnodeSites();
//   label = new char[strlen(_label)];
//   strcpy(label,_label);
// }
// void CorrelationFunction::setNcontractions(const int &n){
//   if(ncontract!=0){
//     sfree(wick[0]);
//     sfree(wick);
//   }
//   ncontract = n;
//   wick = (Rcomplex**)smalloc(n*sizeof(Rcomplex*)); //allocate array of pointers (I wish we could just use std::vector here!)

//   //contiguous memory slot
//   Rcomplex *stack = (Rcomplex*)smalloc(n*time_size*sizeof(Rcomplex));

//   for(int i=0;i<n;i++){
//     wick[i] = stack+i*time_size;
//     for(int t=0;t<time_size;t++){ wick[i][t].real(0.0); wick[i][t].imag(0.0); }
//   }
// }
// CorrelationFunction::~CorrelationFunction(){
//   if(ncontract!=0){
//     sfree(wick[0]); //free the contiguous memory block holding the data
//     sfree(wick); //delete the pointer array holding this
//   }
// }
// Rcomplex & CorrelationFunction::operator()(const int &contraction_idx, const int &t){
//   return wick[contraction_idx][t];
// }
// void CorrelationFunction::sumLattice(){
//   for(int i=0;i<ncontract;i++){
//     for(int t=0;t<time_size;t++){
//       slice_sum( (Float*)&wick[i][t], 2, 99); //2 for re/im, 99 is a *magic* number (we are abusing slice_sum here)
//     }
//   }
// }
// void CorrelationFunction::write(const char *file){
//   FILE *fp;
//   if ((fp = Fopen(file, "w")) == NULL) {
//     ERR.FileW("CorrelationFunction","write(const char *file)",file);
//   }
//   write(fp);
//   Fclose(fp);
// }
// void CorrelationFunction::write(FILE *fp){
//   sumLattice(); //sum the correlation function over all nodes

//   Fprintf(fp,"%s\n",label);
//   Fprintf(fp,"%d contractions\n",ncontract);

//   for(int c=0;c<ncontract;c++){
//     Fprintf(fp,"Contraction %d\n",c);

//     for(int t=0; t<time_size; t++)
//       Fprintf(fp, "%d %.16e %.16e\n",t, wick[c][t].real(), wick[c][t].imag());
//   }
// }


  
// class AlgGparityContract : public Alg{
// private:
//   char *cname;
//   GparityContractArg *args;
// public:
//   AlgGparityContract(Lattice & latt, CommonArg& c_arg, GparityContractArg& arg);
//   void run();
//   void spectrum(const GparityMeasurement &measargs);

//   void pi_plus(const char *q_f0_tag, const char *q_f1_tag, CorrelationFunction &corr);
//   void pi_minus(const char *q_f0_tag, const char *q_f1_tag, CorrelationFunction &corr);
// };
// AlgGparityContract::AlgGparityContract(Lattice & latt, CommonArg& c_arg, GparityContractArg& arg): Alg(latt,&c_arg), args(&arg){ cname = "AlgGparityContract"; }

// void AlgGparityContract::run(){
//   for(int i=0;i<args->meas.meas_len;i++){
//     spectrum(args->meas.meas_val[i]);   
//   }
// }

// void AlgGparityContract::spectrum(const GparityMeasurement &measargs){
//   char * prop_f[3]; bool got_f[3] = {false,false,false};

//   PropagatorContainer &p1 = PropManager::getProp(measargs.prop_1);
//   PropagatorContainer &p2 = PropManager::getProp(measargs.prop_2);
  
//   int f1 = p1.flavor(); int f2 = p2.flavor();
//   prop_f[f1] = measargs.prop_1; got_f[f1] = true;
//   prop_f[f2] = measargs.prop_2; got_f[f2] = true;

//   char out[100];

//   if(got_f[0] && got_f[1]){
//     sprintf(out,"%s pi^+",measargs.label_stub);
//     CorrelationFunction c_piplus(out);
//     sprintf(out,"%s pi^-",measargs.label_stub);
//     CorrelationFunction c_piminus(out);

//     if(UniqueID()==0) printf("Doing pi^+ contractions\n");
//     pi_plus(prop_f[0],prop_f[1],c_piplus);
//     if(UniqueID()==0) printf("Doing pi^- contractions\n");
//     pi_minus(prop_f[0],prop_f[1],c_piminus);
    
//     sprintf(out,"%s_pi_plus.dat",measargs.file_stub);
//     c_piplus.write(out);
//     sprintf(out,"%s_pi_minus.dat",measargs.file_stub);
//     c_piminus.write(out);
//   }else{
//     ERR.General(cname, "spectrum(const char *prop_1, const char *prop_2)", "Could not find any valid contractions for propagators %s and %s\n",measargs.prop_1,measargs.prop_2);
//   }

// }

// void AlgGparityContract::pi_plus(const char *q_f0_tag, const char *q_f1_tag, CorrelationFunction &corr){
//   int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
//   int spatial_vol = GJP.VolNodeSites()/GJP.TnodeSites();

//   corr.setNcontractions(2);

//   WilsonMatrix trace_a, trace_b, tmp;
//   PropagatorContainer &q_f0_pc = PropManager::getProp(q_f0_tag);
//   PropagatorContainer &q_f1_pc = PropManager::getProp(q_f1_tag);

//   QPropW &q_f0_qpw = q_f0_pc.getProp(AlgLattice());
//   QPropW &q_f1_qpw = q_f1_pc.getProp(AlgLattice());
  
//   Rcomplex tmpc;

//   for(int i=0;i<GJP.VolNodeSites();i++){
//     int t = i/spatial_vol + shift_t;

//     // //DEBUG
//     // int shift_x = GJP.XnodeCoor()*GJP.XnodeSites();
//     // int shift_y = GJP.YnodeCoor()*GJP.YnodeSites();
//     // int shift_z = GJP.ZnodeCoor()*GJP.ZnodeSites();
//     // int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
//     // //Local lattice dimensions:
//     // int size_x = GJP.XnodeSites();
//     // int size_y = GJP.YnodeSites();
//     // int size_z = GJP.ZnodeSites();
//     // int size_t = GJP.TnodeSites();
//     // int size_xy = size_x*size_y;
//     // int vol = (GJP.VolNodeSites()/GJP.TnodeSites()); // =size_x*size_y_size_z
    
//     // int pos[4];
//     // pos[3] = i/vol + shift_t;
//     // pos[2] = (i%vol)/size_xy + shift_z;
//     // pos[1] = (i%size_xy)/size_x + shift_y;
//     // pos[0] = i%size_x + shift_x;

//     // printf("Piplus (%d,%d,%d,%d)\n",pos[0],pos[1],pos[2],pos[3]);
//     // //DEBUG


//     //Doing contraction 1 of 2
//     Rcomplex &wick0 = corr(0,t);
//     tmpc.real(1.0); tmpc.imag(0.0);

//     //Doing trace 1 of 1
//     trace_a = q_f0_qpw.SiteMatrix(i,0);
//     trace_a.hconj();
//     trace_a.gl(3).gl(1).gl(-5);
//     trace_a.gr(-5).gr(1).gr(3);
//     trace_b = q_f1_qpw.SiteMatrix(i,1);
//     trace_b.cconj();
//     tmpc*= Trace(trace_a,trace_b);

//     wick0 += tmpc;
//     //printf("Wick0 -> %e,%e\n",wick0.real(),wick0.imag());

//     //Doing contraction 2 of 2
//     Rcomplex &wick1 = corr(1,t);
//     tmpc.real(1.0); tmpc.imag(0.0);

//     //Doing trace 1 of 1
//     trace_a = q_f1_qpw.SiteMatrix(i,0);
//     trace_a.hconj();
//     trace_a.gl(3).gl(1).gl(-5);
//     trace_a.gr(-5).gr(1).gr(3);
//     trace_b = q_f0_qpw.SiteMatrix(i,1);
//     trace_b.cconj();
//     tmpc*= Trace(trace_a,trace_b);

//     wick1 += tmpc;
//     //printf("Wick1 -> %e,%e\n",wick1.real(),wick1.imag());
//   }
// }

// void AlgGparityContract::pi_minus(const char *q_f0_tag, const char *q_f1_tag, CorrelationFunction &corr){
//   int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
//   int spatial_vol = GJP.VolNodeSites()/GJP.TnodeSites();

//   corr.setNcontractions(2);

//   WilsonMatrix trace_a, trace_b, tmp;
//   PropagatorContainer &q_f0_pc = PropManager::getProp(q_f0_tag);
//   PropagatorContainer &q_f1_pc = PropManager::getProp(q_f1_tag);

//   QPropW &q_f0_qpw = q_f0_pc.getProp(AlgLattice());
//   QPropW &q_f1_qpw = q_f1_pc.getProp(AlgLattice());

//   Rcomplex tmpc;
  
//   for(int i=0;i<GJP.VolNodeSites();i++){
//     int t = i/spatial_vol + shift_t;

//     //Doing contraction 1 of 2
//     Rcomplex &wick0 = corr(0,t); 
//     tmpc.real(1.0); tmpc.imag(0.0);

//     //Doing trace 1 of 1
//     trace_a = q_f0_qpw.SiteMatrix(i,0);
//     trace_a.gl(3).gl(1).gl(-5);
//     trace_a.gr(-5).gr(1).gr(3);
//     trace_b = q_f1_qpw.SiteMatrix(i,1);
//     trace_b.transpose();
//     tmpc*= Trace(trace_a,trace_b);

//     wick0 += tmpc;

//     //Doing contraction 2 of 2
//     Rcomplex &wick1 = corr(1,t); 
//     tmpc.real(1.0); tmpc.imag(0.0);

//     //Doing trace 1 of 1
//     trace_a = q_f0_qpw.SiteMatrix(i,1);
//     trace_a.gr(-5).gr(1).gr(3);
//     trace_b = q_f1_qpw.SiteMatrix(i,0);
//     trace_b.transpose();
//     trace_b.gr(-5).gr(1).gr(3);
//     tmpc*= Trace(trace_a,trace_b);

//     wick1 += tmpc;
//   }
// }



// struct AttributeContainer2 {
//   AttrType type;
//   void *ptr;
//   // AttributeContainer2(): ptr(NULL){ printf("AttributeContainer2 constructor\n"); };
//   // void clear();
//   // ~AttributeContainer2();
// };
// typedef struct AttributeContainer2 AttributeContainer2;

// //class VML;
// class PropagatorArg2 {
// public:
// 	 bool Encode(char *filename,char *instance);
// 	 bool Decode(char *filename,char *instance);
// 	 bool Vml(VML *vmls,char *instance);
// 	GenericPropAttrArg generics;
// 	struct {
// 		u_int attributes_len;
// 		AttributeContainer2 *attributes_val;
// 	} attributes;
// 	   ~PropagatorArg2 (  ) ;
// };

// //class VML;
// class JobPropagatorArgs2 {
// public:
// 	 bool Encode(char *filename,char *instance);
// 	 bool Decode(char *filename,char *instance);
// 	 bool Vml(VML *vmls,char *instance);
// 	struct {
// 		u_int props_len;
// 		PropagatorArg2 *props_val;
// 	} props;
// 	   ~JobPropagatorArgs2 (  ) ;
// };

// #ifdef __cplusplus
// extern "C" {
// #endif

// #if defined(__STDC__) || defined(__cplusplus)
// extern  bool_t vml_AttributeContainer2 (VML *, char *instance, AttributeContainer2*);
// extern  bool_t vml_PropagatorArg2 (VML *, char *instance, PropagatorArg2*);
// extern  bool_t vml_JobPropagatorArgs2 (VML *, char *instance, JobPropagatorArgs2*);

// #else /* K&R C */
// extern  bool_t vml_AttributeContainer2 (VML *, char *instance, AttributeContainer2*);
// extern  bool_t vml_PropagatorArg2 (VML *, char *instance, PropagatorArg2*);
// extern  bool_t vml_JobPropagatorArgs2 (VML *, char *instance, JobPropagatorArgs2*);

// #endif /* K&R C */

// #ifdef __cplusplus
// }
// #endif


// // AttributeContainer2::~AttributeContainer2(){
// //   clear();
// // }
// // void AttributeContainer2::clear(){
// //   if(ptr == NULL) return;
  
// //   switch (type) {
// //   case GENERIC_PROP_ATTR:
// //     delete (GenericPropAttrArg*)ptr; ptr = NULL; break;
// //   case POINT_SOURCE_ATTR:
// //     delete (PointSourceAttrArg*)ptr; ptr = NULL; break;
// //   case WALL_SOURCE_ATTR:
// //     delete (PointSourceAttrArg*)ptr; ptr = NULL; break;
// //   case MOMENTUM_ATTR:
// //     delete (MomentumAttrArg*)ptr; ptr = NULL; break;
// //   case PROP_IO_ATTR:
// //     delete (PropIOAttrArg*)ptr; ptr = NULL; break;
// //   case GPARITY_FLAVOR_ATTR:
// //     delete (GparityFlavorAttrArg*)ptr; ptr = NULL; break;
// //   case CG_ATTR:
// //     delete (CGAttrArg*)ptr; ptr = NULL; break;
// //   case GAUGE_FIX_ATTR:
// //     delete (GaugeFixAttrArg*)ptr; ptr = NULL; break;
// //   }
// // }

// bool_t
// vml_AttributeContainer2 (VML *vmls, char *name,AttributeContainer2 *objp)
// {
//   if (!vml_AttrType (vmls, "type", &objp->type)){
//     printf("VML attrtype read fail\n");
//     return FALSE;
//   }

//   if(vmls->x_op == VML_DECODE) printf("Decoding AttributeContainer2, got type %s\n",AttrType_map[objp->type].name);
  
//   printf("ptr has initial value %p\n",objp->ptr);

//   if(objp->ptr == NULL && (vmls->x_op == VML_ENCODE || vmls->x_op == VML_FREE) ) return TRUE;

//   switch (objp->type) {
//   case GENERIC_PROP_ATTR:
//     if(vmls->x_op == VML_DECODE) objp->ptr = (void*)(new GenericPropAttrArg);

//     if (!vml_GenericPropAttrArg (vmls, "generic_prop_attr", (GenericPropAttrArg*)objp->ptr))
//       return FALSE;

//     if(vmls->x_op == VML_FREE){ delete (GenericPropAttrArg*)objp->ptr; objp->ptr = NULL; }

//     break;
//   case POINT_SOURCE_ATTR:
//     if(objp->ptr == NULL && vmls->x_op == VML_DECODE) objp->ptr = (void*)(new PointSourceAttrArg);

//     if (!vml_PointSourceAttrArg (vmls, "point_source_attr", (PointSourceAttrArg*)objp->ptr))
//       return FALSE;

//     if(vmls->x_op == VML_FREE){ delete (PointSourceAttrArg*)objp->ptr; objp->ptr = NULL; }

//     break;
//   case WALL_SOURCE_ATTR:
//     if(objp->ptr == NULL && vmls->x_op == VML_DECODE) objp->ptr = (void*)(new WallSourceAttrArg);

//     if (!vml_WallSourceAttrArg (vmls, "wall_source_attr", (WallSourceAttrArg*)objp->ptr))
//       return FALSE;

//     if(vmls->x_op == VML_FREE){ delete (WallSourceAttrArg*)objp->ptr; objp->ptr = NULL; }

//     break;
//   case MOMENTUM_ATTR:
//     if(objp->ptr == NULL && vmls->x_op == VML_DECODE) objp->ptr = (void*)(new MomentumAttrArg);

//     if (!vml_MomentumAttrArg (vmls, "momentum_attr", (MomentumAttrArg*)objp->ptr))
//       return FALSE;

//     if(vmls->x_op == VML_FREE){ delete (MomentumAttrArg*)objp->ptr; objp->ptr = NULL; }

//     break;
//   case PROP_IO_ATTR:
//     if(objp->ptr == NULL && vmls->x_op == VML_DECODE) objp->ptr = (void*)(new PropIOAttrArg);

//     if (!vml_PropIOAttrArg (vmls, "prop_io_attr", (PropIOAttrArg*)objp->ptr))
//       return FALSE;

//     if(vmls->x_op == VML_FREE){ delete (PropIOAttrArg*)objp->ptr; objp->ptr = NULL; }

//     break;
//   case GPARITY_FLAVOR_ATTR:
//     if(objp->ptr == NULL && vmls->x_op == VML_DECODE) objp->ptr = (void*)(new GparityFlavorAttrArg);
    
//     if (!vml_GparityFlavorAttrArg (vmls, "gparity_flavor_attr", (GparityFlavorAttrArg*)objp->ptr))
//       return FALSE;

//     if(vmls->x_op == VML_FREE){ delete (GparityFlavorAttrArg*)objp->ptr; objp->ptr = NULL; }

//     break;
//   case CG_ATTR:
//     if(objp->ptr == NULL && vmls->x_op == VML_DECODE) objp->ptr = (void*)(new CGAttrArg);

//     if (!vml_CGAttrArg (vmls, "cg_attr", (CGAttrArg*)objp->ptr))
//       return FALSE;

//     if(vmls->x_op == VML_FREE){ delete (CGAttrArg*)objp->ptr; objp->ptr = NULL; }

//     break;
//   case GAUGE_FIX_ATTR:
//     if(objp->ptr == NULL && vmls->x_op == VML_DECODE) objp->ptr = (void*)(new GaugeFixAttrArg);

//     if (!vml_GaugeFixAttrArg (vmls, "gauge_fix_attr", (GaugeFixAttrArg*)objp->ptr))
//       return FALSE;

//     if(vmls->x_op == VML_FREE){ delete (GaugeFixAttrArg*)objp->ptr; objp->ptr = NULL; }

//     break;
//   default:
//     return FALSE;
//   }
//   return TRUE;
// }


// bool PropagatorArg2::Encode(char *filename,char *instance){
//   VML vmls;
//   if ( !vmls.Create(filename,VML_ENCODE)) return false;
//   if ( !Vml(&vmls,instance) ) return false;
//   vmls.Destroy(); return true;
// }

// bool PropagatorArg2::Decode(char *filename,char *instance){
//   VML vmls;
//   if ( !vmls.Create(filename,VML_DECODE)) return false;
//   if ( !Vml(&vmls,instance)) return false;
//   vmls.Destroy(); return true;
// }
// bool PropagatorArg2::Vml(VML *vmls,char *instance){
//   if(!vml_PropagatorArg2(vmls,instance,this)) return false;
//   return true;
// }


// bool_t
// vml_PropagatorArg2 (VML *vmls, char *name,PropagatorArg2 *objp)
// {
//   vml_class_begin(vmls,"PropagatorArg2",name);
//   if (!vml_GenericPropAttrArg (vmls, "generics", &objp->generics))
//     return FALSE;
//   if (!vml_array (vmls, "attributes", (char **)&objp->attributes.attributes_val, (u_int *) &objp->attributes.attributes_len, ~0,
// 		  sizeof (AttributeContainer2), (vmlproc_t) vml_AttributeContainer2))
//     return FALSE;
//   vml_class_end(vmls,"PropagatorArg2",name);
//   return TRUE;
// }
// bool JobPropagatorArgs2::Encode(char *filename,char *instance){
//   VML vmls;
//   if ( !vmls.Create(filename,VML_ENCODE)) return false;
//   if ( !Vml(&vmls,instance) ) return false;
//   vmls.Destroy(); return true;
// }

// bool JobPropagatorArgs2::Decode(char *filename,char *instance){
//   VML vmls;
//   if ( !vmls.Create(filename,VML_DECODE)) return false;
//   if ( !Vml(&vmls,instance)) return false;
//   vmls.Destroy(); return true;
// }
// bool JobPropagatorArgs2::Vml(VML *vmls,char *instance){
//   if(!vml_JobPropagatorArgs2(vmls,instance,this)) return false;
//   return true;
// }


// bool_t
// vml_JobPropagatorArgs2 (VML *vmls, char *name,JobPropagatorArgs2 *objp)
// {
//   vml_class_begin(vmls,"JobPropagatorArgs2",name);
//   if (!vml_array (vmls, "props", (char **)&objp->props.props_val, (u_int *) &objp->props.props_len, ~0,
// 		  sizeof (PropagatorArg2), (vmlproc_t) vml_PropagatorArg2))
//     return FALSE;
//   vml_class_end(vmls,"JobPropagatorArgs2",name);
//   return TRUE;
// }


// PropagatorArg2::~PropagatorArg2(){
//   if(attributes.attributes_val) delete attributes.attributes_val;
// }
  
// JobPropagatorArgs2::~JobPropagatorArgs2(){
//   if(props.props_val) delete props.props_val;
// }


