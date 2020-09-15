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

#include<sstream>
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

EvoArg evo_arg;
DoArg do_arg;
NoArg no_arg;

void checkpoint(int traj);
void checkReproduces(const char *file, const bool &disable_ustar_reconstruction, const Float &tolerance = 1e-10);

int main(int argc, char *argv[])
{ 
  char plaq_file[256];
  char pbp_file[256];
  char wline_file[256];
  char hmc_file[256];

  char *cname=argv[0];
  char *fname="main()";
  Float dtime;
  
  CommonArg common_arg_hmc;
  CommonArg common_arg_plaq;
  CommonArg common_arg_pbp;
  CommonArg common_arg_wline;

  Start(&argc,&argv);

  bool dbl_latt_storemode(false);
  bool load_config(false);
  char *load_config_file;
  bool save_lrg(false);
  char *save_lrg_file;
  bool load_lrg(false);
  char *load_lrg_file;
  
  if ( argc<6 ) { 
    printf("Args:\t <options> do_arg.vml hmc_arg.vml evo_arg.vml\n");
    printf("\tgauge_arg.vml ab1_arg.vml\n");
    exit(-1);
  }

  //Option to disable U* reconstruction for use when loading old lattices that saved both U and U* (does not affect lattices being saved)
  bool disable_ustar_reconstruction(false);

  //When config 'test_conf_reproduce' has been run, load 'test_conf_reproduce_file' and ensure they are equal
  int test_conf_reproduce = -1;
  char* test_conf_reproduce_file;

  int i=1;
  while(1){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-load_config",15) == 0){
      if(i==argc-1){ printf("-save_config requires an argument\n"); exit(-1); }
      load_config=true;
      load_config_file = argv[i+1];
      i+=2;
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
    }else if( strncmp(cmd,"-test_conf_reproduce",30) == 0){
      if(i==argc-2){ printf("-test_conf_reproduce 2 arguments\n"); exit(-1); }
      std::stringstream ss; ss<< argv[i+1]; ss>>test_conf_reproduce;
      test_conf_reproduce_file = argv[i+2];
      i+=3;
      printf("Testing config %d reproduces lattice %s\n",test_conf_reproduce,test_conf_reproduce_file);
    }else if( strncmp(cmd,"-disable_ustar_reconstruction",30) == 0){
      disable_ustar_reconstruction = true;
      printf("Disabling U* reconstruction on gauge field load\n");
      i++;
    }else if( cmd[0] != '-' ){
      break;
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }
  SerialIO::dbl_latt_storemode = dbl_latt_storemode;

  if ( !do_arg.Decode(argv[i++],"do_arg") ) { 
    do_arg.Encode("bum_arg","bum_arg");
    printf("Bum do_arg\n"); 
    exit(-1);
  }

  if ( !hmc_arg.Decode(argv[i++],"hmc_arg")){printf("Bum hmc_arg\n"); exit(-1);}
  if ( !evo_arg.Decode(argv[i++],"evo_arg")){printf("Bum evo_arg\n"); exit(-1);}
  if ( !gauge_arg.Decode(argv[i++],"gauge_arg")){printf("Bum gauge_arg\n"); exit(-1);}
  if ( !ab1_arg.Decode(argv[i++],"ab1_arg")){printf("Bum ab1_arg\n"); exit(-1);}

  if(load_config && !disable_ustar_reconstruction){ //use regular lattice load on first creation
    do_arg.start_conf_kind = START_CONF_FILE;
    do_arg.start_conf_filename = load_config_file;
    load_config = false;
  }
  if(do_arg.start_conf_kind == START_CONF_FILE && disable_ustar_reconstruction){ //take control of load process
    do_arg.start_conf_kind = START_CONF_ORD;
    load_config = true;
    load_config_file = do_arg.start_conf_filename;
  }
  
#ifdef USE_SCU_CHECKSUMS
  ScuChecksum::Initialise(evo_arg.hdw_xcsum,evo_arg.hdw_rcsum);
#endif


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

  if( load_config && disable_ustar_reconstruction){    
    //load the config now where we can set the disable_reconstruction bit
    if(!UniqueID()) printf("Loading config with reconstruction bit disabled\n");
    Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
    ReadLatticeParallel readLat;
    readLat.setSerial();
    readLat.disableGparityReconstructUstarField();
    readLat.read(lat,load_config_file);
    LatticeFactory::Destroy();
  }

  Complex wline[4] CPS_FLOAT_ALIGN;
    wline[0] = 0.0;
    wline[1] = 0.0;
    wline[2] = 0.0;
    wline[3] = 0.0;
  Float temp = UniqueID()*8;

  // Outer config loop
  int traj = evo_arg.traj_start;

  hmc_arg_pass = hmc_arg;

  //!< Create fictitous Hamiltonian (mom + action)
  AlgMomentum mom;
  AlgActionGauge gauge(mom, gauge_arg);
  
  //!< Construct numerical integrators
  AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge, ab1_arg);
  
  for(int conf=0; conf< evo_arg.gauge_configurations; conf ++ ) {

    sprintf(plaq_file,"%s.%d",evo_arg.plaquette_stem,traj);
    FILE * truncate_it = Fopen(plaq_file,"w");
    Fclose(truncate_it);
    common_arg_plaq.set_filename(plaq_file);

    sprintf(wline_file,"%s_wline.%d",evo_arg.plaquette_stem,traj);
    truncate_it = Fopen(wline_file,"w");
    Fclose(truncate_it);
    common_arg_wline.set_filename(wline_file);

    sprintf(hmc_file,"%s.%d",evo_arg.evo_stem,traj);
    common_arg_hmc.set_filename(hmc_file);

    LRGState rng_state;
    // Inner trajectory loop
    for(int i=0;i<evo_arg.gauge_unload_period;i++,traj++ ) {
      {
	// Wilson gauge action used for plaquette measurement
	Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
	AlgPlaq plaq(lat, &common_arg_plaq, &no_arg);
	plaq.run();
	AlgWline wline2(lat, &common_arg_wline, &no_arg);
	wline2.run();
	LatticeFactory::Destroy();
      }    


      {
	if ( (evo_arg.reproduce_interval > 0) &&
	     ((traj+1) % evo_arg.reproduce_interval) == 0 ) { 
	  VRB.Result("","main()","Running traj %d with reproduction\n",traj);
	  hmc_arg_pass.reproduce = REPRODUCE_YES;
	} else { 
	  VRB.Result("","main()","Running traj %d without reproduction\n",traj);
	  hmc_arg_pass.reproduce = REPRODUCE_NO;
	}

	//!< Run hybrid Monte Carlo
	AlgHmc hmc(ab1, common_arg_hmc, hmc_arg_pass);
	Float time = -dclock();
	hmc.run();
	time += dclock();
	print_flops("AlgHmc","run()",0,time);
      }

#ifdef USE_SCU_CHECKSUMS
      if ( ! ScuChecksum::CsumSwap() ) { 
	fprintf(stderr, "Checksum mismatch\n");
	exit(-1);
      }
#endif

    }//End of inter-cfg sweep

    if(traj == test_conf_reproduce) checkReproduces(test_conf_reproduce_file, disable_ustar_reconstruction);
    checkpoint(traj);

  } //End config loop

  AlgIntAB::Destroy(ab1);

  End();

 return(0);
}

void checkpoint(int traj)
{
  char *cname="cps::";
  char *fname="checkpoint()";

  char lat_file[256];
  char rng_file[256];
  
  Float time = -dclock();
  // Save this config to disk
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

  sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
  QioArg wt_arg(lat_file,0.001);
    
  wt_arg.ConcurIONumber=evo_arg.io_concurrency;
  WriteLatticeParallel wl;
  wl.setSerial();
  wl.setHeader(evo_arg.ensemble_id,evo_arg.ensemble_label,traj);
  wl.write(lat,wt_arg);
     
  if(!wl.good()) 
    ERR.General(cname,fname,"Failed write lattice %s",lat_file);

  LatticeFactory::Destroy();
  
  // Save the RNG's
  LRG.setSerial();
  
  sprintf(rng_file,"%s.%d",evo_arg.rng_file_stem,traj);
  if ( !LRG.Write(rng_file) ) 
    ERR.General(cname,fname,"Failed RNG file %s",rng_file);
  
  // Update the parameter files for restart

  do_arg.start_seed_filename = rng_file;
  do_arg.start_seed_kind = START_SEED_FILE;
  do_arg.start_conf_filename = lat_file;
  do_arg.start_conf_kind = START_CONF_FILE;
  evo_arg.traj_start     = traj;

  char vml_file[256];
  sprintf(vml_file,"do_arg.%d",traj);
  if ( !do_arg.Encode(vml_file,"do_arg") ){
    printf("bad do_arg encode\n");
    exit(-1);
  }

  sprintf(vml_file,"hmc_arg.%d",traj);
  if ( !hmc_arg.Encode(vml_file,"hmc_arg") ){
    printf("bad hmc_arg encode\n");
    exit(-1);
  }

  sprintf(vml_file,"evo_arg.%d",traj); 
  if ( !evo_arg.Encode(vml_file,"evo_arg")){
    printf("bad evo_arg encode\n");
    exit(-1);
  }

  sprintf(vml_file,"gauge_arg.%d",traj);
  if ( !gauge_arg.Encode(vml_file,"gauge_arg") ){
    printf("bad gauge_arg encode\n");
    exit(-1);
  }

  sprintf(vml_file,"ab1_arg.%d",traj);
  if ( !ab1_arg.Encode(vml_file,"ab1_arg") ){
    printf("bad ab1_arg encode\n");
    exit(-1);
  }

  time += dclock();
  print_flops("","checkpoint()",0,time);

}


void checkReproduces(const char *file, const bool &disable_ustar_reconstruction, const Float &tolerance){
  //first copy the existing lattice data
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

  Float plaq = lat.SumReTrPlaq() / 18.0 / GJP.VolSites() ;
  if(!UniqueID()) printf("Generated lattice plaquette: %f\n",plaq);  

  int array_size = 2*lat.GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
  Matrix *orig_lattice = (Matrix *) pmalloc(array_size);
  memcpy((void*)orig_lattice, (void*)lat.GaugeField(), array_size);

  //load the config that we want to compare it to
  {
    ReadLatticeParallel readLat;
    readLat.setSerial();
    if(disable_ustar_reconstruction) readLat.disableGparityReconstructUstarField();
    readLat.read(lat,file);
  }
  
  Matrix *loaded_lattice = lat.GaugeField();
  bool fail(false);

  for(int i=0;i<2*lat.GsiteSize() * GJP.VolNodeSites();i++){
    Float *orig = (Float*)orig_lattice + i;
    Float *loaded = (Float*)loaded_lattice + i;

    if( fabs( *orig - *loaded ) > tolerance ){
      printf("Node %d, Reproduction test failed at i=%d: %f %f\n",UniqueID(),i,*orig,*loaded);
      fail=true;
    }
  }
  LatticeFactory::Destroy();

  if(fail) exit(-1);
  else printf("Reproduce test passed\n");
}
