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
ActionRationalQuotientArg rat_quo_l_arg;
ActionRationalQuotientArg rat_quo_h_arg;

IntABArg ab1_arg;
IntABArg ab2_arg;
IntABArg ab3_arg;
IntABArg sum_arg;

EvoArg evo_arg;
DoArg do_arg;
PbpArg pbp_arg;
NoArg no_arg;

void checkpoint(int traj);

int main(int argc, char *argv[])
{ 
  char plaq_file[256];
  char pbp_file[256];
  char wline_file[256];
  char hmc_file[256];

  const char *cname=argv[0];
  const char *fname="main()";
  Float dtime;
  
  CommonArg common_arg_hmc;
  CommonArg common_arg_plaq;
  CommonArg common_arg_pbp;

  Start(&argc,&argv);
  if ( argc!=12 ) { 
    if(!UniqueID()){
      printf("Expected 11 args:\n");
      printf("Args:\tdo_arg.vml hmc_arg.vml evo_arg.vml rat_quo_l_arg.vml\n");
      printf("\trat_quo_h_arg.vml gauge_arg.vml\n");
      printf("\tab1_arg.vml ab2_arg.vml ab3_arg.vml pbp_arg.vml\n");
      printf("\tcurrent_dir \n");
      
      printf("\nGot %d:\n",argc-1);
      for(int i=1;i<argc;i++){
	printf("%d: %s\n",i,argv[i]);
      }      
    }
    exit(-1);
    
  }
  if ( !do_arg.Decode(argv[1],"do_arg") ) { 
    do_arg.Encode("bum_arg","bum_arg");
    printf("Bum do_arg\n"); 
    exit(-1);
  }


  if ( !hmc_arg.Decode(argv[2],"hmc_arg")){printf("Bum hmc_arg\n"); exit(-1);}
  if ( !evo_arg.Decode(argv[3],"evo_arg")){printf("Bum evo_arg\n"); exit(-1);}
  if ( !rat_quo_l_arg.Decode(argv[4],"rat_quo_l_arg")){printf("Bum rat_quo_l_arg\n"); exit(-1);}
  if ( !rat_quo_h_arg.Decode(argv[5],"rat_quo_h_arg")){printf("Bum rat_quo_h_arg\n"); exit(-1);}
  if ( !gauge_arg.Decode(argv[6],"gauge_arg")){printf("Bum gauge_arg\n"); exit(-1);}
  if ( !ab1_arg.Decode(argv[7],"ab1_arg")){printf("Bum ab1_arg\n"); exit(-1);}
  if ( !ab2_arg.Decode(argv[8],"ab2_arg")){printf("Bum ab2_arg\n"); exit(-1);}
  if ( !ab3_arg.Decode(argv[9],"ab3_arg")){printf("Bum ab3_arg\n"); exit(-1);}
  if ( !pbp_arg.Decode(argv[10],"pbp_arg")){printf("Bum pbp_arg\n"); exit(-1);}

#ifdef USE_SCU_CHECKSUMS
  ScuChecksum::Initialise(evo_arg.hdw_xcsum,evo_arg.hdw_rcsum);
#endif

#if TARGET == BGQ
  LRG.setSerial();
#endif

  GJP.Initialize(do_arg);
  LRG.Initialize();

// #if TARGET != BGQ
//   GJP.SetNthreads(8);
//   if(!UniqueID()) printf("Set number of threads to %d\n", GJP.Nthreads());
// #endif

  // Outer config loop
  int traj = evo_arg.traj_start;

  hmc_arg_pass = hmc_arg;
  sum_arg.A_steps = 1;
  sum_arg.B_steps = 1;

  //!< Create fictitous Hamiltonian (mom + action)
  AlgMomentum mom;
  AlgActionGauge gauge(mom, gauge_arg);
  AlgActionRationalQuotient rat_quo_l(mom, rat_quo_l_arg);
  AlgActionRationalQuotient rat_quo_h(mom, rat_quo_h_arg);
  
  //!< Construct numerical integrators
  AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge, ab1_arg);
  AlgIntAB &ab2 = AlgIntAB::Create(ab1, rat_quo_h, ab2_arg);
  AlgIntAB &ab3 = AlgIntAB::Create(ab2, rat_quo_l, ab3_arg);
  
  for(int conf=0; conf< evo_arg.gauge_configurations; conf ++ ) {

    sprintf(plaq_file,"%s.%d",evo_arg.plaquette_stem,traj);
    FILE * truncate_it = Fopen(plaq_file,"w");
    Fclose(truncate_it);
    common_arg_plaq.set_filename(plaq_file);

    sprintf(pbp_file,"%s.%d",evo_arg.pbp_stem,traj);
    truncate_it = Fopen(pbp_file,"w");
    Fclose(truncate_it);
    common_arg_pbp.set_filename(pbp_file);
    pbp_arg.src_u_s = 0;
    pbp_arg.src_l_s = GJP.SnodeSites() * GJP.Snodes() - 1;
    pbp_arg.snk_u_s = GJP.SnodeSites() * GJP.Snodes() - 1;
    pbp_arg.snk_l_s = 0;

    sprintf(hmc_file,"%s.%d",evo_arg.evo_stem,traj);
    common_arg_hmc.set_filename(hmc_file);

    LRGState rng_state;
    // Inner trajectory loop
    for(int i=0;i<evo_arg.gauge_unload_period;i++,traj++ ) {
      { 
	//Setup timestamp
	char timestamp_file[256];
	sprintf(timestamp_file,"%s/timestamp.%d",evo_arg.work_directory,traj);
	TimeStamp::set_file(timestamp_file);
	TimeStamp::stamp_incr("Starting traj %d",traj);
       
	{
	  TimeStamp::stamp_incr("Starting plaquette measurement");
	  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
	  AlgPlaq plaq(lat, &common_arg_plaq, &no_arg);
	  plaq.run();
 	  LatticeFactory::Destroy();
	  TimeStamp::decr_stamp("Finished plaquette measurement");
	}    
	if (evo_arg.measure_pbp){
	  TimeStamp::stamp_incr("Starting psibar-psi measurement");
	  VRB.Result("","main()","Running pbp\n");
	  rng_state.GetStates();
	  Lattice &lat = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_WILSON);
	  AlgPbp pbp(lat,&common_arg_pbp,&pbp_arg);
	  pbp.run();
	  LatticeFactory::Destroy();
	  rng_state.SetStates();
	  TimeStamp::decr_stamp("Finished psibar-psi measurement");
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
	  TimeStamp::stamp_incr("Starting HMC");
	  AlgHmc hmc(ab3, common_arg_hmc, hmc_arg_pass);
          Float time = -dclock();
	  hmc.run();
          time += dclock();
          print_flops("AlgHmc","run()",0,time);
	  TimeStamp::decr_stamp("Finished HMC");
	}

#ifdef USE_SCU_CHECKSUMS
        if ( ! ScuChecksum::CsumSwap() ) { 
	  fprintf(stderr, "Checksum mismatch\n");
	  exit(-1);
	}
#endif

	if(i<evo_arg.gauge_unload_period-1) TimeStamp::close_file();
      }

    }//End of inter-cfg sweep
    TimeStamp::stamp_incr("Checkpointing configuration %d",traj);
    checkpoint(traj);
    TimeStamp::decr_stamp("Finished checkpointing configuration %d",traj);
    TimeStamp::decr_stamp("Finished configuration");
    TimeStamp::close_file();
  } //End config loop

  AlgIntAB::Destroy(ab3);
  AlgIntAB::Destroy(ab2);
  AlgIntAB::Destroy(ab1);

  End();

 return(0);
}

void checkpoint(int traj)
{
  const char *cname="cps::";
  const char *fname="checkpoint()";

  char lat_file[256];
  char rng_file[256];
  
  Float time = -dclock();
  // Save this config to disk
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

  sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
  QioArg wt_arg(lat_file,0.001);
    
  wt_arg.ConcurIONumber=evo_arg.io_concurrency;
  WriteLatticeParallel wl;
#if TARGET == BGQ
  wl.setSerial();
#endif

  wl.setHeader(evo_arg.ensemble_id,evo_arg.ensemble_label,traj);
  wl.write(lat,wt_arg);
     
  if(!wl.good()) 
    ERR.General(cname,fname,"Failed write lattice %s",lat_file);

  LatticeFactory::Destroy();
  
  // Save the RNG's
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

  sprintf(vml_file,"rat_quo_l_arg.%d",traj);
  if ( !rat_quo_l_arg.Encode(vml_file,"rat_quo_l_arg") ){
    printf("bad rat_quo_l_arg encode\n");
    exit(-1);
  }

  sprintf(vml_file,"rat_quo_h_arg.%d",traj);
  if ( !rat_quo_h_arg.Encode(vml_file,"rat_quo_h_arg") ){
    printf("bad rat_quo_h_arg encode\n");
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

  sprintf(vml_file,"ab2_arg.%d",traj);
  if ( !ab2_arg.Encode(vml_file,"ab2_arg") ){
    printf("bad ab2_arg encode\n");
    exit(-1);
  }

  sprintf(vml_file,"ab3_arg.%d",traj);
  if ( !ab3_arg.Encode(vml_file,"ab3_arg") ){
    printf("bad ab3_arg encode\n");
    exit(-1);
  }

  time += dclock();
  print_flops("","checkpoint()",0,time);

}
