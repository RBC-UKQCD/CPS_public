#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/random.h>
#if (VERSION_SUB >17 || VERSION_MAJOR > 4)
#include<util/time_cps.h>
#else
#include<util/time.h>
#endif

#include<alg/alg_hmc.h>
#include<alg/alg_plaq.h>
#include<alg/alg_rect.h>
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
#include<alg/alg_w_spect.h>
#include<alg/array_arg.h>
#include<alg/alg_fix_gauge.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>

#undef DO_SCU_CHECKSUMS
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
#include <qcdocos/scu_checksum.h>
#endif
//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

const int MAX_FILENAME = 256;

HmcArg hmc_arg;
HmcArg hmc_arg_pass;

ActionGaugeArg gauge_arg;
ActionQuotientArg quo_arg;
ActionRationalQuotientArg rat_quo_arg;
ActionRationalQuotientArg rat_quo_arg_new;

IntABArg ab1_arg;
IntABArg ab2_arg;
IntABArg ab3_arg;
IntABArg sum_arg;

EvoArg evo_arg;
DoArg do_arg;
PbpArg pbp_arg;
NoArg no_arg;
FloatArray w_spect_mass_list;

void checkpoint(int traj);

inline int Chdir(const char *dir){
  if (chdir(dir) !=0){
    printf("cannot change directory to %s\n",dir);
    exit(-1);
  }
  return 0;
}

int main(int argc, char *argv[])
{ 

//  char plaq_file[256];
//  char pbp_file[256];
//  char hmc_file[256];

  char *cname=argv[0];
  char *fname="main()";
  Float dtime;
  
  CommonArg common_arg_plaq;
  CommonArg common_arg_pbp;
  CommonArg common_arg_hmc;

  WspectArg w_spect_arg;
  CgArg w_spect_cg_arg;
  WspectOutput w_spect_output;
  FixGaugeArg w_spect_fg_arg;

  Start(&argc,&argv);


  if ( argc < 6 ) { 
  if(!UniqueID()){ 
    printf("Args:\tdo_arg.vml rat_quo_arg.vml gauge_arg.vml current_dir \n");
    exit(-1);
    for(int i =0; i<argc; i++){
      printf("argv[%d]=%s\n",i,argv[i]);
    }
  }
  }

  Chdir (argv[6]);

  if ( !do_arg.Decode(argv[1],"do_arg") ) { 
    do_arg.Encode("bum_arg","bum_arg");
    printf("Bum do_arg\n"); 
    exit(-1);
  }

  if ( !evo_arg.Decode(argv[2],"evo_arg")){printf("Bum evo_arg\n"); exit(-1);}
  if ( !rat_quo_arg.Decode(argv[3],"rat_quo_arg")){printf("Bum rat_quo_arg\n"); exit(-1);}
  if ( !gauge_arg.Decode(argv[4],"gauge_arg")){printf("Bum gauge_arg\n"); exit(-1);}
  if ( !pbp_arg.Decode(argv[5],"pbp_arg")){printf("Bum pbp_arg\n"); exit(-1);}



  Chdir(evo_arg.work_directory);

  GJP.Initialize(do_arg);
  LRG.Initialize();


  // Outer config loop
  int traj = evo_arg.traj_start;
  if (argc > 7 ) sscanf(argv[7],"%d",&traj);
  evo_arg.traj_start = traj;
//  int w_int = evo_arg.measure_w_spect_interval;
  int g_int = evo_arg.gauge_unload_period;
  
  hmc_arg_pass = hmc_arg;
  sum_arg.A_steps = 1;
  sum_arg.B_steps = 1;

  //!< Create fictitous Hamiltonian (mom + action)
  AlgMomentum mom;
  AlgActionGauge gauge(mom, gauge_arg);
//  AlgActionQuotient quotient(mom, quo_arg);
//    AlgActionRationalQuotient rat_quo(mom, rat_quo_arg);
  int n_step = rat_quo_arg.bsn_mass.bsn_mass_len;

  char lat_file[256];
  char rng_file[256];
  char gauge_file[256];
  char plaq_file[256];
  char rect_file[256];
  char vml_file[256];
  sprintf(gauge_file,"gauge.H");
  
  for(int conf=0; conf< evo_arg.gauge_configurations; conf ++ ) {
    sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
    sprintf(plaq_file,"plaq.%d",traj);
    sprintf(rect_file,"rect.%d",traj);

    {
      GwilsonFdwf lat;
      ReadLatticeParallel(lat,lat_file);
      NoArg no_arg;
      CommonArg common_arg;
      common_arg.set_filename(plaq_file);
      AlgPlaq plaq(lat,&common_arg,&no_arg);
      plaq.run();
      common_arg.set_filename(rect_file);
      AlgRect rect(lat,&common_arg,&no_arg);
      rect.run();

    }

    FILE *fp2 = Fopen(gauge_file,"a");
    Float S_g = gauge.energy();
    glb_sum(&S_g);
    Fprintf(fp2,"%d\t%0.14e\n",traj,S_g);
    Fclose(fp2);

//      LRG.Write("reweight_rng.out",32);
//    Fclose(fp);

    traj += g_int;
  } //End config loop


  End();

 return(0);
}
