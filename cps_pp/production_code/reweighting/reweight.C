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


  if ( argc < 5 ) { 
  if(!UniqueID()){ 
    printf("Args:\tdo_arg.vml quo_arg.vml current_dir (start traj)\n");
    for(int i =0; i<argc; i++){
      printf("argv[%d]=%s\n",i,argv[i]);
    }
    exit(-1);
  }
  }

  Chdir (argv[4]);

  if ( !do_arg.Decode(argv[1],"do_arg") ) { 
    do_arg.Encode("bum_arg","bum_arg");
    printf("Bum do_arg\n"); 
    exit(-1);
  }

  if ( !evo_arg.Decode(argv[2],"evo_arg")){printf("Bum evo_arg\n"); exit(-1);}
  if ( !quo_arg.Decode(argv[3],"quo_arg")){printf("Bum quo_arg\n"); exit(-1);}
  if ( !quo_arg.Encode("quo_arg.out","quo_arg")){printf("Bum quo_arg.out\n"); exit(-1);}


  Chdir(evo_arg.work_directory);


  GJP.Initialize(do_arg);
  LRG.Initialize();


  // Outer config loop
  int traj = evo_arg.traj_start;
  if (argc>5) sscanf(argv[5],"%d",&traj);

//  int w_int = evo_arg.measure_w_spect_interval;
  int g_int = evo_arg.gauge_unload_period;
  
  hmc_arg_pass = hmc_arg;
  sum_arg.A_steps = 1;
  sum_arg.B_steps = 1;

  //!< Create fictitous Hamiltonian (mom + action)
  AlgMomentum mom;
//  AlgActionGauge gauge(mom, gauge_arg);
//    AlgActionRationalQuotient rat_quo(mom, rat_quo_arg);

  char lat_file[256];
  char rng_file[256];
  char rw_file[256];
  char vml_file[256];

  int n_step = evo_arg.measure_pbp;
  int n_mass = quo_arg.quotients.quotients_len;
  // for now assumes all masses for quotients are the same
  Float upper = quo_arg.quotients.quotients_val[0].bsn_mass;
  Float lower = quo_arg.quotients.quotients_val[0].frm_mass;
  VRB.Result(cname,fname,"bsn_mass=%g frm_mass=%g n_step=%d n_mass=%d\n",
       upper,lower,n_step,n_mass);
  Float *rw_fac,*norm;
  rw_fac = new Float[n_mass];
  norm = new Float[n_mass];
  
  for(int conf=0; conf< evo_arg.gauge_configurations; conf ++ ) {
    sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
    sprintf(rw_file,"reweight.%d",traj);

    sprintf(vml_file,"quo_arg.%d",traj);
    if ( !quo_arg.Encode(vml_file,"quo_arg") ){
          printf("bad rat_quo_arg encode\n");
          exit(-1);
    }
    FILE *fp = Fopen(rw_file,"a");
//    sprintf(rng_file,"%s.%d",evo_arg.rng_file_stem,traj);
    {
      GwilsonFdwf lat;
      ReadLatticeParallel(lat,lat_file);
    }
//    LRG.Read(rng_file,32);

    for(int i =0;i<n_step;i++){
      Float b_mass = upper - (Float)(i)/(Float)n_step *(upper-lower);
      Float f_mass = upper - (Float)(i+1)/(Float)n_step *(upper-lower);
      for(int j = 0;j<n_mass;j++){
        quo_arg.quotients.quotients_val[j].bsn_mass  = b_mass;
        quo_arg.quotients.quotients_val[j].frm_mass  = f_mass;
      }
      AlgActionQuotient quotient(mom, quo_arg);
      quotient.reweight(rw_fac,norm);
      VRB.Result(cname,fname,"%s:rw_fac=%e norm=%e diff =%e\n",lat_file,*rw_fac,*norm,*rw_fac-*norm);
      Fprintf(fp," %e\t%e\t%0.14e\t%0.14e\n",b_mass,f_mass,*rw_fac,*norm);
    }
    Fclose(fp);
    traj += g_int;
    LRG.Write("reweight_rng.out",32);
  } //End config loop
  delete[] rw_fac;
  delete[] norm;


  End();

 return(0);
}
