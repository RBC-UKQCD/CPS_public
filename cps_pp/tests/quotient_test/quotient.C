#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/random.h>

#include<alg/alg_hmc.h>
#include<alg/common_arg.h>
#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>

#include<alg/alg_int.h>
#include<alg/int_arg.h>

#include<alg/no_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_plaq.h>
#include<alg/alg_remez.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>

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
ActionQuotientArg quo_arg;
ActionRationalQuotientArg rat_quo_arg;

IntABArg ab1_arg;
IntABArg ab2_arg;
IntABArg sum_arg;

EvoArg evo_arg;
DoArg do_arg;
NoArg no_arg;

#if TARGET == QCDOC
extern "C" {
  void _mcleanup(void);
}
#endif

void checkpoint(int traj);

int main(int argc, char *argv[])
{ 

  char plaq_file[256];
  char hmc_file[256];

  char *cname=argv[0];
  char *fname="main()";
  
  CommonArg common_arg_hmc;
  CommonArg common_arg_plaq;

  if ( argc!=10 ) { 
    printf("Args: doarg-file hmcarg-file evoarg-file initial-directory\n");
    exit(-1);
  }

  chdir (argv[9]);

  if ( !do_arg.Decode(argv[1],"do_arg") ) { 
    do_arg.Encode("bum_arg","bum_arg");
    printf("Bum do_arg\n"); 
    exit(-1);
  }

  //Layout the lattice on the machine (without regard to even-odd)
  do_arg.x_nodes = SizeX(); 
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = SizeS();

  do_arg.x_node_sites = do_arg.x_sites/do_arg.x_nodes; 
  do_arg.y_node_sites = do_arg.y_sites/do_arg.y_nodes;
  do_arg.z_node_sites = do_arg.z_sites/do_arg.z_nodes;
  do_arg.t_node_sites = do_arg.t_sites/do_arg.t_nodes;
  do_arg.s_node_sites = do_arg.s_sites/do_arg.s_nodes;

  if (do_arg.x_sites!=do_arg.x_node_sites*do_arg.x_nodes) {printf("Lattice does not fit\n");exit(-1);}
  if (do_arg.y_sites!=do_arg.y_node_sites*do_arg.y_nodes) {printf("Lattice does not fit\n");exit(-1);}
  if (do_arg.z_sites!=do_arg.z_node_sites*do_arg.z_nodes) {printf("Lattice does not fit\n");exit(-1);}
  if (do_arg.t_sites!=do_arg.t_node_sites*do_arg.t_nodes) {printf("Lattice does not fit\n");exit(-1);}
  if (do_arg.s_sites!=do_arg.s_node_sites*do_arg.s_nodes) {printf("Lattice does not fit\n");exit(-1);}

  do_arg.SetupAsqTadU0(do_arg.u0);

  if ( !hmc_arg.Decode(argv[2],"hmc_arg")){printf("Bum hmc_arg\n"); exit(-1);}
  if ( !evo_arg.Decode(argv[3],"evo_arg")){printf("Bum evo_arg\n"); exit(-1);}
  if ( !quo_arg.Decode(argv[4],"quo_arg")){printf("Bum quo_arg\n"); exit(-1);}
  if ( !rat_quo_arg.Decode(argv[5],"rat_quo_arg")){printf("Bum rat_quo_arg\n"); exit(-1);}
  if ( !gauge_arg.Decode(argv[6],"gauge_arg")){printf("Bum gauge_arg\n"); exit(-1);}
  if ( !ab1_arg.Decode(argv[7],"ab1_arg")){printf("Bum ab1_arg\n"); exit(-1);}
  if ( !ab2_arg.Decode(argv[8],"ab2_arg")){printf("Bum ab2_arg\n"); exit(-1);}

  chdir(evo_arg.work_directory);

#ifdef USE_SCU_CHECKSUMS
  ScuChecksum::Initialise(evo_arg.hdw_xcsum,evo_arg.hdw_rcsum);
#endif

  // do_arg.verbose_level=VERBOSE_RESULT_LEVEL;
  GJP.Initialize(do_arg);
  // VRB.Level(VERBOSE_RESULT_LEVEL);
  LRG.Initialize();

  // Outer config loop
  int traj = evo_arg.traj_start;
  
  if ( do_arg.start_conf_kind != START_CONF_FILE ) {
    checkpoint(traj);
  }

  hmc_arg_pass = hmc_arg;
  sum_arg.A_steps = 1;
  sum_arg.B_steps = 1;

  //!< Create fictitous Hamiltonian (mom + action)
  AlgMomentum mom;
  AlgActionGauge gauge(mom, gauge_arg);
  AlgActionQuotient quotient(mom, quo_arg);
  AlgActionRationalQuotient rat_quotient(mom, rat_quo_arg);
  
  AlgIntSum fermions(quotient, rat_quotient, sum_arg);

  //!< Construct numerical integrators
  AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge, ab1_arg);
  AlgIntAB &ab2 = AlgIntAB::Create(ab1, fermions, ab2_arg);
  
  for(int conf=0; conf< evo_arg.gauge_configurations; conf ++ ) {

    sprintf(plaq_file,"%s.%d",evo_arg.plaquette_stem,traj);
    FILE * truncate_it = Fopen(plaq_file,"w");
    Fclose(truncate_it);
    common_arg_plaq.set_filename(plaq_file);

    sprintf(hmc_file,"%s.%d",evo_arg.evo_stem,traj);
    common_arg_hmc.set_filename(hmc_file);

    // Inner trajectory loop
    for(int i=0;i<evo_arg.gauge_unload_period;i++,traj++ ) {
      { 
	
	{
	  // Wilson gauge action used for plaquette measurement
	  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_WILSON);
	  AlgPlaq plaq(lat, &common_arg_plaq, &no_arg);
	  plaq.run();
 	  LatticeFactory::Destroy();
	}    

	{
	  if ( (evo_arg.reproduce_interval > 0) &&
	       (traj % evo_arg.reproduce_interval) == 0 ) { 
	    printf("Running traj %d with reproduction\n",traj);
	    hmc_arg_pass.reproduce = REPRODUCE_YES;
	  } else { 
	    printf("Running traj %d without reproduction\n",traj);
	    hmc_arg_pass.reproduce = REPRODUCE_NO;
	  }

	  //!< Run hybrid Monte Carlo
	  AlgHmc hmc(ab2, common_arg_hmc, hmc_arg_pass);
	  hmc.run();
	}

#ifdef USE_SCU_CHECKSUMS
        if ( ! ScuChecksum::CsumSwap() ) { 
	  fprintf(stderr, "Checksum mismatch\n");
	  exit(-1);
	}
#endif

      }

    }//End of inter-cfg sweep

    checkpoint(traj);

  } //End config loop

  AlgIntAB::Destroy(ab2);
  AlgIntAB::Destroy(ab1);

#if TARGET == QCDOC
  _mcleanup();
#endif

 return(0);
}

void checkpoint(int traj)
{
  char *cname="cps::";
  char *fname="checkpoint()";

  char lat_file[256];
  char rng_file[256];
  
  // Save this config to disk
  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

  sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
  QioArg wt_arg(lat_file,0.001);
    
  wt_arg.ConcurIONumber=evo_arg.io_concurrency;
  WriteLatticeParallel wl;
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
  sprintf(vml_file,"evo_arg.%d",traj); 
  if ( !evo_arg.Encode(vml_file,"evo_arg")){
    printf("bad evo_arg encode\n");
    exit(-1);
  }
  sprintf(vml_file,"hmc_arg.%d",traj);
  if ( !hmc_arg.Encode(vml_file,"hmc_arg") ){
    printf("bad hmc_arg encode\n");
    exit(-1);
  }

  sprintf(vml_file,"quo_arg.%d",traj);
  if ( !quo_arg.Encode(vml_file,"quo_arg") ){
    printf("bad quo_arg encode\n");
    exit(-1);
  }

  sprintf(vml_file,"rat_quo_arg.%d",traj);
  if ( !quo_arg.Encode(vml_file,"quo_arg") ){
    printf("bad quo_arg encode\n");
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

}
