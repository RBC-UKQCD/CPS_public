
/***********************************************************
 * PAB: simplified main programme for DWF RHMC algorithm
 *      using VML for parameter loading
 *  Arguments
 *     DoArg ParameterFile in absolute node-path
 *     HmdArg ParameterFile in absolute node-path
 *     EvoArg ParameterFile in absolute node-path
 *
 *  Future:
 *     Do we pass a pole/residue file to the RHMC (do so through alg_hmd)
 *
 *     I've inserted a placeholder RHMCPoleResStatus class and a
 *     filename entry into the HmdArg class, with the intention of
 *     Mike tweaking the RHMC to use it. At this point  the PoleRes
 *     stuff could be dropped.
 *     
 *     We could use the hmd_arg.rhmc_poles_action to decide how
 *     to modify the HmdArg parameter file at the next checkpoint.
 *****************************************************
 */

#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/random.h>

#include<alg/alg_hmd.h>
#include<alg/common_arg.h>
#include<alg/hmd_arg.h>
#include<alg/no_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_plaq.h>
#include<alg/alg_remez.h>

#include<alg/alg_eig.h>
#include<alg/eig_arg.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>

#include<alg/alg_int.h>
#include<alg/alg_hmc.h>

#define USE_SCU_CHECKSUMS
#ifdef USE_SCU_CHECKSUMS
#include <qcdocos/scu_checksum.h>
#endif
//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

HmdArg hmd_arg;
HmdArg hmd_arg_pass;

void checkpoint(DoArg &do_arg, HmdArg &hmd_arg,EvoArg &evo_arg, int traj);

int main(int argc, char *argv[])
{ 
  char plaq_file[256];
  char hmd_file[256];

  char *cname=argv[0];
  char *fname="main()";
  
  CommonArg common_arg_hmdr;
  CommonArg common_arg_plaq;

  EvoArg evo_arg;
  EigArg eig_arg;
  DoArg do_arg;
  NoArg no_arg;

  if ( argc!=6 ) { 
    printf("Args: doarg-file hmdarg-file evoarg-file eigarg_file initial-directory\n");
    exit(-1);
  }

  chdir (argv[5]);

  if ( !do_arg.Decode(argv[1],"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}

  do_arg.x_nodes = SizeX();/*Layout the lattice on the machine (without regard to even-odd)*/
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

  if ( !hmd_arg.Decode(argv[2],"hmd_arg")){printf("Bum hmd_arg\n"); exit(-1);}
  if ( !evo_arg.Decode(argv[3],"evo_arg")){printf("Bum evo_arg\n"); exit(-1);}
  if ( !eig_arg.Decode(argv[4],"eig_arg")){printf("Bum eig_arg\n"); exit(-1);}

  chdir(evo_arg.work_directory);

#ifdef USE_SCU_CHECKSUMS
  ScuChecksum::Initialise(evo_arg.hdw_xcsum,evo_arg.hdw_rcsum);
#endif

  // do_arg.verbose_level=VERBOSE_RESULT_LEVEL;
  GJP.Initialize(do_arg);
  // VRB.Level(VERBOSE_RESULT_LEVEL);
  LRG.Initialize();

  /************************************************
   * Outer config loop
   ************************************************/
  int traj = evo_arg.traj_start;
  
  
  if ( do_arg.start_conf_kind != START_CONF_FILE )
  {
    /*
     * Expect hot or cold start to be saved
     * but don't overwrite a loaded file
     */
    checkpoint(do_arg, hmd_arg, evo_arg, traj);
  }

  hmd_arg_pass = hmd_arg;

  int g_size = GJP.VolNodeSites()*18*4;
  Matrix *mom = (Matrix*)smalloc(g_size*sizeof(Float),cname,fname,"mom");

  for(int conf=0; conf< evo_arg.gauge_configurations; conf ++ ) {

    sprintf(plaq_file,"%s.%d",evo_arg.plaquette_stem,traj);
    FILE * truncate_it = Fopen(plaq_file,"w");
    Fclose(truncate_it);
    common_arg_plaq.set_filename(plaq_file);


    sprintf(hmd_file,"%s.%d",evo_arg.evo_stem,traj);
    common_arg_hmdr.set_filename(hmd_file);

    /************************************************
     * Inner trajectory loop
     ************************************************/
    for(int i=0;i<evo_arg.gauge_unload_period;i++,traj++ ) {

      {
	Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
	printf("Running plaq\n");
	AlgPlaq plaq(lat,&common_arg_plaq,&no_arg);
	plaq.run();
	LatticeFactory::Destroy();
      }
      
      if ( (traj % 100 ) < evo_arg.reproduce_percent ) { 
	printf("Running traj %d with reproduction\n",traj);
	hmd_arg_pass.reproduce = 1;
      } else { 
	printf("Running traj %d without reproduction\n",traj);
	hmd_arg_pass.reproduce = 0;
      }
      
      AlgIntRHMC rhmc(F_CLASS_DWF, &hmd_arg, mom);
      AlgIntBoson boson(F_CLASS_DWF, &hmd_arg, mom);
      AlgIntGauge gauge(G_CLASS_IMPR_RECT, mom);
      AlgIntMom momentum(mom);
      
      AlgIntSum gauge_boson(&gauge, &boson);
      AlgIntOmelyan omelyan_inner(&momentum, &gauge_boson);
      AlgIntOmelyan omelyan_outer(&omelyan_inner, &rhmc, 6, 1);
      
      AlgHmc hmc(&omelyan_outer, &common_arg_hmdr, &hmd_arg);
      
      hmc.run();
      
      //AlgHmcRHMC rhmc(lat,&common_arg_hmdr,&hmd_arg_pass);
      //rhmc.run();
      
#ifdef USE_SCU_CHECKSUMS
      if ( ! ScuChecksum::CsumSwap() ) { 
	fprintf(stderr, "Checksum mismatch\n");
	exit(-1);
      }
#endif
      
    }/*End of inter-cfg sweep*/
    
    if ( evo_arg.CalcEig ) { 
      Lattice &lat = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_NONE);
      
      // Measure the lowest eigenvalue
      printf("Checking eigen-spectrum\n");
      {
        CommonArg ca_eig;
        char eig_file[256];

        sprintf(eig_file,"%s.%d",evo_arg.eig_lo_stem,traj);

        FILE *truncate_it = Fopen(eig_file,"w");
        Fclose(truncate_it);

        ca_eig.set_filename(eig_file);
        eig_arg.RitzMatOper = MATPCDAG_MATPC;

        AlgEig eig(lat,&ca_eig,&eig_arg);
        eig.run();

      }

      /*
      // Measure the highest eigenvalue
      {
        CommonArg ca_eig;
        char eig_file[256];
  
        sprintf(eig_file,"%s.%d",evo_arg.eig_hi_stem,traj);
  
        FILE * truncate_it = Fopen(eig_file,"w");
        Fclose(truncate_it);
  
        ca_eig.set_filename(eig_file);
        eig_arg.RitzMatOper = NEG_MATPCDAG_MATPC;
  
        AlgEig eig(lat,&ca_eig,&eig_arg);
        eig.run();

      }
      */
      LatticeFactory::Destroy();
    }

    checkpoint(do_arg, hmd_arg, evo_arg, traj);

  } /*End config loop*/

  sfree(mom, cname, fname, "mom");
     
  return(0);
}

void checkpoint(DoArg &do_arg, HmdArg &hmd_arg,EvoArg &evo_arg, int traj)
{
  char *cname="cps::";
  char *fname="checkpoint()";

  char lat_file[256];
  char rng_file[256];

  Lattice &lat = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);

  /***************************************
   * Save this config to disk
   ***************************************/
  sprintf(lat_file,"%s.%d",evo_arg.gauge_file_stem,traj);
  QioArg wt_arg(lat_file,0.001);

  wt_arg.ConcurIONumber=evo_arg.io_concurrency;
  WriteLatticeParallel wl;
  wl.setHeader(evo_arg.ensemble_id,evo_arg.ensemble_label,traj);
  wl.write(lat,wt_arg);
    
  if(!wl.good()) 
    ERR.General(cname,fname,"Failed write lattice %s",lat_file);
  
  /***************************************
   * Save the RNG's
   ***************************************/
  sprintf(rng_file,"%s.%d",evo_arg.rng_file_stem,traj);
  if ( !LRG.Write(rng_file) ) 
    ERR.General(cname,fname,"Failed RNG file %s",rng_file);
  
  /****************************************
   * Update the parameter files for restart
   ******************************************/
  
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
  sprintf(vml_file,"hmd_arg.%d",traj);
  if ( !hmd_arg.Encode(vml_file,"hmd_arg") ){
    printf("bad hmd_arg encode\n");
    exit(-1);
  }

  LatticeFactory::Destroy();
}
