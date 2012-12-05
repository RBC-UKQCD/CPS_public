/***********************************************************
 * PAB: simplified main programme for off line measurements
 *      using VML for parameter loading
 *  Arguments
 *     DoArg ParameterFile in node-path
 *     MeasArg ParameterFile in node-path
 *     Initial workdir for reading the directories
 *
 *****************************************************
 */
#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>

#include<alg/common_arg.h>
#include<alg/do_arg.h>

#include<alg/meas_arg.h>
#include<alg/alg_meas.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>

//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

void ReadGauge(MeasArg *meas);
void WriteGauge(MeasArg *meas);

void ReadRNG(MeasArg *meas);
void WriteRNG(MeasArg *meas);

int main(int argc, char *argv[])
{ 
  char *cname=argv[0];
  char *fname="main()";

  CommonArg common_arg;
  DoArg do_arg;
  MeasArg meas_arg;

  if ( argc!=4 ) { 
    printf("Args: doarg-file measarg-file initial-directory\n");
    exit(-1);
  }

  chdir (argv[3]);

  if ( !do_arg.Decode(argv[1],"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}

  /*Layout the lattice on the machine (without regard to even-odd)*/
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

  if ( !meas_arg.Decode(argv[2],"meas_arg")){printf("Bum meas_arg\n"); exit(-1);}

  chdir (meas_arg.WorkDirectory);

#ifdef USE_SCU_CHECKSUMS
  ScuChecksum::Initialise(meas_arg.HdwXCsum,meas_arg.HdwRCsum);
#endif

  do_arg.verbose_level=VERBOSE_RESULT_LEVEL;
  GJP.Initialize(do_arg);
  VRB.Level(VERBOSE_RESULT_LEVEL);
  LRG.Initialize();
 
  { /*Force a lattice instantiation to read stuff*/
    GnoneFnone lat;
  }

  /************************************************
   * config loop
   ************************************************/
  for(int conf=meas_arg.TrajStart; 
      conf< meas_arg.TrajLessThanLimit; 
      conf +=meas_arg.TrajIncrement ) {

    meas_arg.TrajCur = conf;

    ReadGauge(&meas_arg);
    ReadRNG  (&meas_arg);

    { 
      AlgMeas meas(&common_arg,&meas_arg);
      meas.run();
    }

    WriteGauge(&meas_arg);
    WriteRNG  (&meas_arg);

#ifdef USE_SCU_CHECKSUMS
        if ( ! ScuChecksum::CsumSwap() ) { 
	  fprintf(stderr, "Checksum mismatch\n");
	  exit(-1);
	}
#endif

  } /*End config loop*/
     
 return(0);
}
void ReadGauge(MeasArg *meas)
{
  char *cname = "cps";
  char *fname = "ReadGauge";
  if ( meas->GaugeIO == MeasIOLoad ) { 
    GnoneFnone lat;
    char lat_file[256];
    sprintf(lat_file,"%s.%d",meas->GaugeStem,meas->TrajCur);
    QioArg rd_arg(lat_file,0.001);
    rd_arg.ConcurIONumber=meas->IOconcurrency;

    ReadLatticeParallel rl;
    rl.read(lat,rd_arg);
    if(!rl.good()) 
      ERR.General(cname,fname,"Failed read lattice %s",lat_file);
    
  }
}
void WriteGauge(MeasArg *meas)
{ 
  char *cname = "cps";
  char *fname = "WriteGauge";
  if ( meas->GaugeIO == MeasIOSave ) { 
    GnoneFnone lat;
    char lat_file[256];
    sprintf(lat_file,"%s.%d",meas->GaugeStem,meas->TrajCur);
    QioArg wt_arg(lat_file,0.001);
    wt_arg.ConcurIONumber=meas->IOconcurrency;

    WriteLatticeParallel wl;
    wl.setHeader("MeasurementCodeNotForGaugeProd","Dummy",meas->TrajCur);
    wl.write(lat,wt_arg);
    if(!wl.good()) 
      ERR.General(cname,fname,"Failed write lattice %s",lat_file);
  }
}
void ReadRNG(MeasArg *meas)
{
  char *cname = "cps";
  char *fname = "ReadRNG";
  if ( meas->RNGIO == MeasIOLoad ) { 
    char rng_file[256];
    sprintf(rng_file,"%s.%d",meas->RNGStem,meas->TrajCur);
    if ( !LRG.Read(rng_file) ) 
      ERR.General(cname,fname,"Failed RNG file %s",rng_file);
  }
}
void WriteRNG(MeasArg *meas)
{
  char *cname = "cps";
  char *fname = "WriteRNG";
  if ( meas->RNGIO == MeasIOSave ) { 
    char rng_file[256];
    sprintf(rng_file,"%s.%d",meas->RNGStem,meas->TrajCur);
    if ( !LRG.Write(rng_file) ) 
      ERR.General(cname,fname,"Failed RNG file %s",rng_file);
  }
}

