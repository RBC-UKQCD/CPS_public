/**********************************************************
 * EES: main programme based on PAB's vml_measure/measure.C
 *        to demonstrate usage of qioPropagator read/write functions
 *
 *
 ***********************************************************/



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

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
#define USE_SCU_CHECKSUMS
#else
#undef USE_SCU_CHECKSUMS
#endif

#ifdef USE_SCU_CHECKSUMS
#include <qcdocos/scu_checksum.h>
#endif
//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

void ReadGauge(MeasArg *meas, string &readID, string &readLabel, int &readSeqNum);
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

  CPS_NAMESPACE::Start(&argc,&argv);

  chdir (argv[3]);

  if ( !do_arg.Decode(argv[1],"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}

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

    string readID, readLabel;
    int readSeqNum;


    ReadGauge(&meas_arg, readID, readLabel, readSeqNum);
    ReadRNG  (&meas_arg);

    // calculate a propagator and store it....
    {

      GimprRectFdwf lat;

      QPropWArg qpropw_arg;

      //qpropw_arg.file="my_cg_prop_file";
      sprintf(qpropw_arg.file,"propagator_mass%f.%d", qpropw_arg.cg.mass, readSeqNum);
      sprintf(qpropw_arg.ensemble_label, readLabel.c_str());
      sprintf(qpropw_arg.ensemble_id, readID.c_str());
      qpropw_arg.seqNum = readSeqNum;

      qpropw_arg.x=0;
      qpropw_arg.y=0;
      qpropw_arg.z=0;
      qpropw_arg.t=0;
      qpropw_arg.gauge_fix_src=0;
      qpropw_arg.gauge_fix_snk=0;
      qpropw_arg.store_midprop=0;
      qpropw_arg.save_prop=1;
      qpropw_arg.do_half_fermion=0;	
      
      qpropw_arg.cg.mass=0.5;
      qpropw_arg.cg.max_num_iter=5000;
      qpropw_arg.cg.stop_rsd=1e-5;
      qpropw_arg.cg.true_rsd=1e-5;
      qpropw_arg.cg.RitzMatOper=MAT_HERM;
      qpropw_arg.cg.Inverter=CG;
      qpropw_arg.cg.bicgstab_n=0;
      
      QPropWPointSrc propagator(lat, &common_arg);

      propagator.SetArgs(qpropw_arg);


      printf("\n running \'propagator.Run_saveQIO\'\n  saving propagator -> [work-dir] %s\n", qpropw_arg.file);

      //propagator.Run_saveQIO(outprop, readID.c_str(), readLabel.c_str(), readSeqNum, argc, argv);
      propagator.Run();



    }


    // what else is specified in the meas_arg...
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
void ReadGauge(MeasArg *meas, string &readID, string &readLabel, int &readSeqNum)
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
    
    readID = rl.getEnsembleId();
    readLabel = rl.getEnsembleLabel();
    readSeqNum = rl.getSequenceNumber();


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

