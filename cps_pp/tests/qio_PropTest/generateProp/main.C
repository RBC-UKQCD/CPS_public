

// v4.1


// what to do and what not...

# define VOLFMT QIO_PARTFILE
//undef, QIO_SINGLEFILE QIO_PARTFILE

#include <stdio.h>
#include <stdlib.h>


#include <util/lattice.h>
#include <alg/do_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <util/command_line.h>
#include <sstream>

#include <util/qio_readLattice.h>
#include <util/qio_writeLattice.h>

#include<unistd.h>
#include<config.h>

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include <util/qio_readPropagator.h>
#include <util/qio_writePropagator.h>

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

using namespace std;
using namespace cps;


int main(int argc,char *argv[])
{

  cps::Start(&argc, &argv);

  CommandLine::is(argc,argv);

  DoArg do_arg;
  EvoArg evo_arg;
  QPropWArg qpropw_arg;


  /* get parameter from command line */

  /* 6 parameters: */
  /* x.x  
          DIRECTORY 
	  do_arg evo_arg qpropw_arg
	  randomgauge
	  config_out
  */

 
  chdir(CommandLine::arg());
  
  
  if ( !do_arg.Decode(CommandLine::arg(),"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}
  if ( !evo_arg.Decode(CommandLine::arg(),"evo_arg") ) { printf("Bum evo_arg\n"); exit(-1);}
  if ( !qpropw_arg.Decode(CommandLine::arg(),"qpropq_arg") ) { printf("Bum qpropw_arg\n"); exit(-1);}


  GJP.Initialize(do_arg);

  int gauge_ran(CommandLine::arg_as_int() );

  char out_lat[200];
  char out_qio[200];
  char out_std[200];
  
  sprintf(out_lat, CommandLine::arg() );

  sprintf(out_std,"%s.IEEE64BIG",out_lat);
  sprintf(out_qio,"%s.QIO",out_lat);
  
  GimprRectFdwf lattice;

  

  
  printf("\ncreate random gauge - field (taking # %i) \n",gauge_ran);
 
  for(int ii(0); ii < gauge_ran; ++ii)
    lattice.SetGfieldDisOrd();

  WriteLatticeParallel writeLat;

  int traj(1203);

  writeLat.setHeader(evo_arg.ensemble_id, evo_arg.ensemble_label, traj);
  
  printf("  writing: %s (NERSC-format)\n",out_std);
  
  /* write lattice */
  
  writeLat.write(lattice, out_std, FP_IEEE64BIG);
  
  {
    
    qio_writeLattice writeLatQio;
        
    writeLatQio.setHeader(evo_arg.ensemble_id, evo_arg.ensemble_label, traj);
    
    printf("  writing: %s (QIO-format, double)\n",out_qio);
      
#ifdef VOLFMT
    writeLatQio.write(out_qio, "added ildgLFN", lattice, VOLFMT);
#else
    
    writeLatQio.write( out_qio, "added ildgLFN", lattice);
#endif
  }

  CommonArg common_arg("common_file", "common_label");    

  printf("   **** calculating propagator ****   \n");
  
  QPropWPointSrc propagator(lattice, &qpropw_arg, &common_arg);

  printf("    **** propagator ready ***** \n");


}





