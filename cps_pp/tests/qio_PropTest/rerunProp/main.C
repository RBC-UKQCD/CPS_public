

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

using namespace std;
using namespace cps;



int main(int argc,char *argv[])
{

  cps::Start(&argc,&argv);

  CommandLine::is(argc,argv);

  DoArg do_arg;
  EvoArg evo_arg;
  QPropWArg qpropw_arg;


  /* get parameter from command line */

  /* 6 or 7 parameters: */
  /* x.x  
          DIRECTORY 
	  do_arg evo_arg qpropw_arg 
	  config_type(0=NERSC, 1=QIO) read_config 
	  [ precision ] 
  */

 
  chdir(CommandLine::arg());
  
  
  if ( !do_arg.Decode(CommandLine::arg(),"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}
  if ( !evo_arg.Decode(CommandLine::arg(),"evo_arg") ) { printf("Bum evo_arg\n"); exit(-1);}
  if ( !qpropw_arg.Decode(CommandLine::arg(),"qpropw_arg") ) { printf("Bum qpropw_arg\n"); exit(-1);}

  GJP.Initialize(do_arg);

  int lat_type = CommandLine::arg_as_int();

  char in_lat[200];

  sprintf(in_lat, CommandLine::arg() );
  
  int setPrec(0);
  Float prec(1.);

  if(argc > 7){
    setPrec=1;
    prec=CommandLine::arg_as_Float();
  }

  GimprRectFdwf lattice;
  

  if(lat_type){

    qio_readLattice readLatQio;
    
    printf("  reading: %s (QIO-format)\n",in_lat);
        
#ifdef VOLFMT
    readLatQio.read(in_lat, lattice, VOLFMT);
#else
    readLatQio.read(in_lat, lattice);
#endif
    
  }
  else{

    ReadLatticeParallel readLat;
    
    printf("  reading: %s (NERSC-format)\n",in_lat);
        
    readLat.read(lattice, in_lat);
    
  }
  
  printf("   **** re-calculating propagator ****   \n");
  
  CommonArg common_arg("common_file", "common_label");

  QPropWPointSrc propagator(lattice, &common_arg);

  propagator.SetArgs(qpropw_arg);
	
  if( setPrec){
    printf("  using prec. %e\n",prec);
    propagator.ReRun(prec);
  }
  else{
    printf("  using standard prec.\n");
    propagator.ReRun();
  }

  printf("   **** propagator ready ***** \n");
	


}





