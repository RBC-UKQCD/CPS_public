

/*  v2.0:

moving to CPS-implementation

* Done:
  ======
     

* ToDo:
  ======

     

*/


// what to do and what not...

# define VOLFMT QIO_PARTFILE   
//undef, QIO_SINGLEFILE QIO_PARTFILE

#undef DO_readStandard
#undef DO_writeStandard
#undef DO_readQIO
#define DO_randomGauge
#define DO_writeQIO
#undef DO_writeQIOsingle
#define DO_rereadQIO
 #undef DO_scramble
 #define DO_ordered
#undef DO_measure
#define DO_measure_compare


#include <config.h>
#include <stdio.h>
#include <stdlib.h>

#include <util/lattice.h>
#include <alg/do_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>

#if(0==1)
 #include <ReadLattice.h>
 #include <WriteLattice.h>
#endif

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <util/command_line.h>
#include <sstream>

#include <util/qio_readLattice.h>
#include <util/qio_writeLattice.h>

#include<unistd.h>
#include<config.h>



USING_NAMESPACE_CPS

using namespace std;


#if(0==1)
 GlobalJobParameter GJP;
 Verbose VRB;
 Error ERR;
#endif



int main(int argc,char *argv[])
{


  CommandLine::is(argc,argv);

  DoArg do_arg;
  EvoArg evo_arg;

  /* get parameter from command line */

  /* 8 parameters: */
  /* x.x  
          DIRECTORY do_arg evo_arg 
          read_file write_file (for standard)
	  read_file write_file (for QIO)
	  number, which random gauge to take
 */

 
  chdir(CommandLine::arg());

  
  
  if ( !do_arg.Decode(CommandLine::arg(),"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}
  if ( !evo_arg.Decode(CommandLine::arg(),"evo_arg") ) { printf("Bum evo_arg\n"); exit(-1);}

  GJP.Initialize(do_arg);



  char infile[200], outfile[200];
  char in_qio[200], out_qio[200];
  

  sprintf(infile, CommandLine::arg() );
  sprintf(outfile, CommandLine::arg() );
  
  sprintf(in_qio, CommandLine::arg() );
  sprintf(out_qio, CommandLine::arg() );

  int gauge_ran(CommandLine::arg_as_int() );

  const FP_FORMAT outformat ( 	FP_IEEE64BIG);

			     /* one of
				FP_UNKNOWN
				FP_AUTOMATIC
				FP_TIDSP32
				FP_IEEE32
				FP_IEEE32BIG
				FP_IEEE32LITTLE
				FP_IEEE64
				FP_IEEE64BIG
				FP_IEEE64LITTLE
			     */

#ifdef DO_measure_compare
   Float globPlaq(0), globPlaqOld(0), globPlaqDiff(0);
   Float nodePlaq(0), nodePlaqOld(0), nodePlaqDiff(0);
   int  compPlaq(0);
#endif // DO_measure_compare
  
  GimprRectFdwf lattice;
  //GimprRectFnone lattice;


 



#ifdef DO_readStandard

      ReadLatticeParallel readLat;

      printf("  reading: %s (NERSC-format)\n",infile); 
      
      /* read lattice */

      readLat.read(lattice, infile);

      #ifdef DO_measure_compare

      globPlaq=lattice.SumReTrPlaq();
      nodePlaq=lattice.SumReTrPlaqNode();

      printf("\n measured global plaquette: %f (UID: %i local: %f)\n",globPlaq,UniqueID(),nodePlaq);

      if( compPlaq)
	{
	  globPlaqDiff = globPlaq - globPlaqOld;
	  nodePlaqDiff = nodePlaq - nodePlaqOld;

	  printf("    difference: %f (UID: %i local: %f)\n", globPlaqDiff, UniqueID(), nodePlaqDiff);
	}

      globPlaqOld = globPlaq;
      nodePlaqOld = nodePlaq;
      compPlaq = 1;

      #endif //DO_measure_compare


      #ifdef DO_measure

      printf("\n measured global plaquette: %f (UID: %i local: %f)\n",lattice.SumReTrPlaq(),UniqueID(),lattice.SumReTrPlaqNode());      

      #endif // DO_measure

#endif // DO_readStandard



#ifdef DO_writeStandard

      WriteLatticeParallel writeLat;

      printf("  writing: %s (NERSC-format)\n",outfile);

      /* write lattice */

      writeLat.write(lattice, outfile, outformat);


#endif // DO_writeStandard






      
      printf("   **** starting QIO-part ****   \n");


      #ifdef DO_readQIO

      {

	qio_readLattice readLatQio( argc, argv);

      printf("  reading: %s (QIO-format)\n",in_qio);

      //qio_readField( in_qio, lattice, argc, argv);

	#ifdef VOLFMT
	readLatQio.read(in_qio, lattice, VOLFMT);
	#else
      	readLatQio.read(in_qio, lattice);
      	#endif

      }

      #ifdef DO_measure_compare

      //compPlaq = 0;

	 globPlaq=lattice.SumReTrPlaq();
	 nodePlaq=lattice.SumReTrPlaqNode();

	 printf("\n measured global plaquette: %f (UID: %i local: %f)\n",globPlaq,UniqueID(),nodePlaq);

	 if( compPlaq)
	   {
	     globPlaqDiff = globPlaq - globPlaqOld;
	     nodePlaqDiff = nodePlaq - nodePlaqOld;
	     
	     printf("    difference: %f (UID: %i local: %f)\n", globPlaqDiff, UniqueID(), nodePlaqDiff);
	   }

	 globPlaqOld = globPlaq;
	 nodePlaqOld = nodePlaq;
	 compPlaq = 1;

      #endif //DO_measure_compare


      #ifdef DO_measure

      printf("\n measured global plaquette: %f (UID: %i local: %f)\n",lattice.SumReTrPlaq(),UniqueID(),lattice.SumReTrPlaqNode());      

      #endif // DO_measure



      #endif // DO_readQIO
      

      #ifdef DO_randomGauge

	printf("\ncreate random gauge - field (taking # %i) \n",gauge_ran);
 
	for(int ii(0); ii < gauge_ran; ++ii)
	lattice.SetGfieldDisOrd();

	#ifdef DO_measure_compare
         
	   compPlaq = 0;
	
           globPlaq=lattice.SumReTrPlaq();
	   nodePlaq=lattice.SumReTrPlaqNode();

	   printf("\n measured global plaquette: %f (UID: %i local: %f)\n",globPlaq,UniqueID(),nodePlaq);

	   if( compPlaq)
	     {
	       globPlaqDiff = globPlaq - globPlaqOld;
	       nodePlaqDiff = nodePlaq - nodePlaqOld;

	       printf("    difference: %f (UID: %i local: %f)\n", globPlaqDiff, UniqueID(), nodePlaqDiff);
	     }

	   globPlaqOld = globPlaq;
	   nodePlaqOld = nodePlaq;
	   compPlaq = 1;

        #endif //DO_measure_compare
	
        #ifdef DO_measure

	 printf("\n measured global plaquette: %f (UID: %i local: %f)\n",lattice.SumReTrPlaq(),UniqueID(),lattice.SumReTrPlaqNode());      

        #endif // DO_measure	

      #endif // DO_randomGauge





      #ifdef DO_writeQIO

	{

	  qio_writeLattice writeLatQio(argc, argv);



	int traj(1203);
	writeLatQio.setHeader(evo_arg.ensemble_id, evo_arg.ensemble_label, traj);

      printf("  writing: %s (QIO-format, double)\n",out_qio);

      //qio_writeField( out_qio, lattice, argc, argv);

	#ifdef VOLFMT
        writeLatQio.write(out_qio, "added ildgLFN", lattice, VOLFMT);
        #else

      writeLatQio.write( out_qio, "added ildgLFN", lattice);
	#endif
	}

      #endif // DO_writeQIO

      #ifdef DO_writeQIOsingle

      {
	qio_writeLattice writeLatQio(argc, argv);

      printf("  writing: %s (QIO-format, single)\n",out_qio);

      //qio_writeField( out_qio, lattice, argc, argv, FP_IEEE32);

	#ifdef VOLFMT
      writeLatQio.write( out_qio, "added ildgLFN", lattice, VOLFMT, FP_IEEE32);
	#else
	writeLatQio.write( out_qio, "added ildgLFN", lattice,QIO_PARTFILE,  FP_IEEE32);
	#endif
      }

      #endif // DO_writeQIOsingle


      #ifdef DO_rereadQIO

       #ifdef DO_scramble
        printf("  scramble  gaugefield before re-reading it...\n");
        lattice.SetGfieldDisOrd();
       #endif //DO_scramble

        #ifdef DO_ordered
        printf("  order(=1) gaugefield before re-reading it...\n");
        lattice.SetGfieldOrd();
       #endif //DO_ordered

      #ifdef DO_measure
        printf("\n measured global plaquette: %f (UID: %i local: %f)\n",lattice.SumReTrPlaq(),UniqueID(),lattice.SumReTrPlaqNode());      
      #endif // DO_measure

      {

	qio_readLattice readLatQio(argc,argv);

      printf("  re-reading: %s  (QIO-format)\n",out_qio);

      //qio_readField( out_qio, lattice, argc, argv);
      
	#ifdef VOLFMT
        readLatQio.read(out_qio, lattice, VOLFMT);
        #else

	readLatQio.read(out_qio, lattice);
	#endif
      }

      #ifdef DO_measure_compare

        globPlaq=lattice.SumReTrPlaq();
	nodePlaq=lattice.SumReTrPlaqNode();

	printf("\n measured global plaquette: %f (UID: %i local: %f)\n",globPlaq,UniqueID(),nodePlaq);

	if( compPlaq)
	  {
	    globPlaqDiff = globPlaq - globPlaqOld;
	    nodePlaqDiff = nodePlaq - nodePlaqOld;
	    
	    printf("    difference: %f (UID: %i local: %f)\n", globPlaqDiff, UniqueID(), nodePlaqDiff);
	  }

	globPlaqOld = globPlaq;
	nodePlaqOld = nodePlaq;
	compPlaq = 1;

      #endif //DO_measure_compare


      #ifdef DO_measure

      printf("\n measured global plaquette: %f (UID: %i local: %f)\n",lattice.SumReTrPlaq(),UniqueID(),lattice.SumReTrPlaqNode());      

      #endif // DO_measure

   
      printf("  values: %s %s \n",in_qio, out_qio);



      #endif // DO_rereadQIO

      printf("   ---- leaving QIO-part  ----   \n");








}


