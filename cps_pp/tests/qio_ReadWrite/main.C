

// v4.1


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

#define DO_prop
 #undef DO_readPropQIO
 #define DO_writePropQIO
  #define DO_rereadPropQIO
   #define DO_compareProp   



// if defined writes sg prec. prop.
#define WRITE_propSG

#define WRITE_propTypeA
#define WRITE_propTypeB
#define WRITE_propTypeC
// A= scalarSource + 12 sink
// B= source sink pairs
// C= scalarSource sink pairs


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

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include <util/qio_readPropagator.h>
#include <util/qio_writePropagator.h>

using namespace std;
using namespace cps;

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
  QPropWArg qpropw_arg;


  /* get parameter from command line */

  /* 10 parameters: */
  /* x.x  
          DIRECTORY do_arg evo_arg 
          read_file write_file (for standard)
	  read_file write_file (for QIO)
	  number, which random gauge to take
	  read_prop write_prop
 */

 
  chdir(CommandLine::arg());
  
  
  if ( !do_arg.Decode(CommandLine::arg(),"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}
  if ( !evo_arg.Decode(CommandLine::arg(),"evo_arg") ) { printf("Bum evo_arg\n"); exit(-1);}

  GJP.Initialize(do_arg);



  char infile[200], outfile[200];
  char in_qio[200], out_qio[200];
  char out_prop[200], in_prop[200]; 

  sprintf(infile, CommandLine::arg() );
  sprintf(outfile, CommandLine::arg() );
  
  sprintf(in_qio, CommandLine::arg() );
  sprintf(out_qio, CommandLine::arg() );

  int gauge_ran(CommandLine::arg_as_int() );

  //if ( !qpropw_arg.Decode(CommandLine::arg(),"qpropw_arg") ) { printf("Bum qpropw_arg\n"); exit(-1);} 

  sprintf(in_prop, CommandLine::arg() );
  sprintf(out_prop, CommandLine::arg() );


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



#ifdef DO_prop

      CommonArg common_arg("my_file", "my_label");    



      QPropWPointSrc propagator_read(lattice, &common_arg);
      propagator_read.Allocate(0);


#ifdef WRITE_propTypeA 

      char out_prop_a[220];
      sprintf(out_prop_a,"scS_12sink.%s",out_prop);

      Float dummy_scSource_write[GJP.VolNodeSites()][2];
      
      const int scSourceSize(2*GJP.VolNodeSites());
      
      Float* dummy_scSource_read = (Float*)smalloc(scSourceSize*sizeof(Float));

      for(int ii(0); ii < GJP.VolNodeSites(); ++ii)
	{ 
	  dummy_scSource_write[ii][0] = Float(ii);
	  dummy_scSource_write[ii][1] = -Float(ii);
	  
	}
#endif


#ifdef WRITE_propTypeC 

      char out_prop_c[220];
      sprintf(out_prop_c,"scS_sink_pairs.%s",out_prop);

      Float dummy_scSourcePairs_write[GJP.VolNodeSites()][12][2];

      const int scSourcePairsSize(2*12*GJP.VolNodeSites());
      
      Float* dummy_scSourcePairs_read = (Float*)smalloc(scSourcePairsSize*sizeof(Float));


      
      for(int ii(0); ii < GJP.VolNodeSites(); ++ii)
	for(int jj(0); jj < 12; ++jj)
	  { 
	    dummy_scSourcePairs_write[ii][jj][0] = Float(ii);
	    dummy_scSourcePairs_write[ii][jj][1] = Float(jj);
	  }
	   

#endif

#ifdef WRITE_propTypeB 

      char out_prop_b[220];
      sprintf(out_prop_b,"fS_sink_pairs.%s",out_prop);

      Float dummy_fSourcePairs_write[GJP.VolNodeSites()][12][12][2];

      const int fSourcePairsSize(2*12*12*GJP.VolNodeSites());

      Float* dummy_fSourcePairs_read = (Float*)smalloc(fSourcePairsSize*sizeof(Float));



      for(int ii(0); ii < GJP.VolNodeSites(); ++ii)
	for(int jj(0); jj < 12; ++jj)
	  for(int kk(0); kk < 12; ++kk)
	    { 
	      dummy_fSourcePairs_write[ii][jj][kk][0] = Float(ii);
	      dummy_fSourcePairs_write[ii][jj][kk][1] = Float(kk)+12*Float(jj);
	    }
	

#endif

      #ifdef DO_readPropQIO

      {

#error NOT READY: readPropQIO

        qio_readPropagator readPropQio(argc,argv);

	printf("  reading: %s  (QIO-format)\n",in_prop);


        #ifdef VOLFMT
        readPropQio.read(in_prop, &propagator_read[0], dummy_arbSource_read, 12,  12, 12,VOLFMT);
        #else

        readPropQio.read(in_prop, &propagator_read[0], dummy_arbSource_read, 12, 12, 12 );
        #endif
      }

	#endif // DO_readPropQIO


      #ifdef DO_writePropQIO

	printf("   **** calculating propagator ****   \n");




	qpropw_arg.file="my_cg_prop_file";
	qpropw_arg.x=0;
	qpropw_arg.y=0;
	qpropw_arg.z=0;
	qpropw_arg.t=0;
	qpropw_arg.gauge_fix_src=0;
	qpropw_arg.gauge_fix_snk=0;
	qpropw_arg.store_midprop=0;
	qpropw_arg.save_prop=0;
	qpropw_arg.do_half_fermion=0;	

	qpropw_arg.cg.mass=0.5;
	qpropw_arg.cg.max_num_iter=5000;
	qpropw_arg.cg.stop_rsd=1e-5;
	qpropw_arg.cg.true_rsd=1e-5;
	qpropw_arg.cg.RitzMatOper=MAT_HERM;
	qpropw_arg.cg.Inverter=CG;
	qpropw_arg.cg.bicgstab_n=0;

	QPropWPointSrc propagator(lattice, &qpropw_arg, &common_arg);


	printf("   **** propagator ready ***** \n");


	printf("   ---- starting QIO-part ----- \n");

      

        {

          qio_writePropagator writePropQio(argc, argv);



        int traj(1203);
        writePropQio.setHeader(evo_arg.ensemble_id, evo_arg.ensemble_label, traj);

	printf("  writing: %s (QIO-format, double)\n",out_prop);


	#ifdef WRITE_propTypeA 

	printf("  using format ScalarSource + 12 sinks\n");

	#ifndef WRITE_propSG

        #ifdef VOLFMT
        writePropQio.write_ScS_12sink(out_prop_a, &propagator[0], dummy_scSource_write, VOLFMT);
        #else
	writePropQIO.write_ScS_12sink( out_prop_a,&propagator[0], dummy_scSource_write);
        #endif

	#else

	 #ifdef VOLFMT
        writePropQio.write_ScS_12sink(out_prop_a, &propagator[0], dummy_scSource_write, VOLFMT, FP_IEEE32);
        #else
	writePropQIO.write_ScS_12sink( out_prop_a,&propagator[0], dummy_scSource_write, QIO_PARTFILE, FP_IEEE32 );
        #endif

	#endif

	#endif

	#ifdef WRITE_propTypeB 

	printf("  using format source/sink pairs\n");

	#ifndef WRITE_propSG


        #ifdef VOLFMT
        writePropQio.write_12pairs(out_prop_b, QIO_FULL_SOURCE, &propagator[0], dummy_fSourcePairs_write, VOLFMT);
        #else
	writePropQIO.write_12pairs(out_prop_b, QIO_FULL_SOURCE, &propagator[0], dummy_fSourcePairs_write);
        #endif

	#else

	 #ifdef VOLFMT
        writePropQio.write_12pairs(out_prop_b, QIO_FULL_SOURCE, &propagator[0], dummy_fSourcePairs_write, VOLFMT, FP_IEEE32);
        #else
	writePropQIO.write_12pairs(out_prop_b, QIO_FULL_SOURCE, &propagator[0], dummy_fSourcePairs_write, QIO_PARTFILE,  FP_IEEE32);
        #endif

	#endif


	#endif
	
	#ifdef WRITE_propTypeC

	printf("  using format scalarSource/sink pairs\n");

	#ifndef WRITE_propSG

        #ifdef VOLFMT
        writePropQio.write_12pairs(out_prop_c, QIO_SCALAR_SOURCE, &propagator[0], dummy_scSourcePairs_write, VOLFMT);
        #else
	writePropQIO.write_12pairs(out_prop_c, QIO_SCALAR_SOURCE, &propagator[0], dummy_scSourcePairs_write);
        #endif

	#else

	#ifdef VOLFMT
        writePropQio.write_12pairs(out_prop_c, QIO_SCALAR_SOURCE, &propagator[0], dummy_scSourcePairs_write, VOLFMT, FP_IEEE32);
        #else
	writePropQIO.write_12pairs(out_prop_c, QIO_SCALAR_SOURCE, &propagator[0], dummy_scSourcePairs_write, QIO_PARTFILE, FP_IEEE32);
        #endif

	#endif

	#endif




        }



      #ifdef DO_rereadPropQIO


	
	{

        qio_readPropagator readPropQio(argc,argv);

	printf("  re-reading: %s  (QIO-format)\n",out_prop);

#ifdef WRITE_propTypeA 

	printf("  using format ScalarSource + 12 sinks\n");


        #ifdef VOLFMT
        readPropQio.read_ScS_12sink(out_prop_a, &propagator_read[0], dummy_scSource_read, VOLFMT);
        #else

        readPropQio.read_ScS_12sink(out_prop_a, &propagator_read[0], dummy_scSource_read);
        #endif



#ifdef DO_compareProp
	{
	  Float errCount(0);

	  Float prec(1e-10);

	  #ifdef WRITE_propSG
	  prec = 1e-5;
	  #endif

	  WilsonMatrix matrix_calc, matrix_read;

	  Complex tmp_c, tmp_r;

	  Float diff, diff_source;

	  for(int ii(0); ii < GJP.VolNodeSites(); ++ii)
	    {
	      matrix_calc = propagator[ii];
	      matrix_read = propagator_read[ii];
	      	    
	      
	      Float *source_read_0 = (Float*) dummy_scSource_read + 2*ii;
	      Float *source_read_1 = (Float*) dummy_scSource_read + 2*ii + 1;

	      diff_source = fabs(dummy_scSource_write[ii][0] - *source_read_0) + fabs(dummy_scSource_write[ii][1] - *source_read_1);



	      if(diff_source > prec)
		{
		  errCount += 1.0;
		  printf("mismatch source: index %i: (%f,%f) (%f,%f) \n",ii,*source_read_0,*source_read_1,dummy_scSource_write[ii][0],dummy_scSource_write[ii][1]);
		}

	    

	      for(int s1(0); s1 < 4; ++s1)
		for(int c1(0); c1 < 3; ++c1)
		  for(int s2(0); s2 < 4; ++s2)
		    for(int c2(0); c2 < 3; ++c2)
		      {
			tmp_c = matrix_calc(s1,c1,s2,c2);
			tmp_r = matrix_read(s1,c1,s2,c2);
			
			diff = fabs(tmp_c.real() - tmp_r.real() ) + fabs(tmp_c.imag() - tmp_r.imag() );
		      
			if ( diff > prec)
			  {
			    errCount += 1.0;
			    printf("mismatch: index %i, %i %i %i %i: (%f,%f) (%f,%f) \n",ii,s1,c1,s2,c2,tmp_c.real(),tmp_c.imag(),tmp_r.real(),tmp_r.imag());
			  }
		      }
	    

	    }
	  
	  glb_sum_five(&errCount);
	
	  printf(" total error-count: %f\n", errCount);

	}

#endif // DO_compareProp





#endif // write A

#ifdef WRITE_propTypeB 

	printf("  using 12 fullSource sink pairs\n");


        #ifdef VOLFMT
        readPropQio.read_12pairs(out_prop_b, &propagator_read[0], dummy_fSourcePairs_read, QIO_FULL_SOURCE, VOLFMT);
        #else

        readPropQio.read_12pairs(out_prop_b, &propagator_read[0], dummy_fSourcePairs_read, QIO_FULL_SOURCE);
        #endif


#ifdef DO_compareProp
	{
	  Float errCount(0);

	  Float prec(1e-10);

	  #ifdef WRITE_propSG
	  prec = 1e-5;
	  #endif

	  WilsonMatrix matrix_calc, matrix_read;

	  Complex tmp_c, tmp_r;

	  Float diff, diff_source;

	  for(int ii(0); ii < GJP.VolNodeSites(); ++ii)
	    {
	      matrix_calc = propagator[ii];
	      matrix_read = propagator_read[ii];
		  
	      for(int jj(0); jj < 12; ++jj)
		for(int kk(0); kk <12; ++kk)
		  {

		    
		    Float *source_read_0 = (Float*) dummy_fSourcePairs_read + 288*ii + 24*jj + 2*kk;
		    Float *source_read_1 = (Float*) dummy_fSourcePairs_read + 288*ii + 24*jj + 2*kk + 1;

		    diff_source = fabs(dummy_fSourcePairs_write[ii][jj][kk][0] - *source_read_0) + fabs(dummy_fSourcePairs_write[ii][jj][kk][1] - *source_read_1);
		    
		    		    
		    
		    if(diff_source > prec)
		      {
			errCount += 1.0;
			printf("mismatch source: index %i %i %i: (%f,%f) (%f,%f) \n",ii, jj,kk,
			       *source_read_0,*source_read_1,dummy_fSourcePairs_write[ii][jj][kk][0],dummy_fSourcePairs_write[ii][jj][kk][1]);
		      }
		  }
	    

	      for(int s1(0); s1 < 4; ++s1)
		for(int c1(0); c1 < 3; ++c1)
		  for(int s2(0); s2 < 4; ++s2)
		    for(int c2(0); c2 < 3; ++c2)
		      {
			tmp_c = matrix_calc(s1,c1,s2,c2);
			tmp_r = matrix_read(s1,c1,s2,c2);
			
			diff = fabs(tmp_c.real() - tmp_r.real() ) + fabs(tmp_c.imag() - tmp_r.imag() );
		      
			if ( diff > prec)
			  {
			    errCount += 1.0;
			    printf("mismatch: index %i, %i %i %i %i: (%f,%f) (%f,%f) \n",ii,s1,c1,s2,c2,tmp_c.real(),tmp_c.imag(),tmp_r.real(),tmp_r.imag());
			  }
		      }
	    

	    }
	  
	  glb_sum_five(&errCount);
	
	  printf(" total error-count: %f\n", errCount);

	}

#endif // DO_compareProp


#endif // write B

#ifdef WRITE_propTypeC 

	printf("  using 12 scalarSource sink pairs\n");


        #ifdef VOLFMT
        readPropQio.read_12pairs(out_prop_c, &propagator_read[0], dummy_scSourcePairs_read, QIO_SCALAR_SOURCE, VOLFMT);
        #else

        readPropQio.read_12pairs(out_prop_c, &propagator_read[0], dummy_scSourcePairs_read, QIO_SCALAR_SOURCE);
        #endif


#ifdef DO_compareProp
	{
	  Float errCount(0);

	  Float prec(1e-10);

	  #ifdef WRITE_propSG
	  prec = 1e-5;
	  #endif

	  WilsonMatrix matrix_calc, matrix_read;

	  Complex tmp_c, tmp_r;

	  Float diff, diff_source;

	  for(int ii(0); ii < GJP.VolNodeSites(); ++ii)
	    {
	      matrix_calc = propagator[ii];
	      matrix_read = propagator_read[ii];
	      
	    

	      for(int jj(0); jj < 12; ++jj)
		{


		  
		  Float *source_read_0 = (Float*) dummy_scSourcePairs_read + 24*ii + 2*jj ;
		  Float *source_read_1 = (Float*) dummy_scSourcePairs_read + 24*ii + 2*jj  + 1;

		  diff_source = fabs(dummy_scSourcePairs_write[ii][jj][0] - *source_read_0) + fabs(dummy_scSourcePairs_write[ii][jj][1] - *source_read_1);

		    
		  if(diff_source > prec)
		    {
		      errCount += 1.0;
		      printf("mismatch source: index %i %i : (%f,%f) (%f,%f) \n",ii, jj,
			     *source_read_0,*source_read_1,dummy_scSourcePairs_write[ii][jj][0],dummy_scSourcePairs_write[ii][jj][1]);
		    }
		}
	    

	      for(int s1(0); s1 < 4; ++s1)
		for(int c1(0); c1 < 3; ++c1)
		  for(int s2(0); s2 < 4; ++s2)
		    for(int c2(0); c2 < 3; ++c2)
		      {
			tmp_c = matrix_calc(s1,c1,s2,c2);
			tmp_r = matrix_read(s1,c1,s2,c2);
			
			diff = fabs(tmp_c.real() - tmp_r.real() ) + fabs(tmp_c.imag() - tmp_r.imag() );
		      
			if ( diff > prec)
			  {
			    errCount += 1.0;
			    printf("mismatch: index %i, %i %i %i %i: (%f,%f) (%f,%f) \n",ii,s1,c1,s2,c2,tmp_c.real(),tmp_c.imag(),tmp_r.real(),tmp_r.imag());
			  }
		      }
	    

	    }
	  
	  glb_sum_five(&errCount);
	
	  printf(" total error-count: %f\n", errCount);

	}

#endif // DO_compareProp



#endif // wirte C



      }



	#endif // DO_rereadPropQIO

      #endif // DO_writePropQIO

#endif //DO_prop

}





