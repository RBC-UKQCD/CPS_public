
#define CHECKSOURCE
#define VOLFMT QIO_PARTFILE 
// QIO_PARTFILE, QIO_SINGLEFILE


#include <stdio.h>
#include <stdlib.h>

#include <util/lattice.h>
#include <alg/do_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>

#include <util/command_line.h>
#include <sstream>

#include<unistd.h>
#include<config.h>

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include <util/qio_readPropagator.h>
#include <util/qio_writePropagator.h>

using namespace std;
using namespace cps;

const Float PRECISION = 5e-5;

Float SpinColorIndex( const int spin_src, const int color_src, const int spin_snk, const int color_snk)
{

  Float tmp(0.0);

  tmp = color_src + 3*spin_src + 12*color_snk + 36*spin_snk;
  tmp /= 1000;
  

  return tmp;
}

int SiteIndex( const int mysite)
{
  int tmp(0);
  int vol(1);

  Site site(mysite);

  for(int ii(0); ii <4; ++ii){
    tmp += vol*site.physCoor(ii);
    vol *= GJP.NodeSites(ii)*GJP.Nodes(ii);
  }

  return tmp;
 

}


int main(int argc,char *argv[])
{

  cps::Start(&argc, &argv);

  CommandLine::is(argc,argv);

  DoArg do_arg;

  /* get parameter from command line */

  /* 5 or 6 parameters: */
  /* x.x  
          DIRECTORY 
	  do_arg 
 	  read_prop write_prop
	  mode (r, w, wr=rw)
	  [ precision ]
 */

 
  chdir(CommandLine::arg());
  
  
  if ( !do_arg.Decode(CommandLine::arg(),"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}

  GJP.Initialize(do_arg);

  char out_prop[200], in_prop[200], mode[10]; 

  sprintf(in_prop, CommandLine::arg() );
  sprintf(out_prop, CommandLine::arg() );
  sprintf(mode, CommandLine::arg());

  Float prec(PRECISION);

  if(argc > 6)
    prec = CommandLine::arg_as_Float();

  GimprRectFdwf lattice;
 
  CommonArg common_arg("my_file", "my_label");    

  QPropWPointSrc propagator_read(lattice, &common_arg);
  propagator_read.Allocate(0);


  char out_prop_b[220], out_prop_b_sg[220];
  sprintf(out_prop_b,"fS_sink_pairs_db.%s",out_prop);
  sprintf(out_prop_b_sg,"fS_sink_pairs_sg.%s",out_prop);

  Float fSourcePairs_write[GJP.VolNodeSites()][4][3][4][3][2];
  Float propagator[GJP.VolNodeSites()][4][3][4][3][2];

  const int fSourcePairsSize(2*12*12*GJP.VolNodeSites());

  Float* fSourcePairs_read = (Float*)smalloc(fSourcePairsSize*sizeof(Float));

  for(int site(0); site < GJP.VolNodeSites(); ++site){

    int mySite = SiteIndex(site);

    for(int s_src(0); s_src < 4; ++s_src)
      for(int c_src(0); c_src < 3; ++c_src)
	for(int s_snk(0); s_snk < 4; ++s_snk)
	  for(int c_snk(0); c_snk < 3; ++c_snk){
      
	    Float myIndex = SpinColorIndex(s_src,c_src,s_snk,c_snk);

	    fSourcePairs_write[site][s_snk][c_snk][s_src][c_src][0] = mySite + myIndex;
	    fSourcePairs_write[site][s_snk][c_snk][s_src][c_src][1] = mySite - myIndex;

	    propagator[site][s_snk][c_snk][s_src][c_src][0] = mySite + myIndex;
	    propagator[site][s_snk][c_snk][s_src][c_src][1] = mySite - myIndex;

	  }
  }


  

  printf("   ---- starting QIO-part ----- \n");
	
  


  if( !(strcmp(mode,"w")) || !(strcmp(mode,"rw")) || !(strcmp(mode,"wr")) ) {
    
    qio_writePropagator writePropQIO;
    
    int traj(1203);
    writePropQIO.setHeader("DUMMY_ID", "DEBUG_PROPAGATOR", traj);
    

    printf("  using format source/sink pairs\n");
    printf("  writing: %s (QIO-format, double)\n",out_prop_b);
    
    writePropQIO.write_12pairs(out_prop_b, QIO_FULL_SOURCE, &propagator[0], fSourcePairs_write, VOLFMT);
    
    printf("  writing: %s (QIO-format, single)\n",out_prop_b_sg);
    writePropQIO.write_12pairs(out_prop_b_sg, QIO_FULL_SOURCE, &propagator[0], fSourcePairs_write, VOLFMT,  FP_IEEE32);
    
  }	



  if( !(strcmp(mode,"r")) || !(strcmp(mode,"rw")) || !(strcmp(mode,"wr")) ) {
    
    qio_readPropagator readPropQio;
    
    printf("  using 12 fullSource sink pairs\n");
    printf("  reading: %s (QIO-format)\n",in_prop);
    
    for(int ii(0); ii < fSourcePairsSize; ++ii)
      *(fSourcePairs_read + ii) = 0.0;
    
    {
      Float *zeroIt = (Float *) &propagator_read[0];
      for(int ii(0); ii < 288*GJP.VolNodeSites(); ++ii)
	*(zeroIt + ii) = 0.0;
    }
    
    readPropQio.read_12pairs(in_prop, &propagator_read[0], fSourcePairs_read, QIO_FULL_SOURCE);
    
    {
      Float errCount(0);

      WilsonMatrix matrix_read;
      
      Complex tmp_r;
      Float tmp_c_r, tmp_c_i;
      
      Float diff, diff_source;
      
      
      printf(" checking read propagator with prec. %e\n",prec);
      
      for(int site(0); site < GJP.VolNodeSites(); ++site){
	
	matrix_read = propagator_read[site];
	
	for(int s_src(0); s_src < 4; ++s_src)
	  for(int c_src(0); c_src < 3; ++c_src)
	    for(int s_snk(0); s_snk < 4; ++s_snk)
	      for(int c_snk(0); c_snk < 3; ++c_snk){
		
#ifdef CHECKSOURCE				    
		  Float *source_read_0 = (Float*) fSourcePairs_read + 288*site + 72*s_snk + 24*c_snk  + 6*s_src + 2*c_src;
		  Float *source_read_1 = (Float*) fSourcePairs_read + 288*site + 72*s_snk + 24*c_snk  + 6*s_src + 2*c_src+ 1;

		  diff_source = fabs(fSourcePairs_write[site][s_snk][c_snk][s_src][c_src][0] - *source_read_0) + fabs(fSourcePairs_write[site][s_snk][c_snk][s_src][c_src][1] - *source_read_1);
		    		    
		  if(diff_source > prec){
		    errCount += 1.0;
		    printf("mismatch source: index %i, %i %i %i %i: (%f,%f) (%f,%f) \n",site, s_snk, c_snk, s_src, c_src,
			   *source_read_0,*source_read_1,fSourcePairs_write[site][s_snk][c_snk][s_src][c_src][0],fSourcePairs_write[site][s_snk][c_snk][s_src][c_src][1]);
		  }
	
#endif // CHECKSOURCE

		  tmp_c_r = propagator[site][s_snk][c_snk][s_src][c_src][0];
		  tmp_c_i = propagator[site][s_snk][c_snk][s_src][c_src][1];
		  
		  tmp_r = matrix_read(s_snk,c_snk,s_src,c_src);
		  
		  diff = fabs(tmp_c_r - tmp_r.real() ) + fabs(tmp_c_i - tmp_r.imag() );
		  
		  if ( diff > prec){
		    errCount += 1.0;
		    printf("mismatch prop: index %i, %i %i %i %i: (%f,%f) (%f,%f) \n",site,s_snk,c_snk,s_src,c_src,tmp_r.real(),tmp_r.imag(), tmp_c_r, tmp_c_i);
		  }
	      }
	
	
      }
      
      glb_sum_five(&errCount);
      
      printf(" total error-count: %f\n", errCount);
      
    }





  }
}




