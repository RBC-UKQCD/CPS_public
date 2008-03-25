
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
 	  read_prop 
	  src_float_perSite snk_float_perSite
	  [ precision ]
 */

 
  chdir(CommandLine::arg());
  
  
  if ( !do_arg.Decode(CommandLine::arg(),"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}

  GJP.Initialize(do_arg);

  char in_prop[200]; 

  sprintf(in_prop, CommandLine::arg() );
  
  int src_Float_perSite(0), snk_Float_perSite(0);
  
  src_Float_perSite = CommandLine::arg_as_int();
  snk_Float_perSite = CommandLine::arg_as_int();

  Float prec(PRECISION);

  if(argc > 6)
    prec = CommandLine::arg_as_Float();

  //GimprRectFdwf lattice;
 
  //CommonArg common_arg("my_file", "my_label");    

  //QPropWPointSrc propagator_read(lattice, &common_arg);
  //propagator_read.Allocate(0);

  //Float fSourcePairs_write[GJP.VolNodeSites()][4][3][4][3][2];
  //Float propagator[GJP.VolNodeSites()][4][3][4][3][2];

  Float *src_read = (Float*) smalloc(src_Float_perSite*GJP.VolNodeSites()*sizeof(Float));
  Float *snk_read = (Float*) smalloc(snk_Float_perSite*GJP.VolNodeSites()*sizeof(Float));

  printf("   ---- starting QIO-part ----- \n");
	
  int readSources(0), readProps(0);
  QIO_PROP_SOURCE_TYPES readSourceType; 

  {
    
    qio_readPropagator readPropQio;
    
    printf("  using general reader\n");
    printf("  reading: %s (QIO-format)\n",in_prop);
    
    readPropQio.read(in_prop, snk_read, src_read, snk_Float_perSite*GJP.VolNodeSites(), src_Float_perSite*GJP.VolNodeSites(), VOLFMT );

    readSources = readPropQio.readSources;
    readProps= readPropQio.readProps;
    readSourceType = readPropQio.readSourceType;

    printf("  read %i propagators\n and %i sources of type %i\n", readProps, readSources, readSourceType); 
    
  }

  Float errCount(0);
  
  WilsonMatrix *matrix_read;
  
  Complex tmp_r;
  Float tmp_c_r, tmp_c_i;
  
  Float diff, diff_source;
  
  int checkedSources(0);

  printf(" checking read propagator with prec. %e\n",prec);
  
  for(int site(0); site < GJP.VolNodeSites(); ++site){
    
    matrix_read = (WilsonMatrix*) snk_read + site;
    
    for(int s_src(0); s_src < 4; ++s_src)
      for(int c_src(0); c_src < 3; ++c_src)
	for(int s_snk(0); s_snk < 4; ++s_snk)
	  for(int c_snk(0); c_snk < 3; ++c_snk){
	    
#ifdef CHECKSOURCE	

	    if( checkedSources < readSources){

	      if( readSourceType == QIO_FULL_SOURCE)  {
			    
		Float *source_read_0 = (Float*) src_read + 288*site + 72*s_snk + 24*c_snk  + 6*s_src + 2*c_src;
		Float *source_read_1 = (Float*) src_read + 288*site + 72*s_snk + 24*c_snk  + 6*s_src + 2*c_src+ 1;
	    
		Float exp_0 = SiteIndex(site) + SpinColorIndex(s_src,c_src,s_snk,c_snk);
		Float exp_1 = SiteIndex(site) - SpinColorIndex(s_src,c_src,s_snk,c_snk);

		diff_source = fabs(exp_0 - *source_read_0) + fabs(exp_1 - *source_read_1);
		    		    
		if(diff_source > prec){
		  errCount += 1.0;
		  printf("mismatch source: index %i, %i %i %i %i: (%f,%f) (%f,%f) \n",site, s_snk, c_snk, s_src, c_src,
			 *source_read_0,*source_read_1,exp_0,exp_1);
		}
	      
		++checkedSources;
	      }
	      
	      if( (readSourceType == QIO_SCALAR_SOURCE) && (s_snk==0) && (c_snk==0) ) {
		
		Float *source_read_0 = (Float*) src_read + 24*site + 6*s_src + 2*c_src ;
		Float *source_read_1 = (Float*) src_read + 24*site + 6*s_src + 2*c_src + 1;
	    
		Float exp_0 = SiteIndex(site) + SpinColorIndex(s_src,c_src,0,0);
		Float exp_1 = SiteIndex(site) - SpinColorIndex(s_src,c_src,0,0);

		diff_source = fabs(exp_0 - *source_read_0) + fabs(exp_1 - *source_read_1);
		    		    
		if(diff_source > prec){
		  errCount += 1.0;
		  printf("mismatch source: index %i, %i %i %i %i: (%f,%f) (%f,%f) \n",site, s_snk, c_snk, s_src, c_src,
			 *source_read_0,*source_read_1,exp_0,exp_1);
		}
	      
		++checkedSources;
	      }
		  
	    }

#endif // CHECKSOURCE
	    
	    tmp_c_r = SiteIndex(site) + SpinColorIndex(s_src,c_src,s_snk,c_snk);//propagator[site][s_snk][c_snk][s_src][c_src][0];
	    tmp_c_i = SiteIndex(site) - SpinColorIndex(s_src,c_src,s_snk,c_snk);//propagator[site][s_snk][c_snk][s_src][c_src][1];
	    
	    tmp_r = (*matrix_read)(s_snk,c_snk,s_src,c_src);
	    
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











