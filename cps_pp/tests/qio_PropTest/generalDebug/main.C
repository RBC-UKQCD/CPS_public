
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

  /* 5 parameters: */
  /* x.x  
          DIRECTORY 
	  do_arg 
 	  write_prop
	  type (A,B,C)
	  precision (D,F)
 */

 
  chdir(CommandLine::arg());
  
  
  if ( !do_arg.Decode(CommandLine::arg(),"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}

  GJP.Initialize(do_arg);

  char out_prop[200], type, prec, out_prop_mod[200]; 

  sprintf(out_prop, CommandLine::arg() );
  sprintf(&type, CommandLine::arg());
  sprintf(&prec, CommandLine::arg());

  GimprRectFdwf lattice;
 
  Float scSource[GJP.VolNodeSites()][2];
  Float scSourcePairs[GJP.VolNodeSites()][4][3][2];
  Float fSourcePairs[GJP.VolNodeSites()][4][3][4][3][2];
  Float propagator[GJP.VolNodeSites()][4][3][4][3][2];
  
  for(int site(0); site < GJP.VolNodeSites(); ++site){

    int mySite = SiteIndex(site);

    scSource[site][0]= mySite;
    scSource[site][1]= mySite;
    

    for(int s_src(0); s_src < 4; ++s_src)
      for(int c_src(0); c_src < 3; ++c_src){

	Float myScalIndex = SpinColorIndex(s_src,c_src, 0,0);
	
	scSourcePairs[site][s_src][c_src][0] = mySite + myScalIndex;
	scSourcePairs[site][s_src][c_src][1] = mySite - myScalIndex;


	for(int s_snk(0); s_snk < 4; ++s_snk)
	  for(int c_snk(0); c_snk < 3; ++c_snk){
	    
	    Float myIndex = SpinColorIndex(s_src,c_src,s_snk,c_snk);

	    fSourcePairs[site][s_snk][c_snk][s_src][c_src][0] = mySite + myIndex;
	    fSourcePairs[site][s_snk][c_snk][s_src][c_src][1] = mySite - myIndex;

	    propagator[site][s_snk][c_snk][s_src][c_src][0] = mySite + myIndex;
	    propagator[site][s_snk][c_snk][s_src][c_src][1] = mySite - myIndex;

	  }
	
      }

  }

    

  printf("   ---- starting QIO-part ----- \n");
	
  {

    qio_writePropagator writePropQIO;
    
    int traj(1203);
    writePropQIO.setHeader("DUMMY_ID", "DEBUG_PROPAGATOR", traj);
    

    switch( type ){

    case 'A':

      printf("  using format scalar source + 12 sinks\n");

      if(!strcmp(&prec,"D")){
	sprintf(out_prop_mod,"%s.scS12sinks_db",out_prop);
	printf("  writing: %s (QIO-format, double)\n",out_prop_mod);
	writePropQIO.write_ScS_12sink(out_prop_mod, &propagator[0], scSource, VOLFMT);
      }
      else{
	sprintf(out_prop_mod,"%s.scS12sinks_sg",out_prop);
	printf("  writing: %s (QIO-format, single)\n",out_prop_mod);
	writePropQIO.write_ScS_12sink(out_prop_mod, &propagator[0], scSource, VOLFMT,  FP_IEEE32);
      }



      break;

    case 'B':

      printf("  using format scalar source/sink pairs\n");

      if(!strcmp(&prec,"D")){
	sprintf(out_prop_mod,"%s.scSsinkPairs_db",out_prop);
	printf("  writing: %s (QIO-format, double)\n",out_prop_mod);
	writePropQIO.write_12pairs(out_prop_mod, QIO_SCALAR_SOURCE, &propagator[0], scSourcePairs, VOLFMT);
      }
      else{
	sprintf(out_prop_mod,"%s.scSsinkPairs_sg",out_prop);
	printf("  writing: %s (QIO-format, single)\n",out_prop_mod);
	writePropQIO.write_12pairs(out_prop_mod, QIO_SCALAR_SOURCE, &propagator[0], scSourcePairs, VOLFMT,  FP_IEEE32);
      }

      break;

    case 'C':

      printf("  using format full source/sink pairs\n");

      if(!strcmp(&prec,"D")){
	sprintf(out_prop_mod,"%s.fSsinkPairs_db",out_prop);
	printf("  writing: %s (QIO-format, double)\n",out_prop_mod);
	writePropQIO.write_12pairs(out_prop_mod, QIO_FULL_SOURCE, &propagator[0], fSourcePairs, VOLFMT);
      }
      else{
	sprintf(out_prop_mod,"%s.fSsinkPairs_sg",out_prop);
	printf("  writing: %s (QIO-format, single)\n",out_prop_mod);
	writePropQIO.write_12pairs(out_prop_mod, QIO_FULL_SOURCE, &propagator[0], fSourcePairs, VOLFMT,  FP_IEEE32);
      }

      break;

    default:

      printf(" %s out of type-range\n (A=scalar source + 12sinks, B=12 pairs of scalar source/sink, C=12 pairs of full source/sink)\n",&type);


    }

    
  }	

}




