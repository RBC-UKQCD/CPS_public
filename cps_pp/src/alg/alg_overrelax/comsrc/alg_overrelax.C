#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Definitions of the AlgOverRelax class methods.
  
  $Id: alg_overrelax.C,v 1.3 2004/10/27 10:26:09 mclark Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mclark $
//  $Date: 2004/10/27 10:26:09 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_overrelax/comsrc/alg_overrelax.C,v 1.3 2004/10/27 10:26:09 mclark Exp $
//  $Id: alg_overrelax.C,v 1.3 2004/10/27 10:26:09 mclark Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: alg_overrelax.C,v $
//  $Revision: 1.3 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_overrelax/comsrc/alg_overrelax.C,v $
//  $State: Exp $
//

CPS_END_NAMESPACE
#include <stdlib.h>	// exit()
#include <util/qcdio.h>
#include <math.h>
#include <time.h>
#include <alg/alg_overrelax.h>
#include <alg/common_arg.h>
#include <alg/overrelax_arg.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/random.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>

CPS_START_NAMESPACE

//Uncomment next line to switch on timing
//#define OVERRELAX_TIMING

//------------------------------------------------------------------
/*!
  \param latt The lattice on which to perform the heatbath
  \param c_arg The common argument structure for all algorithms.
  \param arg The parameters specific to this algorithm.
 */
//------------------------------------------------------------------
AlgOverRelax::AlgOverRelax(Lattice& latt, 
	     CommonArg *c_arg,
	     OverRelaxArg *arg) : 
	     Alg(latt, c_arg) 
{
  cname = "AlgOverRelax";
  char *fname = "AlgOverRelax(L&,CommonArg*,OverrelaxArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)   ERR.Pointer(cname,fname, "arg");
  alg_overrelax_arg = arg;
  fGamma = (2./3.) * GJP.Beta();

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgOverRelax::~AlgOverRelax() {
  //char *fname = "~AlgOverRelax()";
  //VRB.Func(cname,fname);
}

//------------------------------------------------------------------
// AlgOverRelax::run()
//  only useful for the naive wilson gauge action. 
//  checkerboard over nodes is not done.

//! Run the overrelaxation algorithm for the Wilson (plaquette) gauge action.
//------------------------------------------------------------------
void AlgOverRelax::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  LRG.SetInterval(1,-1);
int x[4], y[4];
  Matrix mStaple;
  // Run the gauge heat bath
  //----------------------------------------------------------------
  for(int i=0; i< alg_overrelax_arg->num_iter; i++){

     // Over Relaxation
     // Checkerboard everything, do even or odd sites
     // Scan over local subvolume doing all even(odd) links
     // index for "x" is consistant with GsiteOffset: "x,y,z,t" order

	// The y[4] are the local coordinates for the 2^4 cube. 
	// We traverse the points on this cube, and do the corresponding
	// points on every other hypercube before moving on.

  
  for( int checker=0; checker < 2; checker++)  {
    for( y[3] = 0; y[3] < 2; y[3]++) 		
    for( y[2] = 0; y[2] < 2; y[2]++)
    for( y[1] = 0; y[1] < 2; y[1]++) 		
    for( y[0] = 0; y[0] < 2; y[0]++) 		
      if( (y[0]+y[1]+y[2]+y[3])%2 == checker)  {

        for( x[3] = y[3]; x[3] < GJP.TnodeSites(); x[3]+=2) 
        for( x[2] = y[2]; x[2] < GJP.ZnodeSites(); x[2]+=2) 
        for( x[1] = y[1]; x[1] < GJP.YnodeSites(); x[1]+=2) 
        for( x[0] = y[0]; x[0] < GJP.XnodeSites(); x[0]+=2) 
          for( int mu=0; mu<4; mu++ ) {

            LRG.AssignGenerator(x);

#ifdef OVERRELAX_TIMING
	    int clock1 = clock();
#endif

       	    // Build the staple:
            lat.Staple( mStaple, x, mu );

            Matrix * pmLink = lat.GaugeField();
            pmLink += mu + lat.GsiteOffset(x);
	    //   printf("wqw %le\n", (*pmLink).ErrorSU3());
	        UpdateLink(pmLink, mStaple);

#ifdef OVERRELAX_TIMING
    	    printf("OVER: one link update time %d\n", clock()-clock1);
#endif	
          }			// (end mu loop)
        }           		// (end of site loop)
      }				// (end checkerboard loop)
    

    // Increment gauge field counter
    //----------------------------------------------------------------
    lat.GupdCntInc(1);

    // If GJP.Snodes() !=1  the gauge field is spread out
    // accross s-slices of processors. It must be identical
    // on each slice. Check to make sure and exit if it
    // is not identical. A case where this is relevant
    // is the DWF spread-out case.
    //----------------------------------------------------------------
    //     lat.GsoCheck();
    //    cout<<".";
  }

  // Reunitarize
  //----------------------------------------------------------------
  lat.Reunitarize();
  //VRB.FuncEnd(fname, cname);   
}
//Matrix mirror(Matrix  & Link, const  Matrix & Staple)
//{
  //U=Dagger(Staple*Link*Staple);

//}
void AlgOverRelax::UpdateLink(Matrix  * pmLink, const  Matrix & mStaple)
{
  
  Float  old_action,  new_action;
  Float  accept_probability;

  Matrix Trial, Snew,Sold, Temp1, Temp2;
  // remember - you cannot use multi-hit here, it breaks detailed balance
  //generate G0xGxG0
  Trial.Dagger(* pmLink);
  Temp1.DotMEqual(mStaple, Trial);
  Temp2.DotMEqual(Temp1,mStaple);

  //  printf("qqq %le %le %le\n",(* pmLink).ErrorSU3(),mStaple.ErrorSU3(),Trial.ErrorSU3());
  Sold.DotMEqual(*pmLink,mStaple);
  Snew.DotMEqual(Trial,mStaple);

  
  // Now begin the update process. 
  
  old_action = fGamma * (9. - .5 * Sold.ReTr() );
  
  new_action = fGamma * (9. - .5 * Snew.ReTr());
  
  accept_probability = exp( old_action - new_action );
  
  
  if(fabs(LRG.Urand()) < accept_probability ) {    
    for(int i=0;i<18;i++){ 
           *((Float*)pmLink +  i ) =Trial.elem(i);

   }
        
  }
  
}



CPS_END_NAMESPACE
