#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Definitions of the AlgGheatBath class methods.
  
  $Id: alg_ghb.C,v 1.15 2013-04-05 17:51:13 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/comsrc/alg_ghb.C,v 1.15 2013-04-05 17:51:13 chulwoo Exp $
//  $Id: alg_ghb.C,v 1.15 2013-04-05 17:51:13 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_ghb.C,v $
//  $Revision: 1.15 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/comsrc/alg_ghb.C,v $
//  $State: Exp $
//

CPS_END_NAMESPACE
#include <stdlib.h>	// exit()
#include <util/qcdio.h>
#include <math.h>
#include <time.h>
#include <alg/alg_ghb.h>
#include <alg/common_arg.h>
#include <alg/ghb_arg.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/random.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
////#include <mem/p2v.h>
CPS_START_NAMESPACE

//Uncomment next line to switch on timing
//#define GHB_TIMING

#if TARGET == QCDSP
static int dram_grand_seed[6] = {
  0xbeefcafe,
  0xdefaced, 
  0xfacade, 
  0xbabeface,
  1,		// carry bit
  1013904243    // odd const for LCG 
};
#endif
  

//------------------------------------------------------------------
/*!
  \param latt The lattice on which to perform the heatbath
  \param c_arg The common argument structure for all algorithms.
  \param arg The parameters specific to this algorithm.
 */
//------------------------------------------------------------------
AlgGheatBath::AlgGheatBath(Lattice& latt, 
	     CommonArg *c_arg,
	     GhbArg *arg) : 
	     Alg(latt, c_arg) 
{
  cname = "AlgGheatBath";
  char *fname = "AlgGheatBath(L&,CommonArg*,GhbArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_ghb_arg = arg;

  // Setup the seeds for the asm random number generator
  //----------------------------------------------------------------
  {
/*	This won't work because the Assembly code for QCDSP still contains
		references to the old system, and has not yet been changed
    LRG.SetInterval(1,0);
    LRG.AssignGenerator(0);
    dram_grand_seed[0] = (int) (LRG.Urand()*1000000);
    dram_grand_seed[1] = (int) (LRG.Urand()*1000000);
    dram_grand_seed[2] = (int) (LRG.Urand()*1000000);
    dram_grand_seed[3] = (int) (LRG.Urand()*1000000);
*/
  }
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgGheatBath::~AlgGheatBath() {
  //char *fname = "~AlgGheatBath()";
  //VRB.Func(cname,fname);
}

void m_conjugate( Float* );
void cmhb_kernel( Float*, Float* );
void metropolis_kernel( Float*, Float* );

void AlgGheatBath::relocate(){
#if TARGET == QCDSP
  /*{
    // relocate assembly code into CRAM
//    p2vGhb();

    // Special initialization actions required for relocatable 
    // assembly.  Must be preformed _after_ asembly code has been
    // relocated.

    extern	int*	core_iscratch;
    extern	Float*	core_fscratch;
    extern	int*	el_seed_p;

    *( core_fscratch + 57 ) = (2.*GJP.Beta())/3.;

    // Set the "fast" seeds in CRAM to be the same as those
    // in the DRAM area assigned for this job. 

    *(el_seed_p+0) = dram_grand_seed[0];
    *(el_seed_p+1) = dram_grand_seed[1];
    *(el_seed_p+2) = dram_grand_seed[2];
    *(el_seed_p+3) = dram_grand_seed[3];
    *(el_seed_p+4) = dram_grand_seed[4];
    *(el_seed_p+5) = dram_grand_seed[5];
  }*/
#endif
}

void AlgGheatBath::preserve_seed(){
#if TARGET == QCDSP
/*    {
      extern    int*    el_seed_p;
      // put the "fast" seeds back into dram 
      dram_grand_seed[0] = *(el_seed_p+0);
      dram_grand_seed[1] = *(el_seed_p+1);
      dram_grand_seed[2] = *(el_seed_p+2);
      dram_grand_seed[3] = *(el_seed_p+3);
      dram_grand_seed[4] = *(el_seed_p+4);
      dram_grand_seed[5] = *(el_seed_p+5);
    }*/
#endif
}

//------------------------------------------------------------------
// AlgGheatBath::run()
//  only useful for the naive wilson gauge action. 
//  checkerboard over nodes is not done.

//! Run the heatbath algorithm for the Wilson (plaquette) gauge action.
//------------------------------------------------------------------
void AlgGheatBath::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();
  if (lat.Gclass() != G_CLASS_WILSON){
    ERR.General(cname,fname," Only correct for Wilson gauge action\n");
  }
  int i; 

  // relocate the heatbath kernal and set up rand seeds
  //--------------------------------------------------
  relocate();

  LRG.SetInterval(1,-1);

  // Run the gauge heat bath
  //----------------------------------------------------------------
  for(i=0; i< alg_ghb_arg->num_iter; i++){

     // Heat bath
     // Checkerboard everything, do even or odd sites
     // Scan over local subvolume doing all even(odd) links
     // index for "x" is consistant with GsiteOffset: "x,y,z,t" order

	// The y[4] are the local coordinates for the 2^4 cube. 
	// We traverse the points on this cube, and do the corresponding
	// points on every other hypercube before moving on.

  int x[4], y[4];
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

#ifdef GHB_TIMING
	    int clock1 = clock();
#endif

       	    // Build the staple:
            Matrix mStaple;
            lat.Staple( mStaple, x, mu );

            Matrix * pmLink = lat.GaugeField();
            pmLink += mu + lat.GsiteOffset(x);

            UpdateLink(pmLink, mStaple);

#ifdef GHB_TIMING
    	    printf("one link update time %d\n", clock()-clock1);
#endif	
          }			// (end mu loop)
        }           		// (end of site loop)
      }				// (end checkerboard loop)
    preserve_seed();	

    // Increment gauge field counter
    //----------------------------------------------------------------
    lat.GupdCntInc(1);

    // If GJP.Snodes() !=1  the gauge field is spread out
    // accross s-slices of processors. It must be identical
    // on each slice. Check to make sure and exit if it
    // is not identical. A case where this is relevant
    // is the DWF spread-out case.
    //----------------------------------------------------------------
    lat.GsoCheck();
  }

  // Reunitarize
  //----------------------------------------------------------------
  lat.Reunitarize();
  //VRB.FuncEnd(fname, cname);   
}

//------------------------------------------------------------------
// AlgGheatBath::NoCheckerBoardRun()
//  only useful for the naive wilson gauge action. 
//  checkerboard over nodes is not done.
//! Run the heatbath algorithm for the Wilson (plaquette) gauge action.
/*!
  This is identical to AlgGheatBath::run apart from a call to
  Lattice::ClearAllBufferedLink at the end.
*/
  
//------------------------------------------------------------------
void AlgGheatBath::NoCheckerBoardRun()
{
  char *fname = "NoCheckerBoardRun()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();
  if (lat.Gclass() != G_CLASS_WILSON){
    ERR.General(cname,fname," Only correct for Wilson gauge action\n");
  }

  relocate(); 
  LRG.SetInterval(1,-1);

  // Run the gauge heat bath
  //----------------------------------------------------------------
  for(int i=0; i< alg_ghb_arg->num_iter; i++){

    // Heat bath

    // Checkerboard everything, do even or odd sites
    for( int checker=0; checker <=1; checker++ ) {

      // Scan over local subvolume doing all even(odd) links
      // index for "x" is consistant with GsiteOffset: "x,y,z,t" order
    int x[4], y[4];
    for( y[3] = 0; y[3] < 2; y[3]++)
    for( y[2] = 0; y[2] < 2; y[2]++)
    for( y[1] = 0; y[1] < 2; y[1]++)
    for( y[0] = 0; y[0] < 2; y[0]++)
      if( (y[0]+y[1]+y[2]+y[3])%2 == checker)  {

        for( x[3] = y[3]; x[3] < GJP.TnodeSites(); x[3]+=2)
        for( x[2] = y[2]; x[2] < GJP.ZnodeSites(); x[2]+=2)
        for( x[1] = y[1]; x[1] < GJP.YnodeSites(); x[1]+=2)
        for( x[0] = y[0]; x[0] < GJP.XnodeSites(); x[0]+=2)
          for( int mu=0; mu<4; mu++ )   {

#ifdef GHB_TIMING
	  int clock1 = clock();
#endif
          // Build the staple:
          Matrix mStaple;
          lat.AllStaple( mStaple, x, mu );
          lat.ClearBufferedLink(x,mu);

          Matrix * pmLink = lat.GaugeField();
          pmLink += mu + lat.GsiteOffset(x);
          UpdateLink(pmLink, mStaple);

#ifdef GHB_TIMING
	  printf("one link update time %d\n", clock()-clock1);
#endif
        }			// (end mu loop)
      }  			// (end site loop)
    }  				// (end checkerboard loop)
  
    preserve_seed();

    // Increment gauge field counter
    //----------------------------------------------------------------
    lat.GupdCntInc(1);

    // If GJP.Snodes() !=1  the gauge field is spread out
    // accross s-slices of processors. It must be identical
    // on each slice. Check to make sure and exit if it
    // is not identical. A case where this is relevant
    // is the DWF spread-out case.
    //----------------------------------------------------------------
    lat.GsoCheck();
  }
  // Reunitarize
  lat.Reunitarize();
  lat.ClearAllBufferedLink();
}



//------------------------------------------------------------------
//AlgGheatBath::NodeCheckerBoardRun()
// Do checkerboard over the nodes when at least one of the direction
// perpendicular to the updated link has node_sites of 2. 
// This is necessary for the correctness of the heatbath algorithm, when
// the rectangle action is involved.
//! Run the heatbath algorithm 
//------------------------------------------------------------------
void AlgGheatBath::NodeCheckerBoardRun()
{
  char *fname = "NodeCheckerBoardRun()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  int do_checker[4]= {0, 0, 0, 0};

  for(int u=0;u<4;u++)
    if(GJP.NodeSites(u)==2) 
      for(int v=0;v<4;v++)
        if (v!= u)
          do_checker[v]= 1;

  relocate();
  LRG.SetInterval(1,-1);

  int node_checker = (GJP.XnodeCoor()+GJP.YnodeCoor()
                      +GJP.ZnodeCoor()+GJP.TnodeCoor())&1;
  
  //VRB.Flow(cname, fname, 
   //        "do_checker0=%d,do_checker1=%d,do_checker2=%d,do_checker3=%d\n", 
    //        do_checker[0],do_checker[1],do_checker[2],do_checker[3]);

  // Run the gauge heat bath
  //----------------------------------------------------------------
  for(int i=0; i< alg_ghb_arg->num_iter; i++){

    // Heat bath
    // Checkerboard everything, do even or odd sites
    for( int checker=0; checker <=1; checker++ ) {

      // Scan over local subvolume doing all even(odd) links
      // index for "x" is consistant with GsiteOffset: "x,y,z,t" order
      int x[4], y[4];
      for( y[3] = 0; y[3] < 2; y[3]++)
      for( y[2] = 0; y[2] < 2; y[2]++)
      for( y[1] = 0; y[1] < 2; y[1]++)
      for( y[0] = 0; y[0] < 2; y[0]++)
      if( (y[0]+y[1]+y[2]+y[3])%2 == checker)  {

 //     printf("y = %d %d %d %d\n",y[0],y[1],y[2],y[3]);
      for( x[3] = y[3]; x[3] < GJP.TnodeSites(); x[3]+=2)
      for( x[2] = y[2]; x[2] < GJP.ZnodeSites(); x[2]+=2)
      for( x[1] = y[1]; x[1] < GJP.YnodeSites(); x[1]+=2)
      for( x[0] = y[0]; x[0] < GJP.XnodeSites(); x[0]+=2)



        for( int mu=0; mu<4; mu++ )   {
//      printf("x = %d %d %d %d\n",x[0],x[1],x[2],x[3]);fflush(stdout);


#ifdef GHB_TIMING  
	  int clock1 = clock();
#endif
        for( int node_ck =0; node_ck<2; node_ck++){

         // Build the staple:

         //if don't do checkerboard over nodes, then only get the staple
         //once. if we need to do checkerboard over the nodes, we need to 
         //get the staple twice, for every node must be doing exactly same 
         //transmission of links. 

         //VRB.Flow(cname, fname, "x=%d y=%d z=%d t=%d mu=%d\n", 
           //       x[0], x[1], x[2], x[3], mu);
         Matrix mStaple;
         if((do_checker[mu]==0 && node_ck==0) || do_checker[mu] == 1) {
          //VRB.Flow(cname, fname, "mu=%d do_checker=%d node_ck=%d\n", 
            //       mu, do_checker[mu], node_ck);
          //lat.ClearBufferedLink(x,mu);
          lat.AllStaple( mStaple, x, mu );
          lat.ClearBufferedLink(x,mu);
         }
         
         //if don't do the checkerboard over nodes, then all nodes 
         //update the link together once. if we need to do checkerboard over 
         //nodes do the update of the link only when node_checker==node_ck 
         //
         if((do_checker[mu]==0 && node_ck==0) || 
            (do_checker[mu]==1 && node_checker == node_ck)){  

            Matrix * pmLink = lat.GaugeField();
            pmLink += mu + lat.GsiteOffset(x);
            UpdateLink(pmLink, mStaple);


         }      //end inner if
        }       //end of node_ck loop


#ifdef GHB_TIMING
	printf("one link update time %d\n", clock()-clock1);
#endif
    
       }			// (end mu loop )
      }  			// (end site loop)
    }  				// (end checkerboard loop)

    preserve_seed();

    // Increment gauge field counter
    //----------------------------------------------------------------
    lat.GupdCntInc(1);

    // If GJP.Snodes() !=1  the gauge field is spread out
    // accross s-slices of processors. It must be identical
    // on each slice. Check to make sure and exit if it
    // is not identical. A case where this is relevant
    // is the DWF spread-out case.
    //----------------------------------------------------------------
    lat.GsoCheck();
  }
  // Reunitarize
  lat.Reunitarize();
  lat.ClearAllBufferedLink();
  //VRB.FuncEnd(cname, fname);
}


void
AlgGheatBath::UpdateLink(Matrix * pmLink, const Matrix & mStaple){

          Float* pfStapleCMHBorder;

/*  # ifdef _TARTAN
	  extern Float* fast_sigma; 
	  pfStapleCMHBorder =  fast_sigma; 
# else
*/
	  Float staple_buf[18];
	  pfStapleCMHBorder = (Float*)staple_buf;

// # endif

          *(pfStapleCMHBorder +  0) = *( (Float*)&mStaple +  0 );
          *(pfStapleCMHBorder +  9) = *( (Float*)&mStaple +  1 );
          *(pfStapleCMHBorder +  1) = *( (Float*)&mStaple +  2 );
          *(pfStapleCMHBorder + 10) = *( (Float*)&mStaple +  3 );
          *(pfStapleCMHBorder +  2) = *( (Float*)&mStaple +  4 );
          *(pfStapleCMHBorder + 11) = *( (Float*)&mStaple +  5 );
          *(pfStapleCMHBorder +  3) = *( (Float*)&mStaple +  6 );
          *(pfStapleCMHBorder + 12) = *( (Float*)&mStaple +  7 );
          *(pfStapleCMHBorder +  4) = *( (Float*)&mStaple +  8 );
          *(pfStapleCMHBorder + 13) = *( (Float*)&mStaple +  9 );
          *(pfStapleCMHBorder +  5) = *( (Float*)&mStaple + 10 );
          *(pfStapleCMHBorder + 14) = *( (Float*)&mStaple + 11 );
          *(pfStapleCMHBorder +  6) = *( (Float*)&mStaple + 12 );
          *(pfStapleCMHBorder + 15) = *( (Float*)&mStaple + 13 );
          *(pfStapleCMHBorder +  7) = *( (Float*)&mStaple + 14 );
          *(pfStapleCMHBorder + 16) = *( (Float*)&mStaple + 15 );
          *(pfStapleCMHBorder +  8) = *( (Float*)&mStaple + 16 );
          *(pfStapleCMHBorder + 17) = *( (Float*)&mStaple + 17 );

          Float* pfLinkCMHBorder;

/*  # ifdef _TARTAN
	  extern Float* fast_link; 
	  pfLinkCMHBorder =  fast_link; 
# else
*/
	  Float link_buf[18];
	  pfLinkCMHBorder = (Float*)link_buf;

//	  printf("link_buf = %p  pfLinkCMHBorder = %p pmLink = %p\n",link_buf,pfLinkCMHBorder,pmLink);
// # endif

          *(pfLinkCMHBorder +  0) = *((Float*)pmLink +  0 );
          *(pfLinkCMHBorder +  9) = *((Float*)pmLink +  1 );
          *(pfLinkCMHBorder +  1) = *((Float*)pmLink +  2 );
          *(pfLinkCMHBorder + 10) = *((Float*)pmLink +  3 );
          *(pfLinkCMHBorder +  2) = *((Float*)pmLink +  4 );
          *(pfLinkCMHBorder + 11) = *((Float*)pmLink +  5 );
          *(pfLinkCMHBorder +  3) = *((Float*)pmLink +  6 );
          *(pfLinkCMHBorder + 12) = *((Float*)pmLink +  7 );
          *(pfLinkCMHBorder +  4) = *((Float*)pmLink +  8 );
          *(pfLinkCMHBorder + 13) = *((Float*)pmLink +  9 );
          *(pfLinkCMHBorder +  5) = *((Float*)pmLink + 10 );
          *(pfLinkCMHBorder + 14) = *((Float*)pmLink + 11 );
          *(pfLinkCMHBorder +  6) = *((Float*)pmLink + 12 );
          *(pfLinkCMHBorder + 15) = *((Float*)pmLink + 13 );
          *(pfLinkCMHBorder +  7) = *((Float*)pmLink + 14 );
          *(pfLinkCMHBorder + 16) = *((Float*)pmLink + 15 );
          *(pfLinkCMHBorder +  8) = *((Float*)pmLink + 16 );
          *(pfLinkCMHBorder + 17) = *((Float*)pmLink + 17 );

          // Call the heat bath
          {
	    metropolis_kernel( pfStapleCMHBorder, pfLinkCMHBorder );
//            cmhb_kernel( pfStapleCMHBorder, pfLinkCMHBorder );
          }
          
          // Copy the link back into the lattice, and
          // arrange the internal storage order back to
          // the "canonical" order.

          *((Float*)pmLink +  0 ) = *(pfLinkCMHBorder +  0);
          *((Float*)pmLink +  1 ) = *(pfLinkCMHBorder +  9);
          *((Float*)pmLink +  2 ) = *(pfLinkCMHBorder +  1);
          *((Float*)pmLink +  3 ) = *(pfLinkCMHBorder + 10);
          *((Float*)pmLink +  4 ) = *(pfLinkCMHBorder +  2);
          *((Float*)pmLink +  5 ) = *(pfLinkCMHBorder + 11);
          *((Float*)pmLink +  6 ) = *(pfLinkCMHBorder +  3);
          *((Float*)pmLink +  7 ) = *(pfLinkCMHBorder + 12);
          *((Float*)pmLink +  8 ) = *(pfLinkCMHBorder +  4);
          *((Float*)pmLink +  9 ) = *(pfLinkCMHBorder + 13);
          *((Float*)pmLink + 10 ) = *(pfLinkCMHBorder +  5);
          *((Float*)pmLink + 11 ) = *(pfLinkCMHBorder + 14);
          *((Float*)pmLink + 12 ) = *(pfLinkCMHBorder +  6);
          *((Float*)pmLink + 13 ) = *(pfLinkCMHBorder + 15);
          *((Float*)pmLink + 14 ) = *(pfLinkCMHBorder +  7);
          *((Float*)pmLink + 15 ) = *(pfLinkCMHBorder + 16);
          *((Float*)pmLink + 16 ) = *(pfLinkCMHBorder +  8);
          *((Float*)pmLink + 17 ) = *(pfLinkCMHBorder + 17);

}














CPS_END_NAMESPACE
