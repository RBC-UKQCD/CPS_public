#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:03 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_50MHz_nos/glb_sum_five.C,v 1.3 2004-06-04 21:14:03 chulwoo Exp $
//  $Id: glb_sum_five.C,v 1.3 2004-06-04 21:14:03 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_sum_five.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_50MHz_nos/glb_sum_five.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// glb_sum_five.C
//
// Wrapper for Ping's optimized global sum. This will
// eventually be part of the qos.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<comms/glb.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/gjp.h>
#include<global_sum.h>
#include<comms/glb_sum_init.h>
#include <comms/sysfunc.h>                // CoorT(), ...
#include <stdlib.h>                 // void *malloc(unsigned size)
CPS_START_NAMESPACE

static int is_initialized = 0;


void glb_sum_five(Float *float_p)
{
  int masterCoord[4];
  int thisCoord[4];
  int size[4];
  char *cname = " ";
  char *fname = "glb_sum_five";

  //----------------------------------------------------------------
  // If this is the first call initialize the global sum parameters
  //----------------------------------------------------------------
  if(is_initialized == 0){

    //--------------------------------------------------------------
    // remap the physics coordinates into machine coordinates.
    //--------------------------------------------------------------
    {
      int physSize[4] =  {SizeT(), SizeX(), SizeY(), SizeZ()};
      int masterPhysCoord[4] = {0, 0, 0, 0};
      int thisPhysCoord[4] = { CoorT(), CoorX(), CoorY(), CoorZ()};
      
      for (int physDim = 0; physDim < 4; ++physDim) {
	int machDim = SCURemap((SCUDir)(2 * physDim)) / 2; 
	size[machDim]        = physSize[physDim];    
	masterCoord[machDim] = masterPhysCoord[physDim];
	thisCoord[machDim]   = thisPhysCoord[physDim];
      }
    }

    //--------------------------------------------------------------
    // Initialize global sum parameters
    //--------------------------------------------------------------
    glb_sum_init(masterCoord, thisCoord, size, (unsigned int)malloc(6), 0);
    VRB.Flow(cname, fname,"Initialized global sum\n");
    is_initialized = 1;
  }


  //----------------------------------------------------------------
  // Do the global sum over the full machine 
  //----------------------------------------------------------------
  Float sum;
  IFloat sum_f;
  int partitions;
  sum = *float_p;
  sum_f = sum;
  int exp = global_sum(&sum_f);

  //----------------------------------------------------------------
  // Exit if overflow
  //----------------------------------------------------------------
  if (exp > 127) {
    ERR.General(cname, fname, "Global sum overflow.\n");
  }

  //----------------------------------------------------------------
  // Divide the sum by the number of partitions.
  //----------------------------------------------------------------
  sum = sum_f;
  partitions = SizeX() / GJP.Xnodes();
  partitions = partitions * ( SizeY() / GJP.Ynodes() );
  partitions = partitions * ( SizeZ() / GJP.Znodes() );
  partitions = partitions * ( SizeT() / GJP.Tnodes() );
  sum = sum / Float(partitions);
  *float_p = sum;

}


CPS_END_NAMESPACE
