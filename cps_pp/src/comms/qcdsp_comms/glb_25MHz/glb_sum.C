#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_25MHz/glb_sum.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: glb_sum.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2001/08/16 12:53:58  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.4  2001/08/16 10:49:45  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:11:45  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:01  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: glb_sum.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_25MHz/glb_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// glb_sum.C
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
#include <sysfunc.h>                // CoorT(), ...
#include <stdlib.h>                 // void *malloc(unsigned size)
CPS_START_NAMESPACE

static int is_initialized = 0;


void glb_sum(Float *float_p)
{
  int masterCoord[4];
  int thisCoord[4];
  int size[4];
  char *cname = " ";
  char *fname = "glb_sum";

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
