#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:45 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_hmd/alg_hmd.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Id: alg_hmd.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.3  2002/12/04 17:16:27  zs
//  Merged the new 2^4 RNG into the code.
//  This new RNG is implemented in the LatRanGen class.
//  The following algorithm and utility classes are affected:
//
//  AlgEig                  Fdwf
//  AlgGheatBath            Fstag
//  AlgHmd                  GlobalJobParameter
//  AlgNoise                Lattice
//  AlgPbp                  Matrix
//  AlgThreept              RandomGenerator
//                          Vector
//
//  Revision 1.2  2001/06/19 18:11:27  anj
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
//  Revision 1.2  2001/05/25 06:15:59  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: alg_hmd.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_hmd/alg_hmd.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_hmd.C
//
// AlgHmd is derived from Alg and is relevant to the Hybrid  
// Molecular Dynamics algorithms.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <string.h>
#include<alg/alg_hmd.h>
#include<alg/common_arg.h>
#include<alg/hmd_arg.h>
#include<alg/cg_arg.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
AlgHmd::AlgHmd(Lattice& latt, 
		     CommonArg *c_arg,
		     HmdArg *arg) : 
		     Alg(latt, c_arg) 
{
  cname = "AlgHmd";
  char *fname = "AlgHmd(L&,CommonArg*,HmdArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  hmd_arg = arg;
  
  // Print out arguments
  //----------------------------------------------------------------
  //??? Must be done.

  // Calculate the gauge field size.
  //----------------------------------------------------------------
  g_size = GJP.VolNodeSites() * latt.GsiteSize();

  // Allocate memory for the conjugate momenta.
  //----------------------------------------------------------------
  mom = (Matrix *) smalloc(g_size * sizeof(Float));
  if(mom == 0)
    ERR.Pointer(cname,fname, "mom");
  VRB.Smalloc(cname,fname,
	      "mom",mom, g_size * sizeof(Float));

  // Initialize the Number of Checkerboards for fermion field (Ncb)
  //----------------------------------------------------------------
  if(latt.FchkbEvl() == 1)
    Ncb = 1;                    // Half Checkerboard
  else if(latt.FchkbEvl() == 0)
    Ncb = 2;                    // Full Checkerboard
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgHmd::~AlgHmd() {
  char *fname = "~AlgHmd()" ;
  VRB.Func(cname,fname);


  // Free memory for the conjugate momenta.
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "mom",mom);
  sfree(mom);

}






CPS_END_NAMESPACE
