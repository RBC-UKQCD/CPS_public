#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Definitions of the AlgHmd constructor and destructor.

  $Id: alg_hmd.C,v 1.10 2005/02/18 19:56:00 mclark Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mclark $
//  $Date: 2005/02/18 19:56:00 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_hmd/alg_hmd.C,v 1.10 2005/02/18 19:56:00 mclark Exp $
//  $Id: alg_hmd.C,v 1.10 2005/02/18 19:56:00 mclark Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: alg_hmd.C,v $
//  $Revision: 1.10 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_hmd/alg_hmd.C,v $
//  $State: Exp $
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
#include <alg/alg_hmd.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>
#include <alg/cg_arg.h>
#include <util/lattice.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
/*!
  Allocates memory for all the fields used and copies parameters
  to internal locations.
  \param latt The lattice on which run the HMD algorithm.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
 */
//------------------------------------------------------------------
AlgHmd::AlgHmd(Lattice& latt, CommonArg *c_arg, HmdArg *arg) : 
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
  if(mom == 0) ERR.Pointer(cname,fname, "mom");
  VRB.Smalloc(cname,fname,"mom",mom, g_size * sizeof(Float));

  // Initialize the Number of Checkerboards for fermion field (Ncb)
  //----------------------------------------------------------------
  if(latt.FchkbEvl() == 1)
    Ncb = 1;                    // Half Checkerboard
  else if(latt.FchkbEvl() == 0)
    Ncb = 2;                    // Full Checkerboard

}


//----------------------------------------------------------------
/*!  Frees allocated memory.*/
//----------------------------------------------------------------

AlgHmd::~AlgHmd() {
  char *fname = "~AlgHmd()" ;
  VRB.Func(cname,fname);

  // Free memory for the conjugate momenta.
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "mom",mom);
  sfree(mom);

}



    







CPS_END_NAMESPACE
