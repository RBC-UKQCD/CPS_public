#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FwilsonTypes class.

  $Id: f_wilson_t.C,v 1.7 2005/10/04 05:35:47 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005/10/04 05:35:47 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/f_wilson_types/comsrc/f_wilson_t.C,v 1.7 2005/10/04 05:35:47 chulwoo Exp $
//  $Id: f_wilson_t.C,v 1.7 2005/10/04 05:35:47 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: f_wilson_t.C,v $
//  $Revision: 1.7 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/f_wilson_types/comsrc/f_wilson_t.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_wilson_types.C
//
// FwilsonTypes is derived from Lattice and is relevant to
// all fermion classes with Wilson type fermions 
// (e.g Fwilson, Fclover, Fdwf, ...). These classes are derived
// from FwilsonTypes
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/enum.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/sproj_tr.h>
CPS_START_NAMESPACE

extern void gamma_5(IFloat *v_out, IFloat *v_in, int num_sites) ;

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
FwilsonTypes::FwilsonTypes()
{
  cname = "FwilsonTypes";
  char *fname = "FwilsonTypes()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Fill in the array of sproj_tr functions
  //----------------------------------------------------------------
  sproj_tr[SPROJ_XM] = sprojTrXm;
  sproj_tr[SPROJ_YM] = sprojTrYm;
  sproj_tr[SPROJ_ZM] = sprojTrZm;
  sproj_tr[SPROJ_TM] = sprojTrTm;
  sproj_tr[SPROJ_XP] = sprojTrXp;
  sproj_tr[SPROJ_YP] = sprojTrYp;
  sproj_tr[SPROJ_ZP] = sprojTrZp;
  sproj_tr[SPROJ_TP] = sprojTrTp;

  //----------------------------------------------------------------
  // Fill in the array of Sigma proj_tr functions
  //----------------------------------------------------------------

  Sigmaproj_tr[SIGMAPROJ_XY] = SigmaprojTrXY;
  Sigmaproj_tr[SIGMAPROJ_XZ] = SigmaprojTrXZ;
  Sigmaproj_tr[SIGMAPROJ_XT] = SigmaprojTrXT;
  Sigmaproj_tr[SIGMAPROJ_YX] = SigmaprojTrYX;
  Sigmaproj_tr[SIGMAPROJ_YZ] = SigmaprojTrYZ;
  Sigmaproj_tr[SIGMAPROJ_YT] = SigmaprojTrYT;
  Sigmaproj_tr[SIGMAPROJ_ZX] = SigmaprojTrZX;
  Sigmaproj_tr[SIGMAPROJ_ZY] = SigmaprojTrZY;
  Sigmaproj_tr[SIGMAPROJ_ZT] = SigmaprojTrZT;
  Sigmaproj_tr[SIGMAPROJ_TX] = SigmaprojTrTX;
  Sigmaproj_tr[SIGMAPROJ_TY] = SigmaprojTrTY;
  Sigmaproj_tr[SIGMAPROJ_TZ] = SigmaprojTrTZ;

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
FwilsonTypes::~FwilsonTypes()
{
  char *fname = "~FwilsonTypes()";
  VRB.Func(cname,fname);

}


//------------------------------------------------------------------
    /*!
      The vectors should have 4x3 complex components per site.
      gamma_5 is in the chiral basis:
    
               [ 1  0  0  0]
     gamma_5 = [ 0  1  0  0]
               [ 0  0 -1  0]
               [ 0  0  0 -1]
          
      \param v_out The product vector.
      \param v_in The vector to be multiplied,
      \param num_sites The number of lattice sites at which to do the
      multiplication.
    */
//------------------------------------------------------------------
void FwilsonTypes::Gamma5(Vector *v_out, Vector *v_in, int num_sites)
{
  char *fname = "Gamma5(V*,V*,i)" ;
 
  if (Colors() != 3)
    ERR.General(cname,fname,"Wrong nbr of colors.") ;
 
  if (SpinComponents() != 4)
    ERR.General(cname,fname,"Wrong nbr of spin comp.") ;
 
  gamma_5((IFloat *)v_out, (IFloat *)v_in, num_sites) ;
}

//------------------------------------------------------------------
// int ExactFlavors() : 
// Returns the number of exact flavors of the matrix that
// is inverted during a molecular dynamics evolution.
//------------------------------------------------------------------
int FwilsonTypes::ExactFlavors() const
{
  return 2;
}

int FwilsonTypes::FsiteSize() const
{
  return 2 * Colors() * SpinComponents();  
  // re/im * colors * spin_components
}

//------------------------------------------------------------------
// Float FhamiltonNode(Vector *phi, Vector *chi):
// The fermion Hamiltonian of the node sublattice.
// chi must be the solution of Cg with source phi.	       
//------------------------------------------------------------------
Float FwilsonTypes::FhamiltonNode( Vector *phi,  Vector *chi) {
  char *fname = "FhamiltonNode(V*,V*)";
  VRB.Func(cname,fname);

  if (phi == 0)
    ERR.Pointer(cname,fname,"phi") ;

  if (chi == 0)
    ERR.Pointer(cname,fname,"chi") ;

  return phi->ReDotProductNode(chi,((GJP.VolNodeSites()*FsiteSize())>>1)*(2-FchkbEvl())) ;
}

//------------------------------------------------------------------
// int FsiteOffsetChkb(const int *x):
// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is not the canonical one but it is particular
// to the fermion type. x[i] is the 
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
//------------------------------------------------------------------
int FwilsonTypes::FsiteOffsetChkb(const int *x) const {

  return ( x[0] - x[0]%2 + GJP.XnodeSites() *
           ( x[1] + GJP.YnodeSites() *
             (x[2] + GJP.ZnodeSites() * x[3])) +
           ((x[0]+x[1]+x[2]+x[3]+1)%2)*GJP.VolNodeSites()
         )>>1 ;
}

int FwilsonTypes::FsiteOffset(const int *x) const {

  return ( x[0] + GJP.XnodeSites() *
	   ( x[1] + GJP.YnodeSites() *
	     ( x[2] + GJP.ZnodeSites() * x[3])));
}

int FwilsonTypes::SpinComponents() const
{
  return 4;
}

//!< Dummy routine for Wilson fermions
void FwilsonTypes::BforceVector(Vector *in, CgArg *cg_arg) {

}

CPS_END_NAMESPACE
