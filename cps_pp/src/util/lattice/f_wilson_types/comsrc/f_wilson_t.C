#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/comsrc/f_wilson_t.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: f_wilson_t.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:35  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:25  anj
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
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: f_wilson_t.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/comsrc/f_wilson_t.C,v $
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
#include<util/lattice.h>
#include<util/enum.h>
#include<util/verbose.h>
#include<util/gjp.h>
#include<util/error.h>
#include<util/sproj_tr.h>
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
// Gamma5(Vector *v_out, Vector *v_in, int num_sites):
//
// v_out = Gamma5 * v_in. Gamme5 is in the chiral basis
//
//          [ 1  0  0  0]
// Gamma5 = [ 0  1  0  0]
//          [ 0  0 -1  0]
//          [ 0  0  0 -1]
//
// num_sites is the number of sites. It is assumed
// that each site has 24 components.
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


CPS_END_NAMESPACE
