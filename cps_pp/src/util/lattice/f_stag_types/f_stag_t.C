#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Definition of FstagTypes methods.

  $Id: f_stag_t.C,v 1.4 2004-01-13 20:39:50 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:50 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_stag_types/f_stag_t.C,v 1.4 2004-01-13 20:39:50 chulwoo Exp $
//  $Id: f_stag_t.C,v 1.4 2004-01-13 20:39:50 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.3.2.1  2003/11/05 16:32:16  mike
//  Attempting to create new working branch!
//
//  Revision 1.2  2003/07/24 16:53:54  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.2  2001/06/19 18:13:23  anj
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
//  $RCSfile: f_stag_t.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_stag_types/f_stag_t.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_stag_types.C
//
// FstagTypes is derived from Lattice and is relevant to
// all fermion classes with Staggered type fermions 
// These classes are derived from FstagTypes
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

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
FstagTypes::FstagTypes()
{
  cname = "FstagTypes";
  char *fname = "FstagTypes()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
FstagTypes::~FstagTypes()
{
  char *fname = "~FstagTypes()";
  VRB.Func(cname,fname);
}

CPS_END_NAMESPACE
