#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_nn.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: lattice_nn.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:13:30  anj
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
//  Revision 1.2  2001/05/25 06:16:10  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: lattice_nn.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_types/lattice_nn.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_nn.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// Gnone and Fnone.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// No gauge action -- no fermion action
//------------------------------------------------------------------
GnoneFnone::GnoneFnone()
{
  cname = "GnoneFnone";
  char *fname = "GnoneFnone()";
  VRB.Func(cname,fname);

  //???
}

GnoneFnone::~GnoneFnone()
{
  char *fname = "~GnoneFnone()";
  VRB.Func(cname,fname);

  //???
}


CPS_END_NAMESPACE
