#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:45 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_base/alg_base.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Id: alg_base.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:11:21  anj
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
//  $RCSfile: alg_base.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_base/alg_base.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_base.C
//
// Alg is the base abstract class
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>	// exit()
#include<alg/alg_base.h>
#include<alg/common_arg.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE


Alg::Alg(Lattice & latt, 
	 CommonArg *c_arg) :
	 alg_lattice(latt)
{
  cname = "Alg";
  char *fname = "Alg(L&,CommonArg*)";
  VRB.Func(cname,fname);
  
  // Set the common argument pointer
  //----------------------------------------------------------------
  if(c_arg == 0)
    ERR.Pointer(cname,fname, "common_arg");
  common_arg = c_arg;
}


Alg::~Alg() {
  char *fname = "~Alg()";
  VRB.Func(cname,fname);

  //???
}


Lattice& Alg::AlgLattice()
{
  return alg_lattice;
}



CPS_END_NAMESPACE
