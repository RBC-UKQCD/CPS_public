#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Alg class methods.
  
  $Id: alg_base.C,v 1.2 2003-07-24 16:53:53 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_base/alg_base.C,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: alg_base.C,v 1.2 2003-07-24 16:53:53 zs Exp $
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
//  $Revision: 1.2 $
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
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE

/*!
  \param latt The lattice and actions with which the algorithm is to act.
  \param c_arg The common argument structure for all algorithms, apparently.
*/
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
