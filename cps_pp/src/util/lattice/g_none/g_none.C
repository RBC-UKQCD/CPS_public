#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Gnone class.

  $Id: g_none.C,v 1.3 2003-10-23 13:38:59 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-10-23 13:38:59 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_none/g_none.C,v 1.3 2003-10-23 13:38:59 zs Exp $
//  $Id: g_none.C,v 1.3 2003-10-23 13:38:59 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2.6.1  2003/09/23 21:39:20  zs
//  Gnone needs an AllStaple method
//
//  Revision 1.2  2003/07/24 16:53:54  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.2  2001/06/19 18:13:27  anj
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
//  $RCSfile: g_none.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/g_none/g_none.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// g_none.C
//
// Gnone is derived from Lattice. Its functions act
// as if there is no gauge action i.e. beta = 0.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE



//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
Gnone::Gnone()
{
  cname = "Gnone";
  char *fname = "Gnone()";
  VRB.Func(cname,fname);

  //???
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Gnone::~Gnone()
{
  char *fname = "~Gnone()";
  VRB.Func(cname,fname);

  //???
}


//------------------------------------------------------------------
// GclassType Gclass(void):
// It returns the type of gauge class
//------------------------------------------------------------------
GclassType Gnone::Gclass(void){
  return G_CLASS_NONE;
}


//------------------------------------------------------------------
/*!
  In this case, of course, the force is zero.
  \param force The computed force from the gauge action.
  \param x the lattice site coordinates.
  \param mu The direction mu.
  \todo Should this not be be a virtual Lattice method?
*/
//------------------------------------------------------------------
void Gnone::GforceSite(Matrix& force, int *x, int mu)
{
  char *fname = "GforceSite(M&,i*,i)";
  VRB.Func(cname,fname);
  force.ZeroMatrix();
}


//------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float step_size):
// It evolves the canonical momentum mom by step_size
// using the pure gauge force.
//------------------------------------------------------------------
void Gnone::EvolveMomGforce(Matrix *mom, Float step_size){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Float GhamiltonNode(void):
// The pure gauge Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float Gnone::GhamiltonNode(void){
  char *fname = "GhamiltonNode()";
  VRB.Func(cname,fname);

  return 0.0;
}


//-------------------------------------------------------------------------//
//void GactionGradient(Matrix &grad, int *x, int mu)
//   Calculates the partial derivative of the gauge action
//   w.r.t. the link U_mu(x).
//----------------------------------------------------------------------------
void Gnone::GactionGradient(Matrix &grad, int *x, int mu)
{
  ERR.NotImplemented(cname, "GactionGradient(M&,i*,i)") ;
}

void Gnone::AllStaple(Matrix &grad, const int *x, int mu)
{
  ERR.NotImplemented(cname, "AllStaple(M&,i*,i)") ;
}

CPS_END_NAMESPACE
