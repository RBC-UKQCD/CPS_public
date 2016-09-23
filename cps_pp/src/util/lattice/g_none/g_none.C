#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Gnone class.

  $Id: g_none.C,v 1.7 2006/04/13 18:21:51 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006/04/13 18:21:51 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/g_none/g_none.C,v 1.7 2006/04/13 18:21:51 chulwoo Exp $
//  $Id: g_none.C,v 1.7 2006/04/13 18:21:51 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.7 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/g_none/g_none.C,v $
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
ForceArg Gnone::EvolveMomGforce(Matrix *mom, Float step_size){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);

  return ForceArg(0.0,0.0,0.0);
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
