#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fnone class.

  $Id: f_none.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_none/f_none.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: f_none.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:13:22  anj
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
//  $RCSfile: f_none.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_none/f_none.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_none.C
//
// Fnone is derived from Lattice. Its functions do nothing
// and return values as if there is no fermion action or
// fermion fields. The number of spin components is zero
// The site size of the fermion array FsiteSize() is
// set to 1 so that memory allocation would proceed normally.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Initialize static variables.
//------------------------------------------------------------------


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
Fnone::Fnone()
{
  cname = "Fnone";
  char *fname = "Fnone()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Fnone::~Fnone()
{
  char *fname = "~Fnone()";
  VRB.Func(cname,fname);
}

//------------------------------------------------------------------
// FclassType Fclass(void):
// It returns the type of fermion class.
//------------------------------------------------------------------
FclassType Fnone::Fclass(void){
  return F_CLASS_NONE;
}

//------------------------------------------------------------------
// int ExactFlavors() : 
// Returns the number of exact flavors of the matrix that
// is inverted during a molecular dynamics evolution.
//------------------------------------------------------------------
int Fnone::ExactFlavors(void)
{
  return 0;
}


//------------------------------------------------------------------
// int SpinComponents() : 
// Returns the number of spin components.
//------------------------------------------------------------------
int Fnone::SpinComponents(void)
{
  return 0;
}


//------------------------------------------------------------------
/*!
  \return 1, in order for memory allocation to proceed normally.
*/
//------------------------------------------------------------------
int Fnone::FsiteSize(void)
{
  return 1;
}

//------------------------------------------------------------------
// int FchkbEvl() :
// returns 1 => The fermion fields in the evolution
//      or the CG that inverts the evolution matrix
//      are defined on a single checkerboard (half the 
//      lattice).
//------------------------------------------------------------------
int Fnone::FchkbEvl(void)
{
  return 1;
}


//------------------------------------------------------------------
// int FmatEvlInv(Vector *f_out, Vector *f_in, 
//                CgArg *cg_arg, 
//                Float *true_res,
//		  CnvFrmType cnv_frm = CNV_FRM_YES):
// It does nothing and returns 0.
//------------------------------------------------------------------
int Fnone::FmatEvlInv(Vector *f_out, Vector *f_in, 
		      CgArg *cg_arg, 
		      Float *true_res,
		      CnvFrmType cnv_frm)
{
  char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  // Return the number of iterations
  return 0;
}


//------------------------------------------------------------------
// Overloaded function is same as original but with true_res=0;
//------------------------------------------------------------------
int Fnone::FmatEvlInv(Vector *f_out, Vector *f_in, 
		      CgArg *cg_arg, 
		      CnvFrmType cnv_frm)
{ return FmatEvlInv(f_out, f_in, cg_arg, 0, cnv_frm); }


//------------------------------------------------------------------
// int FmatInv(Vector *f_out, Vector *f_in, 
//             CgArg *cg_arg, 
//             Float *true_res,
//             CnvFrmType cnv_frm = CNV_FRM_YES,
//             PreserveType prs_f_in = PRESERVE_YES):
// It does nothing and returns 0.
//------------------------------------------------------------------
int Fnone::FmatInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm,
		   PreserveType prs_f_in)
{
  char *fname = "FmatInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);
  
  // Return the number of iterations
  return 0;
}


//------------------------------------------------------------------
// Overloaded function is same as original but with true_res=0;
//------------------------------------------------------------------
int Fnone::FmatInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   CnvFrmType cnv_frm,
		   PreserveType prs_f_in)
{ return FmatInv(f_out, f_in, cg_arg, 0, cnv_frm, prs_f_in); }


//------------------------------------------------------------------
// int FeigSolv(Vector **f_eigenv, Float *lambda, int valid_eig[],
//              EigArg *eig_arg, 
//              CnvFrmType cnv_frm = CNV_FRM_YES):
//------------------------------------------------------------------
int Fnone::FeigSolv(Vector **f_eigenv, Float lambda[], 
		    Float chirality[], int valid_eig[],
		    Float **hsum,
		    EigArg *eig_arg, 
		    CnvFrmType cnv_frm)
{
  char *fname = "FeigSolv(V*,F*,I*,E)";
  VRB.Func(cname,fname);

  return 0;
}

//------------------------------------------------------------------
// SetPhi(Vector *phi, Vector *frm1, Vector *frm2, Float mass):
//! Does nothing.
//------------------------------------------------------------------
void Fnone::SetPhi(Vector *phi, Vector *frm1, Vector *frm2, 
		   Float mass){
  char *fname = "SetPhi(V*,V*,V*,F)";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// FforceSite(Matrix& force, Vector *frm, int *x, int mu):
// It calculates the fermion force at site x and direction mu.
// frm is the fermion field that resulted from the application
// of the inverter on the pseudofermion field.
//------------------------------------------------------------------
void Fnone::FforceSite(Matrix& force, Vector *frm, int *x, int mu)
{
  char *fname = "FforceSite(M&,V*,i*,i)";
  VRB.Func(cname,fname);

  force.ZeroMatrix();
}


//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *frm, Float mass, 
//                 Float step_size):
// It evolves the canonical momentum mom by step_size
// using the fermion force.
//------------------------------------------------------------------
void Fnone::EvolveMomFforce(Matrix *mom, Vector *frm, 
			    Float mass, Float step_size){
  char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Float FhamiltonNode(Vector *phi, Vector *chi):
// The fermion Hamiltonian of the node sublattice.
// chi must be the solution of Cg with source phi.	       
//------------------------------------------------------------------
Float Fnone::FhamiltonNode(Vector *phi, Vector *chi){
  char *fname = "FhamiltonNode(V*,V*)";
  VRB.Func(cname,fname);

  return 0.0;
}


//------------------------------------------------------------------
// Float BhamiltonNode(Vector *boson, Float mass):
// The boson Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float Fnone::BhamiltonNode(Vector *boson, Float mass){
  char *fname = "BhamiltonNode(V*,F)";
  VRB.Func(cname,fname);

  //???
  return 0.0;
}




//------------------------------------------------------------------
// int FsiteOffsetChkb(const int *x):
// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is not the canonical one but it is particular
// to the fermion type. x[i] is the 
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
//------------------------------------------------------------------
int Fnone::FsiteOffsetChkb(const int *x) const {
// ???
  ERR.NotImplemented(cname, "FsiteOffsetChkb");
  return 0; 
}


//------------------------------------------------------------------
// int FsiteOffset(const int *x):
// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is the canonical one. X[I] is the
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
//------------------------------------------------------------------
int Fnone::FsiteOffset(const int *x) const {
// ???
  ERR.NotImplemented(cname, "FsiteOffset");
  return 0; 
}


CPS_END_NAMESPACE
