#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:45 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_threept/SpinMatrix.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Id: SpinMatrix.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:49:41  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:11:31  anj
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
//  Revision 1.2  2001/05/25 06:16:00  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: SpinMatrix.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_threept/SpinMatrix.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// SpinMatrix.C
//
// The spin matrix (4x4 complex) class functions.
//
// This file contains the definitions of the SpinMatrix 
// class function.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<alg/spin_matrix.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
//------------------------------------------------------------------
// The SpinMatrix class.
//------------------------------------------------------------------
//------------------------------------------------------------------

SpinMatrix::SpinMatrix()
{}

//------------------------------------------------------------------
SpinMatrix::SpinMatrix(IFloat c) 
{ *this = c; }


//------------------------------------------------------------------
SpinMatrix::SpinMatrix(const Complex& c) 
{ *this = c; }


//------------------------------------------------------------------
SpinMatrix::SpinMatrix(const SpinMatrix& m) 
{ *this = m; }


//------------------------------------------------------------------
SpinMatrix& SpinMatrix::operator = (IFloat c)
{
  this -> ZeroSpinMatrix();
  u[0] = u[10] = u[20] = u[30] = c;
  return *this;
}

//------------------------------------------------------------------
SpinMatrix& SpinMatrix::operator = (const Complex& c)
{
  this -> ZeroSpinMatrix();
  ((Complex*)u)[0] = ((Complex*)u)[5] = ((Complex*)u)[10] =
	((Complex*)u)[15] = c;
  return *this;
}


//------------------------------------------------------------------
void SpinMatrix::UnitSpinMatrix(void)
{
  IFloat *p = (IFloat *)u;

  for(int i = 0; i < 32; ++i) {
    *p++ = 0.0;
  }
  p = (IFloat *)u;
  *p = 1.0;
  *(p+10) = 1.0;
  *(p+20) = 1.0;
  *(p+30) = 1.0;
}

//------------------------------------------------------------------
void SpinMatrix::ZeroSpinMatrix(void)
{
  IFloat *p = (IFloat *)u;

  for(int i = 0; i < 32; ++i) {
    *p++ = 0.0;
  }
}


//------------------------------------------------------------------

//------------------------------------------------------------------
Complex& SpinMatrix::operator()(int i, int j)
{ return ((Complex*)u)[i*SPINS+j]; }


//------------------------------------------------------------------
const Complex& SpinMatrix::operator()(int i, int j) const
{ return ((Complex*)u)[i*SPINS+j]; }


//------------------------------------------------------------------
Complex SpinMatrix::Tr() const
{ return ((Complex*)u)[0] + ((Complex*)u)[5] + 
		((Complex*)u)[10] + ((Complex*)u)[15];
}

CPS_END_NAMESPACE
