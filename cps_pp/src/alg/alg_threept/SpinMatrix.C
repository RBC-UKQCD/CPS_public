#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:40 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_threept/SpinMatrix.C,v 1.6 2004-08-18 11:57:40 zs Exp $
//  $Id: SpinMatrix.C,v 1.6 2004-08-18 11:57:40 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: SpinMatrix.C,v $
//  $Revision: 1.6 $
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
#include <alg/spin_matrix.h>
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
