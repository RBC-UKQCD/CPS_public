#include <config.h>
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

#include <alg/spin_matrix.h>

CPS_START_NAMESPACE

//------------------------------------------------------------------
// The SpinMatrix class.
//------------------------------------------------------------------

//------------------------------------------------------------------
SpinMatrix::SpinMatrix(IFloat c) { *this = c; }

//------------------------------------------------------------------
SpinMatrix::SpinMatrix(const Complex& c) { *this = c; }

//------------------------------------------------------------------
SpinMatrix& SpinMatrix::operator=(IFloat c) {
  ZeroSpinMatrix();
  u[2*(0+0*SPINS)] = u[2*(1+1*SPINS)] = u[2*(2+2*SPINS)] = u[2*(3+3*SPINS)] = c;
  return *this;
}

//------------------------------------------------------------------
SpinMatrix& SpinMatrix::operator=(const Complex& c) {
  ZeroSpinMatrix();
  ((Complex*)u)[0+0*SPINS] = 
	((Complex*)u)[1+1*SPINS] = 
	((Complex*)u)[2+2*SPINS] = 
	((Complex*)u)[3+3*SPINS] = c;
  return *this;
}


//------------------------------------------------------------------
void SpinMatrix::UnitSpinMatrix(void) {
  IFloat *p = (IFloat*)u;

  for(int i=0; i<2*SPINS*SPINS; i++) *p++ = 0.0;
  p = (IFloat *)u;
  p[2*(0+0*SPINS)] = 1.0;
  p[2*(1+1*SPINS)] = 1.0;
  p[2*(2+2*SPINS)] = 1.0;
  p[2*(3+3*SPINS)] = 1.0;
}

//------------------------------------------------------------------
Complex SpinMatrix::Tr() const
{ return ((Complex*)u)[0+0*SPINS] + ((Complex*)u)[1+1*SPINS] + 
	((Complex*)u)[2+2*SPINS] + ((Complex*)u)[3+3*SPINS];
}

CPS_END_NAMESPACE
