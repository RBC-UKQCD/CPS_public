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

//------------------------------------------------------------------
//------------------------------------------------------------------
// The SpinMatrix class.
//------------------------------------------------------------------
//------------------------------------------------------------------

SpinMatrix::SpinMatrix()
{}

//------------------------------------------------------------------
SpinMatrix::SpinMatrix(Float c) 
{ *this = c; }


//------------------------------------------------------------------
SpinMatrix::SpinMatrix(const Complex& c) 
{ *this = c; }


//------------------------------------------------------------------
SpinMatrix::SpinMatrix(const SpinMatrix& m) 
{ *this = m; }


//------------------------------------------------------------------
SpinMatrix& SpinMatrix::operator = (Float c)
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
  Float *p = (Float *)u;

  for(int i = 0; i < 32; ++i) {
    *p++ = 0.0;
  }
  p = (Float *)u;
  *p = 1.0;
  *(p+10) = 1.0;
  *(p+20) = 1.0;
  *(p+30) = 1.0;
}

//------------------------------------------------------------------
void SpinMatrix::ZeroSpinMatrix(void)
{
  Float *p = (Float *)u;

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

