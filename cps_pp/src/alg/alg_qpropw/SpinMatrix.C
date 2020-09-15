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
#include <util/error.h>

CPS_START_NAMESPACE

//------------------------------------------------------------------
// The SpinMatrix class.
//------------------------------------------------------------------

#ifndef INLINE_SPIN_MATRIX
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
#endif


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


SpinMatrix SpinMatrix::One()
{
    SpinMatrix ret;
    ret.UnitSpinMatrix();
    return ret;
}

SpinMatrix SpinMatrix::Gamma5()
{
    SpinMatrix ret;
    ret.ZeroSpinMatrix();
    ret(0, 0) = ret(1, 1) = 1.0;
    ret(2, 2) = ret(3, 3) = -1.0;
    return ret;
}

SpinMatrix SpinMatrix::Gamma(int mu)
{
    SpinMatrix ret;
    ret.ZeroSpinMatrix();

    switch (mu) {
	case 0:
	    ret(0, 3) = ret(1, 2) = Rcomplex(0, 1);
	    ret(2, 1) = ret(3, 0) = Rcomplex(0, -1);
	    break;

	case 1:
	    ret(0, 3) = ret(3, 0) = -1;
	    ret(1, 2) = ret(2, 1) = 1;
	    break;

	case 2:
	    ret(0, 2) = ret(3, 1) = Rcomplex(0, 1);
	    ret(1, 3) = ret(2, 0) = Rcomplex(0, -1);
	    break;

	case 3:
	    ret(0, 2) = ret(1, 3) = ret(2, 0) = ret(3, 1) = 1;
	    break;

	default:
	    ERR.General("SpinMatrix", "Gamma", "Unknown mu = %d\n", mu);
	    break;
    }

    return ret;
}

SpinMatrix SpinMatrix::GammaMuGamma5(int mu)
{
    return Gamma(mu) * Gamma5();
}

SpinMatrix SpinMatrix::Sigma(int mu, int nu)
{
    return (Gamma(mu) * Gamma(nu) - Gamma(nu) * Gamma(mu)) * Rcomplex(0.5, 0);
}

SpinMatrix SpinMatrix::SigmaGamma5(int mu, int nu)
{
    return Sigma(mu, nu) * Gamma5();
}


SpinMatrix SpinMatrix::OnePlusGamma(int mu)
{
    return One() + Gamma(mu);
}

SpinMatrix SpinMatrix::OneMinusGamma(int mu)
{
    return One() - Gamma(mu);
}

CPS_END_NAMESPACE
