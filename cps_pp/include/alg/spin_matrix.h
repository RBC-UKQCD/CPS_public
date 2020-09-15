// -*- mode:c++; c-basic-offset:4 -*-
//------------------------------------------------------------------
//
// Header file for the SpinMatrix class.
//
// This file contains the declarations of the SpinMatrix 
//
//------------------------------------------------------------------


#ifndef INCLUDED_SPINMARTRIX_H
#define INCLUDED_SPINMARTRIX_H

#include <string.h>
#include <util/data_types.h>
CPS_START_NAMESPACE
class SpinMatrix;

#define INLINE_SPIN_MATRIX



//------------------------------------------------------------------
// The SpinMatrix class.
//------------------------------------------------------------------
class SpinMatrix
{
    enum{SPINS=4};

  public:
   Float u[2*SPINS*SPINS];	// The matrix
    // CREATORS
    SpinMatrix() {}
    SpinMatrix(Float c){ *this = c; }
    SpinMatrix(const Complex& c){ *this = c; }
    SpinMatrix(Float *c){
	for(int i =0;i<2*SPINS*SPINS;i++)
	u[i] = *(c+i);
    }

#ifndef INLINE_SPIN_MATRIX
    SpinMatrix& operator=(Float c);
    SpinMatrix& operator=(const Complex& c);
#else
SpinMatrix& operator=(IFloat c) {
  ZeroSpinMatrix();
  u[2*(0+0*SPINS)] = u[2*(1+1*SPINS)] = u[2*(2+2*SPINS)] = u[2*(3+3*SPINS)] = c;
  return *this;
}
SpinMatrix& operator=(const Complex& c) {
  ZeroSpinMatrix();
  ((Complex*)u)[0+0*SPINS] = 
	((Complex*)u)[1+1*SPINS] = 
	((Complex*)u)[2+2*SPINS] = 
	((Complex*)u)[3+3*SPINS] = c;
  return *this;
}
#endif

    void ZeroSpinMatrix(void) {
        for(int i=0; i<2*SPINS*SPINS; i++) u[i] = 0;
    }
    void UnitSpinMatrix(void);

    // ACCESSORS
    Complex& operator()(int i, int j){ return ((Complex*)u)[i*SPINS+j]; }
    const Complex& operator()(int i, int j) const { return ((Complex*)u)[i*SPINS+j]; }
    inline Complex& operator[](int i) { return ((Complex*)u)[i]; }
    inline const Complex& operator[](int i) const { return ((Complex*)u)[i]; }
   

    SpinMatrix operator*(const Complex &c)const {
        SpinMatrix ret;
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                ret(i, j) = operator()(i, j) * c;
            }
        }
        return ret;
    }

    const SpinMatrix &operator*=(const Complex &c) {
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                operator()(i, j) *= c;
            }
        }
        return *this;
    }

    SpinMatrix operator+(const SpinMatrix &s)const {
        SpinMatrix ret;
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                ret(i, j) = operator()(i, j) + s.operator()(i, j);
            }
        }
        return ret;
    }

    const SpinMatrix &operator+=(const SpinMatrix &s) {
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                operator()(i, j) += s.operator()(i, j);
            }
        }
        return *this;
    }

    SpinMatrix operator-(const SpinMatrix &s)const {
        SpinMatrix ret;
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                ret(i, j) = operator()(i, j) - s.operator()(i, j);
            }
        }
        return ret;
    }

    const SpinMatrix &operator-=(const SpinMatrix &s) {
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                operator()(i, j) -= s.operator()(i, j);
            }
        }
        return *this;
    }

    SpinMatrix operator*(const SpinMatrix &s)const {
        SpinMatrix ret(0.0);
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                for(int k = 0; k < 4; ++k) {
                    ret(i, j) += operator()(i, k) * s.operator()(k, j);
                }
            }
        }
        return ret;
    }




    Complex Tr() const;

    SpinMatrix transpose(void)const {
        SpinMatrix ret;
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                ret(j, i) = operator()(i, j);
            }
        }
        return ret;
    }

    static SpinMatrix One();
    static SpinMatrix Gamma5();
    static SpinMatrix Gamma(int mu);
    static SpinMatrix GammaMuGamma5(int mu);
    static SpinMatrix Sigma(int mu, int nu);
    static SpinMatrix SigmaGamma5(int mu, int nu);
    static SpinMatrix OnePlusGamma(int mu);
    static SpinMatrix OneMinusGamma(int mu);
};

static inline Rcomplex Trace(const SpinMatrix &a, const SpinMatrix &b) {
    Rcomplex ret = 0;
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            ret += a.operator()(i, j) * b.operator()(j, i);
        }
    }
    return ret;
}

// Tr [A * B^T]
static inline Rcomplex TraceTranspose(const SpinMatrix &a,
                                      const SpinMatrix &b)
{
    Rcomplex ret = 0;
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            ret += a.operator()(i, j) * b.operator()(i, j);
        }
    }
    return ret;
}

CPS_END_NAMESPACE
#endif

