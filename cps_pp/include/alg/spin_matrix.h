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

class SpinMatrix;


enum{SPINS=4};

//------------------------------------------------------------------
// The SpinMatrix class.
//------------------------------------------------------------------
class SpinMatrix
{
    Float u[2*SPINS*SPINS];	// The matrix


  public:
    // CREATORS
    SpinMatrix();
    SpinMatrix(Float c);
    SpinMatrix(const Complex& c);
    SpinMatrix(const SpinMatrix& m);

    SpinMatrix& operator=(Float c);
    SpinMatrix& operator=(const Complex& c);

    void ZeroSpinMatrix(void);
    void UnitSpinMatrix(void);

    // ACCESSORS
    Complex& operator()(int i, int j);
    const Complex& operator()(int i, int j) const;
    Complex& operator[](int i) { return ((Complex*)u)[i]; }
    const Complex& operator[](int i) const { return ((Complex*)u)[i]; }
    Complex Tr() const;
};


#endif

