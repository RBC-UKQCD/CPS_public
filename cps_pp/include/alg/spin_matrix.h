#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:13:57 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/spin_matrix.h,v 1.3 2004-06-04 21:13:57 chulwoo Exp $
//  $Id: spin_matrix.h,v 1.3 2004-06-04 21:13:57 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/spin_matrix.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// Header file for the SpinMatrix class.
//
// This file contains the declarations of the SpinMatrix 
//
//------------------------------------------------------------------


#ifndef INCLUDED_SPINMARTRIX_H
#define INCLUDED_SPINMARTRIX_H

CPS_END_NAMESPACE
#include <string.h>
#include <util/data_types.h>
CPS_START_NAMESPACE

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
    SpinMatrix(IFloat c);
    SpinMatrix(const Complex& c);
    SpinMatrix(const SpinMatrix& m);

    SpinMatrix& operator=(IFloat c);
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


CPS_END_NAMESPACE
