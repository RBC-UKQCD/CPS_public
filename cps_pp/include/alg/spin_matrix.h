#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/spin_matrix.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: spin_matrix.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:49:42  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:11:33  anj
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
//  $RCSfile: spin_matrix.h,v $
//  $Revision: 1.1.1.1 $
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
#include<util/data_types.h>
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
