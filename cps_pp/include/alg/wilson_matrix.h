#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/wilson_matrix.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: wilson_matrix.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
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
//  $RCSfile: wilson_matrix.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/wilson_matrix.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// wilson_matrix.h
//
// Header file WilsonMatrix class for 
// Wilson-like quarks.
//
// For now this is specific to three colors. The constructor
// will exit if the number of colors is not equal to three.
//
//------------------------------------------------------------------


#ifndef INCLUDED_WILSONMATRIX_H
#define INCLUDED_WILSONMATRIX_H

CPS_END_NAMESPACE
#include <math.h>
#include <string.h>
#include<util/rcomplex.h>
#include<alg/spin_matrix.h>
#include<util/vector.h>
#include<util/wilson.h>
CPS_START_NAMESPACE

typedef struct { Float real; Float imag; } complex;
typedef struct { Rcomplex c[3]; } su3_vector;
typedef struct { su3_vector d[4]; } wilson_vector;
typedef struct { wilson_vector c[3]; } color_wilson_vector;
typedef struct { color_wilson_vector d[4]; } wilson_matrix;


class WilsonMatrix
{
    wilson_matrix p;

  public:

    // CONSTRUCTORS
    WilsonMatrix();
    WilsonMatrix(const WilsonMatrix& rhs);
    WilsonMatrix(const wilson_matrix& rhs);
    WilsonMatrix(const Float& rhs);
    WilsonMatrix(const Rcomplex& rhs);
    WilsonMatrix(int source_spin, int source_color, const wilson_vector&);
    WilsonMatrix(int source_spin, int source_color, int sink_spin, 
	int sink_color, const Rcomplex&);

    void hconj(); // hermitean conjugate the WilsonMatrix
    WilsonMatrix& gl(int dir); //mult the prop by gamma_dir on the left
    WilsonMatrix& gr(int dir); //mult the prop by gamma_dir on the left
    WilsonMatrix conj_cp(); // make a copy of the hermitean conjugate
    wilson_vector& sol(int source_spin, int source_color); // get a sol. vector
    void load_vec(int source_spin, int source_color, const wilson_vector&);
    Rcomplex Trace();
    const wilson_matrix& wmat() const; // get p 

    // operator functions
    WilsonMatrix& operator=(const WilsonMatrix& rhs);
    WilsonMatrix& operator=(const wilson_matrix& rhs);
    WilsonMatrix& operator=(const Float& rhs);
    const WilsonMatrix& operator[](int i);
    WilsonMatrix& operator+=(const WilsonMatrix& rhs);
    WilsonMatrix& operator-=(const WilsonMatrix& rhs);
    WilsonMatrix& operator*=(const WilsonMatrix& rhs);
    WilsonMatrix& operator*=(const Float& rhs);
    WilsonMatrix& operator*=(const Rcomplex& rhs);

    friend Rcomplex Trace(const WilsonMatrix& p1, const WilsonMatrix& p2);

};

// some proto-types for functions that operate on WilsonMatrices
extern WilsonMatrix operator*(const WilsonMatrix& lhs, const WilsonMatrix& rhs);
	// times operator
extern WilsonMatrix operator*(const Float& num, const WilsonMatrix& mat);
	// times operator
extern WilsonMatrix operator*(const WilsonMatrix& mat, const Float& num);
	// times operator
extern WilsonMatrix operator*(const Rcomplex& num, const WilsonMatrix& mat);
	// times operator
extern WilsonMatrix operator*(const WilsonMatrix& mat, const Rcomplex& num);
	// times operator
extern WilsonMatrix operator+(const WilsonMatrix& lhs, const WilsonMatrix& rhs);
	// plus operator
extern WilsonMatrix operator-(const WilsonMatrix& lhs, const WilsonMatrix& rhs);
	// minus operator
extern void mult_by_gamma_left( int dir, const wilson_matrix& src, wilson_matrix& dest );
	// left multiply by gamma_dir
extern void mult_by_gamma_right(int dir, const wilson_matrix& src, wilson_matrix& dest );
	// right multiply by gamma_dir
extern Rcomplex Trace(const WilsonMatrix& p1, const WilsonMatrix& p2);
	// Spin and Color trace of a WilsonMatrix
extern Matrix SpinTrace(const WilsonMatrix& Wmat); 
	// Spin trace of a WilsonMatrix
extern Matrix SpinTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2); 
	// Spin trace of two WilsonMatrices
extern Matrix SpinTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2,
						const WilsonMatrix& Wmat3); 
	// Spin trace of two WilsonMatrices
extern SpinMatrix ColorTrace(const WilsonMatrix& Wmat);
	// Color trace of a WilsonMatrix
extern SpinMatrix ColorTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2);
	// Color trace of three WilsonMatrices
extern SpinMatrix ColorTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2,
						const WilsonMatrix& Wmat3);
	// Color trace of three WilsonMatrices
extern Rcomplex Tr(const Matrix& a, const Matrix& b);
        // trace of two (color) Matrices
extern Rcomplex Tr(const SpinMatrix& a, const SpinMatrix& b);
        // trace of two SpinMatrices


#endif

CPS_END_NAMESPACE
