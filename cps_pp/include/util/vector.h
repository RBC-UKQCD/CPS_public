#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/vector.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: vector.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2002/12/04 17:16:27  zs
//  Merged the new 2^4 RNG into the code.
//  This new RNG is implemented in the LatRanGen class.
//  The following algorithm and utility classes are affected:
//
//  AlgEig                  Fdwf
//  AlgGheatBath            Fstag
//  AlgHmd                  GlobalJobParameter
//  AlgNoise                Lattice
//  AlgPbp                  Matrix
//  AlgThreept              RandomGenerator
//                          Vector
//
//  Revision 1.4  2001/08/16 10:50:31  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:19  anj
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
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: vector.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/vector.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// vector.h
//
// Header file for the vector class.
//
// For now this is specific to three colors. The constructor
// will exit if the number of colors is not equal to three.
//
// This file contains the declarations of the Vector and Matrix 
// classes as well as the declarations of some genaral c-style
// functions that perform operations on vectors of general
// length.
//
//------------------------------------------------------------------


#ifndef INCLUDED_VECTOR_H
#define INCLUDED_VECTOR_H

CPS_END_NAMESPACE
#include <string.h>
#include<util/data_types.h>
CPS_START_NAMESPACE


class Vector; 	// forward declaration
class Matrix;


//------------------------------------------------------------------
// Declarations of some genaral c-style functions that perform
// operations on vectors. For these functions there exists 
// optimized assembly code.
//------------------------------------------------------------------
extern "C" 
{
void moveMem(void *b, const void *a, int len); // b=a, REQUIRE:len >= 4

void mDotMEqual(IFloat* c, const IFloat* a, const IFloat* b); // C = A*B
                                         // Assumption: c!=a and c!=b 

void mDotMPlus(IFloat* c, const IFloat* a, const IFloat* b); // C += A*B
                                         // Assumption: c!=a and c!=b 

void uDotXEqual(IFloat* y, const IFloat* m, const IFloat* x); // y = M*x

IFloat dotProduct(const IFloat *a, const IFloat *b, int);  // len >= 3

void vecAddEquVec(IFloat *a, const IFloat *b, int); 	// A += B

void vecMinusEquVec(IFloat *a, const IFloat *b, int); 	// A -= B

void vecNegative(IFloat *a, const IFloat *b, int); 	// A = -B

void vecTimesEquFloat(IFloat *a, IFloat b, int); // A *= b

void fTimesV1PlusV2(IFloat *a, IFloat b, const IFloat *c,
			const IFloat *d, int size); 	// A = b*C+D

void fTimesV1MinusV2(IFloat *a, IFloat b, const IFloat *c,
			const IFloat *d, int size); 	// A = b*C-D

void compDotProduct(IFloat *c_r, IFloat *c_i, 
        	    const IFloat *a, const IFloat *b, int);  // len >= 3

void cTimesV1PlusV2(IFloat *a, IFloat b_re, IFloat b_im, const IFloat *c,
                    const IFloat *d, int size);      // A = b*C+D

void cTimesV1MinusV2(IFloat *a, IFloat b_re, IFloat b_im, const IFloat *c,
	             const IFloat *d, int size);      // A = b*C-D

void oneMinusfTimesMatrix(IFloat *a, IFloat b, const IFloat *c,
			int n);     // A = 1-b*C, A and C are Matrix

}


//------------------------------------------------------------------
// Declarations of some genaral c-style functions that perform
// operations on vectors. For these functions there is
// no optimized assembly.
//------------------------------------------------------------------

void uDotXMinus(IFloat* y, const IFloat* u, const IFloat* x);

void uDagDotXEqual(IFloat* y, const IFloat* u, const IFloat* x);

void uDagDotXPlus(IFloat* y, const IFloat* u, const IFloat* x);

//------------------------------------------------------------------
// Declarations of some genaral c-style functions that compute
// re/im parts of the su3 characters of various representations
// of su3 matrices.
//------------------------------------------------------------------

extern IFloat reChar6(IFloat *p) ;
extern IFloat imChar6(IFloat *p) ;
extern IFloat reChar8(IFloat *p) ;
extern IFloat reChar10(IFloat *p) ;
extern IFloat imChar10(IFloat *p) ;


//------------------------------------------------------------------
// For now the number of colors for the Matrix and Vector classes
// is fixed. If the number of colors is not 3 the appropriate 
// constructors will exit.
//------------------------------------------------------------------
enum{COLORS=3};

//------------------------------------------------------------------
// The Matrix class.
// For now Matrix is a class of general 3x3 complex matrices.
// If the number of colors is not 3 the constructor will exit.
//------------------------------------------------------------------
class Matrix
{
    Float u[2*COLORS*COLORS];	// The matrix

    friend class Vector;

  public:
    // CREATORS
    Matrix();
    Matrix(IFloat c);
    Matrix(const Complex& c);
    Matrix(const Matrix& m);

    Matrix& operator=(IFloat c);
    Matrix& operator=(const Complex& c);

    Matrix& operator=(const Matrix& m)
    { moveMem(u, m.u, COLORS*COLORS*2*sizeof(Float)); 
      return *this; }

    // MANIPULATORS
    Matrix& operator+=(const Matrix& m)
    { vecAddEquVec((IFloat *)u, (IFloat *)m.u, COLORS*COLORS*2);
      return *this; }

    Matrix& operator+=(IFloat c)
     { u[0] += c;  u[8] += c;  u[16] += c;  return *this; }

    Matrix& operator-=(const Matrix& m)
    { vecMinusEquVec((IFloat *)u, (IFloat *)m.u, COLORS*COLORS*2);
      return *this; }

    Matrix& operator-=(IFloat c)
      { u[0] -= c;  u[8] -= c;  u[16] -= c;  return *this; }

    Matrix& operator*=(IFloat c)
    { vecTimesEquFloat((IFloat *)u, c, COLORS*COLORS*2); return *this; }

    void DotMEqual(const Matrix& a, const Matrix& b)
    { mDotMEqual((IFloat *)u, (IFloat *) a.u, (IFloat *) b.u);}
       //u = a.u * b.u

    void DotMPlus(const Matrix& a, const Matrix& b)
    { mDotMPlus((IFloat *)u, (IFloat *)a.u, (IFloat *)b.u);}
       //u += a.u * b.u

    void Trans(const IFloat* m);
       // *this = transpose of m.                 Assume: this != m

    void Trans(const Matrix& m)
        { Trans((const IFloat *)(m.u)); }

    void Dagger(const IFloat* m);
    	// put dag(m) and save it to this object

    void Dagger(const Matrix& m)
    	{ Dagger((const IFloat *)&m); }
   
    void TrLessAntiHermMatrix(const Matrix& this_dag);
    	// translate this to traceless anti-hermitian matrix
	// passing this_dag as an argument in oder to make it
	// more efficient.

    void Cross2(const Vector& v1, const Vector& v2);
    	// *this = 2 * (V1 x V2^*)

    void AntiHermMatrix(const IFloat *a);
    	// a points to an array of 8 real numbers
	    // *this = i \lamda^i * a^i
	    // \lambda^i are 8 Gellmann matrices

    void Unitarize(void);
    	// force this matrix to be a unitary SU(3) matrix

    void UnitMatrix(void);

    void ZeroMatrix(void);

    void NegMatrix(const Matrix& m)
    	{ vecNegative((IFloat *)u, (IFloat *)&m, COLORS*COLORS*2); }

    void OneMinusfTimesM(IFloat x, const Matrix& m)
    	{ oneMinusfTimesMatrix((IFloat *)u, x, (IFloat *)&m,
	  COLORS*COLORS*2); }


    // ACCESSORS
    Complex& operator()(int i, int j);

    const Complex& operator()(int i, int j) const;

    Complex& operator[](int i) { return ((Complex*)u)[i]; }

    const Complex& operator[](int i) const { return ((Complex*)u)[i]; }

    void Det(IFloat* c) const;
    	// get determinant and store it to c

    IFloat ReTr() const;		// real part of the trace

    Complex Tr() const;

    IFloat NegHalfTrSquare() const;
    	// -0.5 tr(*this)^2, this matrix must be anti-hermitian

    IFloat ErrorSU3() const;
    	// return |U~dag U - I|^2

    // SU(3) Characters

    Complex Char3() const { return Tr() ; } ;

    Complex Char6() const ;

    Complex Char8() const ;

    Complex Char10() const ;

};



//------------------------------------------------------------------
// The Vector class.
// For now Vector is a class of general 3 component column
// vectors. If the number of colors is not 3 the constructor
// will exit.
// Vector is not automatically normalized when constructed.
//------------------------------------------------------------------
class Vector
{
    Float v[2*COLORS];	// Vector

    friend class Matrix;

  public:
    // CREATORS
    Vector();

    Vector& operator=(const Vector& x)
    { moveMem(v, x.v, COLORS*2*sizeof(Float)); 
      return *this; }

    // MANIPULATORS
    Vector& operator*=(IFloat p)
    { vecTimesEquFloat((IFloat *)v, p, COLORS*2); return *this; }

    Vector& operator+=(const Vector& x)
    { vecAddEquVec((IFloat *)v, (IFloat *)x.v, COLORS*2);
      return *this; }

    Vector& operator-=(const Vector& x)
    { vecMinusEquVec((IFloat *)v, (IFloat *)x.v, COLORS*2);
      return *this; }

    void DotXEqual(const Matrix& m, const Vector& x)
    { uDotXEqual((IFloat *)v, (IFloat *) m.u, (IFloat *) x.v); }
       // v = m.u * x.v, m should be in CRAM, x MUST be in DRAM */

    void Normalize(void);
       // Normalize v to 1.

    void Orthogonalize(const Vector& p1);



    //--------------------------------------------------------------
    // Functions that act on arrays of vectors of general length.
    // The array of vectors is treated as an array of IFloating 
    // numbers pointed to by &v and having length len.
    // This set of functions does not really fit the way
    // Vector is currently defined (as an array with re/im and 
    // color indeces only. It extends the notion of Vector to
    // a general array of IFloating numbers.
    //--------------------------------------------------------------

    void CopyVec(const Vector *b, int len)
    { moveMem(&v, b, len*sizeof(Float)); }
       // len is the number of real numbers in the array.
       // Copy b to v

    Float NormSqNode(int len)
    {return Float( dotProduct((IFloat *)&v, (IFloat *)&v, len) ); }
       // len is the number of real numbers in the array.
       // Returns the dot product (v*, v) on the node.

    Float NormSqGlbSum(int len);
       // len is the number of real numbers in the array.
       // Returns the dot product (v*, v) summed over all nodes.

    Float ReDotProductNode(const Vector *b, int len)
    {return Float( dotProduct((IFloat *)&v, (IFloat *)b, len) ); }
       // len is the number of real numbers in the array.
       // Returns the real part of the dot product (v, b)

    Float ReDotProductGlbSum(const Vector *b, int len);
       // len is the number of real numbers in the array.
       // Returns the real part of the dot product (v, b) 
       // summed over all nodes.

    void Vector::NormSqArraySliceSum(Float *f_out, const int size, const int dir);
       // Slice sum norm square the vector  v  of element size "size" perp. 
       // to direction dir

    void Vector::SliceArraySum(Float *sum, const Float *f_in, const int dir);
       // Slice sum a lattice of Floats f_in to form an array of length L[dir]
       // lattice size of virtual grid in direction dir

    void Vector::SliceArraySumFive(Float *sum, const Float *f_in, const int dir);
       // 5D Slice sum a lattice of Floats f_in to form an array of length L[dir]
       // lattice size of virtual grid in direction dir

    void VecNegative(const Vector *b, int len)
    {vecNegative((IFloat *)&v, (IFloat *)b, len);}
       // len is the number of real numbers in the array.
       //  v = -v 

    void VecTimesEquFloat(const Float &fb, int len)
    {vecTimesEquFloat((IFloat *)&v, IFloat(fb), len);}
       // len is the number of real numbers in the array.
       //  v = v * fb 

    void VecAddEquVec(const Vector *b, int len)
    { vecAddEquVec((IFloat *)&v, (IFloat *)b, len);}
       // len is the number of real numbers in the array.
       // v = v + b

    void VecMinusEquVec(const Vector *b, int len)
    { vecMinusEquVec((IFloat *)&v, (IFloat *)b, len);}
       // v = v - b

    void FTimesV1PlusV2(const Float &fb, const Vector *c,
			const Vector *d, int len)
    { fTimesV1PlusV2((IFloat *)&v, IFloat(fb), (IFloat *)c, 
                        (IFloat *)d, len); }
       // len is the number of real numbers in the array.
       // v = fb * c + d

    void FTimesV1MinusV2(const Float &fb, const Vector *c,
			 const Vector *d, int len)
    { fTimesV1MinusV2((IFloat *)&v, IFloat(fb), (IFloat *)c, 
                        (IFloat *)d, len); }
       // len is the number of real numbers in the array.
       // v = fb * c - d

    Complex CompDotProductNode(const Vector *b, int len)
    {
       IFloat c_r, c_i;
       compDotProduct(&c_r, &c_i, (IFloat *)&v, (IFloat *)b, len);
       return Complex(c_r,c_i);
     }
       // len is the number of complex numbers in the array.
       // Returns the complex part of the dot product (v*, b)

    Complex CompDotProductGlbSum(const Vector *b, int len);
       // len is the number of real numbers in the array.
       // Returns the complex part of the dot product (v*, b) 
       // summed over all nodes.

    void CTimesV1PlusV2(const Complex &fb, const Vector *c,
                         const Vector *d, int len)
    { cTimesV1PlusV2((IFloat *)&v, real(fb), imag(fb), (IFloat *)c, 
	              (IFloat *)d, len); }
       // len is the number of real numbers in the array.
       // v = fb * c + d
 
    void CTimesV1MinusV2(const Complex &fb, const Vector *c,
                           const Vector *d, int len)
    { cTimesV1MinusV2((IFloat *)&v, real(fb), imag(fb), (IFloat *)c, 
	                (IFloat *)d, len); }
       // len is the number of real numbers in the array.
       // v = fb * c - d

};




#endif

CPS_END_NAMESPACE
