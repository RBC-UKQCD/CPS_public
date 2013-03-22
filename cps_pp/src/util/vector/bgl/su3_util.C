#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Utility routines for SU(3) matrices.

  $Id: su3_util.C,v 1.6 2013-03-22 05:40:05 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-03-22 05:40:05 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/bgl/su3_util.C,v 1.6 2013-03-22 05:40:05 chulwoo Exp $
//  $Id: su3_util.C,v 1.6 2013-03-22 05:40:05 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/bgl/su3_util.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <math.h>
#include <util/vector.h>
CPS_START_NAMESPACE


//---------------------------------------------------------------//
// pre-calculated 1/3
//---------------------------------------------------------------//
IFloat Matrix::inv3 = 1./3.;


//---------------------------------------------------------------//
//  Matrix
//---------------------------------------------------------------//
/*!
  \param a A linear array representation of a 3x3 complex matrix, such that 
  real part of the (i,j) element is at array position [6i+2j] 
  and the imaginary part of the (i,j) element is at array position [6i+2j+1].
  \post This matrix is the hermitian conjugate of \a m.

  \a a must not be an alias of this matrix.
 */
void Matrix::Dagger(const IFloat* a)
{
    u[0]  = a[0];   u[1]  = -a[1];
    u[6]  = a[2];   u[7]  = -a[3];
    u[12] = a[4];   u[13] = -a[5];
    u[2]  = a[6];   u[3]  = -a[7];
    u[8]  = a[8];   u[9]  = -a[9];
    u[14] = a[10];  u[15] = -a[11];
    u[4]  = a[12];  u[5]  = -a[13];
    u[10] = a[14];  u[11] = -a[15];
    u[16] = a[16];  u[17] = -a[17];
}

/*!
  \param dag A matrix \a A.
  \post This matrix is set to\n
  <em>1/2(M-A) - 1/6 Trace M-A)</em>
  \n where \a M is the original value of this matrix.
*/
void Matrix::TrLessAntiHermMatrix(const Matrix &dag)
{
    // get 1/2(A - dag(A)) =  1/2A - dag(1/2A)
    *this -= dag;

    IFloat *p = (IFloat *)u;
    vecTimesEquFloat(p, 0.5, 18);

    IFloat c = inv3 * (*(p+1) + *(p+9) + *(p+17));
    *(p+1) -= c;
    *(p+9) -= c;
    *(p+17) -= c;
}

/*!
  \param v1 A complex 3-vector \a u
  \param v2 A complex 3-vector \a v
  \post This matrix is assigned the values M(i,j) = 2 u(i)v(j)*
*/
void Matrix::Cross2(const Vector& v1, const Vector& v2)
{
    const IFloat *p1 = (const IFloat *)v1.v;
    const IFloat *p2 = (const IFloat *)v2.v;
    IFloat *p = (IFloat *)u;
    for(int i = 0; i < 3; ++i ){
	int m1 = 2*i;
	int n1 = m1+1;
        for(int j = 0; j < 3; ++j ) {
	    int m = 2*(3*i+j);
	    int n = m+1;
	    int m2 = 2*j;
	    int n2 = m2+1;

	    *(p+m) = 2.0*(*(p1+m1) * *(p2+m2) + *(p1+n1) * *(p2+n2));
	    *(p+n) = 2.0*(*(p1+n1) * *(p2+m2) - *(p1+m1) * *(p2+n2));
	}
    }
}

/*!
  \param a an array of 8 real numbers.
  \post This matrix is assigned the value\n
  <em>i L<sub>i</sub> * a[i]</em>
  \n where <em>L<sub>i</sub></em> is the \a i 'th Gell-Mann matrix.
*/
void Matrix::AntiHermMatrix(const IFloat *a)
{
    IFloat s3 = 0.57735 * *(a+7);       // 1/sqrt(3) = 0.57735;
    IFloat *p = (IFloat *)u;

    *p = 0.0;
    *(p+8) = 0.0;
    *(p+16) = 0.0;

    *(p+1) = *(a+2) + s3;
    *(p+9) = -*(a+2) + s3;
    *(p+17) = -2.0 * s3;

    *(p+2) = *(a+1);
    *(p+3) = *a;
    *(p+4) = *(a+4);
    *(p+5) = *(a+3);
    *(p+6) = -*(a+1);
    *(p+7) = *a;
    *(p+10) = *(a+6);
    *(p+11) = *(a+5);
    *(p+12) = -*(a+4);
    *(p+13) = *(a+3);
    *(p+14) = -*(a+6);
    *(p+15) = *(a+5);
}


/*!
  \param q An array of length at least 2.
  \post The real and imaginarey parts of the determinant of this matrix
  are written to the 0 and 1 elements of \a q.
*/
void Matrix::Det(IFloat* q) const
{
    const IFloat *p = (const IFloat *)u;

    //------------------------------------------------------------
    //  (re0, im0) = u[4]*u[8] - u[5]*u[7];
    //------------------------------------------------------------
    IFloat re0 =   *(p+8)  * *(p+16) - *(p+9)  * *(p+17)
    	        - *(p+10) * *(p+14) + *(p+11) * *(p+15) ;
    IFloat im0 =   *(p+8)  * *(p+17) + *(p+9)  * *(p+16)
    		- *(p+10) * *(p+15) - *(p+11) * *(p+14) ;

    //------------------------------------------------------------
    //  (re1, im1) = u[5]*u[6] - u[3]*u[8];
    //------------------------------------------------------------
    IFloat re1 =   *(p+10) * *(p+12) - *(p+11) * *(p+13)
    		- *(p+6)  * *(p+16) + *(p+7)  * *(p+17) ;
    IFloat im1 =   *(p+10) * *(p+13) + *(p+11) * *(p+12)
    		- *(p+6)  * *(p+17) - *(p+7)  * *(p+16) ;

    //------------------------------------------------------------
    //  (re2, im2) = u[3]*u[7] - u[4]*u[6];
    //------------------------------------------------------------
    IFloat re2 =   *(p+6)  * *(p+14) - *(p+7)  * *(p+15)
    		- *(p+8)  * *(p+12) + *(p+9)  * *(p+13) ;
    IFloat im2 =   *(p+6)  * *(p+15) + *(p+7)  * *(p+14)
    		- *(p+8)  * *(p+13) - *(p+9)  * *(p+12) ;

    *q = *p * re0 - *(p+1) * im0 + *(p+2) * re1 - *(p+3) * im1
         + *(p+4) * re2 - *(p+5) * im2 ;
    *(q+1) = *p * im0 + *(p+1) * re0 + *(p+2) * im1 + *(p+3) * re1
            + *(p+4) * im2 + *(p+5) * re2 ;
}

/*!
  This method is specifically for an anti-hermitian matrix only.
*/
IFloat Matrix::NegHalfTrSquare() const
{
    IFloat *p = (IFloat *)u;

    return *(p+3) * *(p+3) + *(p+2) * *(p+2) + 
           (*(p+1)- *(p+9))*(*(p+1)- *(p+9))*0.25 +
	   *(p+5) * *(p+5) + *(p+4) * *(p+4) +
	   *(p+11) * *(p+11) + *(p+10) * *(p+10) +
	   *(p+17) * *(p+17)*0.75 ;    // 0.75 = (sqrt(3)/2)^2
}


//---------------------------------------------------------------//
//  Vector
//---------------------------------------------------------------//
/*!
  \post The L2 norm of this vector is 1
*/
void Vector::Normalize()
{
    IFloat norm = dotProduct((IFloat *)v, (IFloat *)v, 6);

    if( !(norm == 1.0) ){
        norm = 1.0/sqrt(norm);
	vecTimesEquFloat((IFloat *)v, norm, 6);
    }
}


//
//	v2' = v2 - v1 * (v1, v2)
// 	then	(v1, v2') = 0
// 
/*!
  \param p1 A vector
  \post This vector is orthogonalised against the vector \a p1 so
  that its inner product with p1 vanishes.
*/
void Vector::Orthogonalize(const Vector& p1)
{
    IFloat re = 0.0;
    IFloat im = 0.0;
    int i;
    IFloat *q = (IFloat *)v;
    const IFloat *q1 = (const IFloat *)p1.v;

    for(i = 0; i < 3; ++i) {
	int m = 2*i;
	int n = m+1;
        re += *(q1+m) * *(q+m) + *(q1+n) * *(q+n);
	im += *(q1+m) * *(q+n) - *(q1+n) * *(q+m);
    }
    if( !(re == 0.0 && im == 0.0) ) {
        for(i = 0; i < 3; ++i){
	    int m = 2*i;
	    int n = m+1;
	    *(q+m) -= re * *(q1+m) - im * *(q1+n);
	    *(q+n) -= re * *(q1+n) + im * *(q1+m);
	}
    }
}








CPS_END_NAMESPACE
