#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of Matrix class method Matrix::Unitarize

  $Id: unitarize.C,v 1.4 2004-08-18 11:58:08 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/comsrc/unitarize.C,v 1.4 2004-08-18 11:58:08 zs Exp $
//  $Id: unitarize.C,v 1.4 2004-08-18 11:58:08 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/comsrc/unitarize.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/vector.h>
#include <math.h>	// sqrt()
CPS_START_NAMESPACE

//! Utility routine used by Matrix::Unitarize
/*!\todo Since these routines are only used by a Matrix class method they 
  could be made private methods of the Matrix class.
*/
 void normalize(Float *p);
//! Utility routine used by Matrix::Unitarize
void orthogonalize(Float *p2, const Float *p1); 
//! Utility routine used by Matrix::Unitarize
void crossProductConj(Float *v3, const Float *v1, const Float *v2);

//------------------------------------------------------------------
// The unitarize routine
//------------------------------------------------------------------
/*!
  \post This matrix is a member of the SU(3) group.
*/
void Matrix::Unitarize()
{
    Float *p1 = &u[0];	// row 1
    Float *p2 = &u[6] ;	// row 2
    Float *p3 = &u[12];	// row 3

    //------------------------------------------------------------
    //  first step: Normalize the first row.
    //------------------------------------------------------------
    normalize(p1);

    //------------------------------------------------------------
    // 2nd step: make the second row orthogonal to the first row
    //  	 and Normalize it
    //
    //		p2' = p2 - p1 * (p1^*, p2)
    //	then	(p1^*, p2') = 0
    //------------------------------------------------------------
    orthogonalize(p2, p1);
    normalize(p2);

    //------------------------------------------------------------
    // the conjugate of last row is the cross product of the first 2
    // p3 = (p1^* X p2^*)
    //------------------------------------------------------------

    crossProductConj(p3, p1, p2);
}


//------------------------------------------------------------------
// Few routines that are needed by the Unitarize routine
//------------------------------------------------------------------
void normalize(Float *p)
{
    Float norm = 0;
    int i;
    for(i = 0; i < 6; i++){
      norm += p[i] * p[i];
    }
    if( !(norm == 1.0) ){
        norm = 1.0/sqrt(norm);
	for(i = 0; i < 6; i++)
	  p[i] *= norm;
    }
}

//	v2' = v2 - v1 * (v1^*, v2)
// 	then	(v1^*, v2') = 0


void orthogonalize(Float *p2, const Float *p1)
{
    Float re = 0.0;
    Float im = 0.0;
    int i;

    for(i = 0; i < 3; ++i) {
	int m = 2*i;
	int n = m+1;
        re += *(p1+m) * *(p2+m) + *(p1+n) * *(p2+n);
	im += *(p1+m) * *(p2+n) - *(p1+n) * *(p2+m);
    }
    if( !(re == 0.0 && im == 0.0) ) {
        for(i = 0; i < 3; ++i){
	    int m = 2*i;
	    int n = m+1;
	    *(p2+m) -= re * *(p1+m) - im * *(p1+n);
	    *(p2+n) -= re * *(p1+n) + im * *(p1+m);
	}
    }
}


// v3 = ( v1 x v2 )^*
void crossProductConj(Float *v3, const Float *v1, const Float *v2)
{
  v3[0] = v1[2]*v2[4] - v1[3]*v2[5] -v1[4]*v2[2] + v1[5]*v2[3];
  v3[1] = -v1[2]*v2[5] - v1[3]*v2[4] + v1[4]*v2[3] + v1[5]*v2[2];

  v3[2] = -v1[0]*v2[4] + v1[1]*v2[5] +v1[4]*v2[0] - v1[5]*v2[1];
  v3[3] = v1[0]*v2[5] + v1[1]*v2[4] - v1[4]*v2[1] - v1[5]*v2[0];

  v3[4] = v1[0]*v2[2] - v1[1]*v2[3] -v1[2]*v2[0] + v1[3]*v2[1];
  v3[5] = -v1[0]*v2[3] - v1[1]*v2[2] + v1[2]*v2[1] + v1[3]*v2[0];
}



CPS_END_NAMESPACE
