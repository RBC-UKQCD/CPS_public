#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of class rfloat.

  $Id: rfloat.h,v 1.5 2004-12-15 07:32:07 chulwoo Exp $
 */
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-12-15 07:32:07 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/rfloat.h,v 1.5 2004-12-15 07:32:07 chulwoo Exp $
//  $Id: rfloat.h,v 1.5 2004-12-15 07:32:07 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: rfloat.h,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/rfloat.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// rfloat.h
//

#ifndef INCLUDED_RFLOAT_H
#define INCLUDED_RFLOAT_H           //!< Prevent multiple inclusion

class rfloat;
CPS_END_NAMESPACE
#include <util/enum.h>
//#include <util/data_types.h>        
CPS_START_NAMESPACE

//! A floating point type.
/*!
  This appears to be a fancy class wrapper for the IFloat type.
  Why, I do not know.
*/

class rfloat {
    IFloat x;
public:
    rfloat(IFloat a = 0);       /*!< Initialised to 0.0 by default */ 
    rfloat(const rfloat& a);    
    ~rfloat();

    rfloat& operator=(IFloat a)
    	{ x = a; return *this; }

    //! Conversion from  rfloat to IFloat.
    operator IFloat() const { return x; }
    /*!<
      Even though the user-defined conversion is a bad
      practice in many situation, it's OK here, because the
      behaviour of rfloat is exactly the same as IFloat except for
      the overloaded operators.
    */

    // Overloaded = operator to allow printf to print rfloats as IFloats:
    //IFloat operator=( const rfloat& a ) { return a.x; }

    //---------------------------------------------------------
    //  overloaded operators
    //---------------------------------------------------------
    //! overloaded binary plus
    friend rfloat operator+(const rfloat&, const rfloat&);
    //! overloaded binary plus
    friend rfloat operator+(double, const rfloat&);
    //! overloaded binary plus
    friend rfloat operator+(const rfloat&, double);

    //! overloaded binary minus
    friend rfloat operator-(const rfloat&, const rfloat&);
    //! overloaded binary minus
    friend rfloat operator-(double, const rfloat&);
    //! overloaded binary minus
    friend rfloat operator-(const rfloat&, double);

    //! overloaded binary multiply
    friend rfloat operator*(const rfloat&, const rfloat&);
    //! overloaded binary multiply
    friend rfloat operator*(double, const rfloat&);
    //! overloaded binary multiply    
    friend rfloat operator*(const rfloat&, double);

    //! overloaded binary division
    friend rfloat operator/(const rfloat&, const rfloat&);
    //! overloaded binary division
    friend rfloat operator/(double, const rfloat&);
    //! overloaded binary division
    friend rfloat operator/(const rfloat&, double);


    	//! overloaded prefix unary minus
    friend rfloat operator-(const rfloat&);

    //! overloaded sum
    rfloat& operator+=(IFloat a);
    //! overloaded sum
    rfloat& operator+=(const rfloat& a)
   	{ *this += a.x;  return *this; }

    //! overloaded subtraction
    rfloat& operator-=(IFloat a);
    //! overloaded subtraction
    rfloat& operator-=(const rfloat& a)
    	{ *this -= a.x;  return *this; }

    //! overloaded multiplication
    rfloat& operator*=(IFloat a);
    //! overloaded multiplication
    rfloat& operator*=(const rfloat& a)
    	{ *this *= a.x;  return *this; }
    //! overloaded division
    rfloat& operator/=(IFloat a);
    //! overloaded division
    rfloat& operator/=(const rfloat& a)
    	{ *this /= a.x;  return *this; }
};


//-----------------------------------------------------------------
// Free friend functions
//-----------------------------------------------------------------
rfloat operator-(const rfloat& a);

rfloat operator+(const rfloat& a, const rfloat& b);
rfloat operator+(double a, const rfloat& b);
rfloat operator+(const rfloat& a, double b);

rfloat operator-(const rfloat& a, const rfloat& b);
rfloat operator-(double a, const rfloat& b);
rfloat operator-(const rfloat& a, double b);

rfloat operator*(const rfloat& a, const rfloat& b);
rfloat operator*(double a, const rfloat& b);
rfloat operator*(const rfloat& a, double b);

rfloat operator/(const rfloat& a, const rfloat& b);
rfloat operator/(double a, const rfloat& b);
rfloat operator/(const rfloat& a, double b);





#endif

CPS_END_NAMESPACE
