#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of Double64 class/type.

  $Id: double64.h,v 1.4 2004-08-18 11:57:36 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:36 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/double64.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Id: double64.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/double64.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//
// WARNING WARNING WARNING
//-------------------------
//
// Do not use this file for other purposes except for the glb_sum
// If you want to use it for other applications you must 
// allocate memory and set the pointer  D64CRAM below.

/****************************************************************
 * Double64, written by Roy 					*
 ****************************************************************/

/* header file to support 64-bit (Double64) mathematical functions */

#ifndef INCLUDED_DOUBLE64_H_
#define INCLUDED_DOUBLE64_H_

//Get the global options:
CPS_END_NAMESPACE
CPS_START_NAMESPACE

//! A type used (solely?) to accumulate global sums in double precision.
/*!
  On platforms other than QCDSP, this should default to \c double.
  On QCDSP this is a custom class implementing 64-bit arithmetic.
*/

#ifdef GLOBALSUM_TYPE

typedef GLOBALSUM_TYPE Double64;

//Otherwise, on QCDSP, use the special 64bit double implementation:
#else

CPS_END_NAMESPACE
#include <flotar.h>		// float64
CPS_START_NAMESPACE

class Double64;

extern "C" double double64round(const Double64 *);

class Double64 {
   float64 dvalue;
 public:
   //-------------------------------------------------------------
   // CTORs
   //-------------------------------------------------------------
   Double64() {}
   Double64(const Double64& dbl) { _cpyf64(&dvalue, &(dbl.dvalue)); }
   Double64(const float64& f64) { _cpyf64(&dvalue, &f64); }
   Double64(const double d) { convf64f32(&dvalue, d); }
   Double64(const int i) { convf64i32(&dvalue, i); }


   //-------------------------------------------------------------
   // For setting bits explicitly:
   //-------------------------------------------------------------
   Double64(const int eu,  const double du) {
	dvalue.ext = eu; dvalue.flt = du;
   } 


   //-------------------------------------------------------------
   // Conversions
   //-------------------------------------------------------------
   //   operator double() const {return dvalue.flt;} /* floor */
   operator double() const {return double64round(this);} 


   //-------------------------------------------------------------
   // Arithmetic
   //-------------------------------------------------------------
   Double64& operator+=(const Double64& x)
   	{ addf64(&dvalue, &dvalue, &(x.dvalue));  return *this; }
   Double64& operator+=(const double x)
   	{ addf32f64(&dvalue, x, &dvalue);  return *this; }

   Double64& operator-=(const Double64& x)
   	{ subf64(&dvalue, &dvalue, &(x.dvalue));  return *this; }
   Double64& operator-=(const double x)
   	{ subf64f32(&dvalue, &dvalue, x);  return *this; }

   Double64& operator*=(const Double64& x)
   	{ mpyf64(&dvalue, &dvalue, &(x.dvalue));  return *this; }
   Double64& operator*=(const double x)
   	{ mpyf32f64(&dvalue, x, &dvalue);  return *this; }

   Double64& operator/=(const Double64& x)
   	{ divf64(&dvalue, &dvalue, &(x.dvalue));  return *this; }
   Double64& operator/=(const double x)
   	{ divf64f32(&dvalue, &dvalue, x);  return *this; }

   Double64& operator=(const Double64& x)
   	{ _cpyf64(&dvalue, &(x.dvalue)); return *this; }
   Double64& operator=(const double d)
   	{ convf64f32(&dvalue, d);  return *this; }


   //-------------------------------------------------------------
   // Friend arithmetic functions
   //-------------------------------------------------------------
   friend Double64 operator-(const Double64&); 	/* negation */

   friend Double64 operator+(const Double64&, const Double64&);
   friend Double64 operator+(const Double64&, const double);
   friend Double64 operator+(const double, const Double64&);
   
   friend Double64 operator-(const Double64&, const Double64&);
   friend Double64 operator-(const Double64&, const double);
   friend Double64 operator-(const double, const Double64&);
   
   friend Double64 operator*(const Double64&, const Double64&);
   friend Double64 operator*(const Double64&, const double);
   friend Double64 operator*(const double, const Double64&);
   
   friend Double64 operator/(const Double64&, const Double64&);
   friend Double64 operator/(const Double64&, const double);
   friend Double64 operator/(const double, const Double64&);
   
   friend Double64 fabs(const Double64&);

   //-------------------------------------------------------------
   // Friend comparisons
   //-------------------------------------------------------------
   friend int operator==(const Double64& left,  const Double64& right)
       { return left.dvalue.flt == right.dvalue.flt &&
                left.dvalue.ext == right.dvalue.ext; }
   friend int operator==(const Double64& left,   const double right)
       { return left.dvalue.flt == right && left.dvalue.ext == 0; }
   friend int operator==(const double  right, const Double64& left)
       { return left.dvalue.flt == right && left.dvalue.ext == 0; }

   friend int operator>=(const Double64& left,  const Double64& right)
       { return left.dvalue.flt > right.dvalue.flt ||
	       (left.dvalue.flt == right.dvalue.flt &&
	        left.dvalue.ext >= right.dvalue.ext); }
   friend int operator>=(const Double64& left,   const double right)
       { return left.dvalue.flt >= right; }
   friend int operator>=(const double  left, const Double64& right)
       { return left > right.dvalue.flt ||
	       (left == right.dvalue.flt && right.dvalue.ext == 0) ; }

   friend int operator>(const Double64& left,  const Double64& right)
       { return left.dvalue.flt > right.dvalue.flt ||
	       (left.dvalue.flt == right.dvalue.flt &&
	        left.dvalue.ext > right.dvalue.ext); }
   friend int operator>(const Double64& left,   const double right)
       { return left.dvalue.flt > right ||
	       (left.dvalue.flt == right &&
	        left.dvalue.ext != 0); }
   friend int operator>(const double  left, const Double64& right)
       { return left > right.dvalue.flt; }

   friend int operator<(const Double64& left,  const Double64& right)
       { return left.dvalue.flt < right.dvalue.flt ||
	       (left.dvalue.flt == right.dvalue.flt &&
	        left.dvalue.ext < right.dvalue.ext); }
   friend int operator<(const Double64& left,   const double right)
       { return left.dvalue.flt < right; }
   friend int operator<(const double  left, const Double64& right)
       { return left < right.dvalue.flt ||
	       (left == right.dvalue.flt &&
	        right.dvalue.ext != 0) ; }

   friend int operator<=(const Double64& left,  const Double64& right)
       { return left.dvalue.flt < right.dvalue.flt ||
	       (left.dvalue.flt == right.dvalue.flt &&
	        left.dvalue.ext <= right.dvalue.ext); }
   friend int operator<=(const Double64& left,   const double right)
       { return left.dvalue.flt < right ||
	       (left.dvalue.flt == right &&
	        left.dvalue.ext == 0); }
   friend int operator<=(const double  left, const Double64& right)
       { return left <= right.dvalue.flt; }

   /* compare functions return 1 if the left operand is greater than
      right operand, -1 if left operand is less than the right operand,
      and 0 if the operands are equal.	*/
   friend int compare(const Double64& l, const Double64& r)
       { return cmpf64(&(l.dvalue), &(r.dvalue)); }
   friend int compare(const Double64& l, const double r)
       { return cmpf64f32(&(l.dvalue), r); }
   friend int compare(const double l,  const Double64& r)
       { return -cmpf64f32(&(r.dvalue), l); }



   friend Double64 acos(const Double64& x);	/* inverse cosine */
   friend Double64 acosh(const Double64& x);
   					/* inverse hyperbolic cosine */
   friend Double64 acot(const Double64& x);	/* inverse cotangent */
   friend Double64 acot2(const Double64& x, const Double64& y);
					/* inverse cotangent (2 args);*/
   friend Double64 acoth(const Double64& x);
   				/* inverse hyperbolic cotangent */

   friend Double64 asin(const Double64& x);	/* inverse sine */
   friend Double64 asinh(const Double64& x);
   					/* inverse hyperbolic sine */
   friend Double64 atan(const Double64& x);	/* inverse tangent */
   friend Double64 atan2(const Double64& x, const Double64& y);
					/* inverse tangent (2 args); */
   friend Double64 atanh(const Double64& x);
   					/* inverse hyperbolic tangnet */

   friend long ceili(const Double64& x);	/* ceiling functions */
   friend unsigned long ceilu(const Double64& x);
   friend double ceilf(const Double64& x);

   friend Double64 cos(const Double64& x);	/* cosine */
   friend Double64 cosh(const Double64& x);	/* hyperbolic cosine */
   friend Double64 cot(const Double64& x);	/* cotangent */
   friend Double64 coth(const Double64& x);/* hyperbolic cotangent */

   friend Double64 exp(const Double64& x);/* e raised to some power */
   friend Double64 expl(const Double64& x);/* e raised to some power */

   friend Double64 inv(const Double64& x);/* computes 1/argument */

   friend Double64 log(const Double64& x);/* natural logarithm */
   friend Double64 log10(const Double64& x);	/* base 10 logarithm */

   friend Double64 pow(const Double64& x, const Double64& y); /* x**y */

   friend double roundf(const Double64& x);
   					/* round to nearest integer */
   friend long roundi(const Double64& x);
   friend unsigned long roundu(const Double64& x);

   friend Double64 rsqrt(const Double64& x);/* reciprocal square root */

   friend Double64 sin(const Double64& x);	/* sine */
   friend Double64 sinh(const Double64& x);	/* hyperbolic sine */
   friend Double64 sqr(const Double64& x);/* square the argument */
   friend Double64 sqrt(const Double64& x);
   				/* square root of the argument */

   friend Double64 tan(const Double64& x);	/* tangent */
   friend Double64 tanh(const Double64& x);	/* hyperbolic tangent */
   friend Double64 tento(const Double64& x);
   					/* ten raised to some power */

   friend double truncf(const Double64& x);
   				/* fractional part removed, if any */
   friend long trunci(const Double64& x);
   friend unsigned long truncu(const Double64& x);


} /*class Double64 */;


//-------------------------------------------------------------------
// make a scratch area in CRAM
//-------------------------------------------------------------------
extern float64 *const D64CRAM;	/* = (float64 *)0x80fbe */


//-------------------------------------------------------------
// Friend arithmetic functions
//-------------------------------------------------------------
inline Double64 operator-(const Double64& x) 	/* negation */
    { negf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }

inline Double64 operator+(const Double64& l, const Double64& r)
    { addf64(D64CRAM, &(l.dvalue), &(r.dvalue)); return *D64CRAM; }

inline Double64 operator+(const Double64& l, const double r)
    { addf32f64(D64CRAM, r, &(l.dvalue)); return *D64CRAM; }

inline Double64 operator+(const double l, const Double64& r)
    { addf32f64(D64CRAM, l, &(r.dvalue)); return *D64CRAM; }


inline Double64 operator-(const Double64& l, const Double64& r)
    { subf64(D64CRAM, &(l.dvalue), &(r.dvalue)); return *D64CRAM; }
inline Double64 operator-(const Double64& l, const double r)
    { subf64f32(D64CRAM, &(l.dvalue), r); return *D64CRAM; }
inline Double64 operator-(const double l, const Double64& r)
    { subf32f64(D64CRAM, l, &(r.dvalue)); return *D64CRAM; }


inline Double64 operator*(const Double64& l, const Double64& r)
    { mpyf64(D64CRAM, &(l.dvalue), &(r.dvalue)); return *D64CRAM; }
inline Double64 operator*(const Double64& l, const double r)
    { mpyf32f64(D64CRAM, r, &(l.dvalue)); return *D64CRAM; }
inline Double64 operator*(const double l, const Double64& r)
    { mpyf32f64(D64CRAM, l, &(r.dvalue)); return *D64CRAM; }

inline Double64 operator/(const Double64& l, const Double64& r)
    { divf64(D64CRAM, &(l.dvalue), &(r.dvalue)); return *D64CRAM; }
inline Double64 operator/(const Double64& l, const double r)
    { divf64f32(D64CRAM, &(l.dvalue), r); return *D64CRAM; }
inline Double64 operator/(const double l, const Double64& r)
    { divf32f64(D64CRAM, l, &(r.dvalue)); return *D64CRAM; }

inline Double64 fabs(const Double64& x)
    { absf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }



//-------------------------------------------------------------------
// Trigonometric Functions
//-------------------------------------------------------------------
inline Double64 acos(const Double64& x)
    { acosf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 acosh(const Double64& x)
    { acoshf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 acot(const Double64& x)
    { acotf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 acot2(const Double64& x, const Double64& y)
    {acot2f64(D64CRAM, &(x.dvalue), &(y.dvalue)); return *D64CRAM;}
				
inline Double64 acoth(const Double64& x)
    { acothf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }

inline Double64 asin(const Double64& x)
    { asinf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 asinh(const Double64& x)
    { asinhf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 atan(const Double64& x)
    { atanf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 atan2(const Double64& x, const Double64& y)
    {atan2f64(D64CRAM, &(x.dvalue), &(y.dvalue)); return *D64CRAM;}
				
inline Double64 atanh(const Double64& x)
    { atanhf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }

inline long ceili(const Double64& x)
    { return ceili32f64(&(x.dvalue)); }
inline unsigned long ceilu(const Double64& x)
    { return ceilu32f64(&(x.dvalue)); }
inline double ceilf(const Double64& x)
    { return ceilf32f64(&(x.dvalue)); }

inline Double64 cos(const Double64& x)
    { cosf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 cosh(const Double64& x)
    { coshf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 cot(const Double64& x)
    { cotf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 coth(const Double64& x)
    { cothf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }

inline Double64 exp(const Double64& x)
    { expf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 expl(const Double64& x)
    { expf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }

inline Double64 inv(const Double64& x)
    { invf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }

inline Double64 log(const Double64& x)
    { logf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 log10(const Double64& x)
    { log10f64(D64CRAM, &(x.dvalue)); return *D64CRAM; }

inline Double64 pow(const Double64& x, const Double64& y)
    { powf64(D64CRAM, &(x.dvalue), &(y.dvalue)); return *D64CRAM; }

inline double roundf(const Double64& x)
    { return roundf32f64(&(x.dvalue)); }
inline long roundi(const Double64& x)
    { return roundi32f64(&(x.dvalue)); }
inline unsigned long roundu(const Double64& x)
    { return roundu32f64(&(x.dvalue)); }

inline Double64 rsqrt(const Double64& x)
    { rsqrtf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }

inline Double64 sin(const Double64& x)
    { sinf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 sinh(const Double64& x)
    { sinhf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 sqr(const Double64& x)
    { sqrf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 sqrt(const Double64& x)
    { sqrtf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
				

inline Double64 tan(const Double64& x)
    { tanf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 tanh(const Double64& x)
    { tanhf64(D64CRAM, &(x.dvalue)); return *D64CRAM; }
inline Double64 tento(const Double64& x)
    { tentof64(D64CRAM, &(x.dvalue)); return *D64CRAM; }

inline double truncf(const Double64& x)
    { return truncf32f64(&(x.dvalue)); }
inline long trunci(const Double64& x)
    { return trunci32f64(&(x.dvalue)); }
inline unsigned long truncu(const Double64& x)
    { return truncu32f64(&(x.dvalue)); }



#endif /* End of Double64 precision choice */

#endif /* _INCLUDED_DOUBLE64_H_ */

CPS_END_NAMESPACE
