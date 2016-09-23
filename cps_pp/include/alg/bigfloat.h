//------------------------------------------------------------------
/*!\file
  \brief  Definitions of the bigfloat wrapper class.

  $Id: bigfloat.h,v 1.9 2007/11/07 05:24:41 chulwoo Exp $
*/
//------------------------------------------------------------------
#include<config.h>
#include <gmp.h>

#ifdef USE_MPFR
#include <mpfr.h>
#include <mpf2mpfr.h>
#endif

CPS_START_NAMESPACE

#ifndef INCLUDED_BIGFLOAT_H
#define INCLUDED_BIGFLOAT_H

//------------------------------------------------------------------
//
// Simple C++ wrapper for multiprecision datatype used for Remez
// algorithm
//
//------------------------------------------------------------------
//! Arbitrary precision arithmetic.
/*!
  This is used by the AlgRemez class
  It is a wrapper around the GNU
  Multiprecision Library (GMP).GMP library.
  and therefore can only be used if the CPS is built with GMP.
  This is achieved at CPS configure time by using the option
  --enable-gmp=GMP_PATH with the configure script, where GMP_PATH is the 
  directory where GMP has been installed. 	
*/
class bigfloat {

private:

  mpf_t x;

public:

  bigfloat() { mpf_init(x); }
  bigfloat(const bigfloat& y) { mpf_init_set(x, y.x); }
  bigfloat(const unsigned long u) { mpf_init_set_ui(x, u); }
  bigfloat(const long i) { mpf_init_set_si(x, i); }
  bigfloat(const int i) {mpf_init_set_si(x,(long)i);}
  bigfloat(const float d) { mpf_init_set_d(x, (double)d); }
  bigfloat(const double d) { mpf_init_set_d(x, d); }  
  ~bigfloat(void) { mpf_clear(x); }
  operator const Float (void) const { return (Float)mpf_get_d(x); }
  static void setDefaultPrecision(unsigned long dprec) {
    unsigned long bprec =  (unsigned long)(3.321928094 * (double)dprec);
    mpf_set_default_prec(bprec);
  }

  void setPrecision(unsigned long dprec) {
    unsigned long bprec =  (unsigned long)(3.321928094 * (double)dprec);
    mpf_set_prec(x,bprec);
  }
  
  unsigned long getPrecision(void) const { return mpf_get_prec(x); }

  unsigned long getDefaultPrecision(void) const { return mpf_get_default_prec(); }

  bigfloat& operator=(const bigfloat& y) {
    mpf_set(x, y.x); 
    return *this;
  }

  bigfloat& operator=(const unsigned long y) { 
    mpf_set_ui(x, y);
    return *this; 
  }
  
  bigfloat& operator=(const signed long y) {
    mpf_set_si(x, y); 
    return *this;
  }
  
  bigfloat& operator=(const float y) {
    mpf_set_d(x, (double)y); 
    return *this;
  }

  bigfloat& operator=(const double y) {
    mpf_set_d(x, y); 
    return *this;
  }

  size_t write(void);
  size_t read(void);

  /* Arithmetic Functions */

  bigfloat& operator+=(const bigfloat& y) { return *this = *this + y; }
  bigfloat& operator-=(const bigfloat& y) { return *this = *this - y; }
  bigfloat& operator*=(const bigfloat& y) { return *this = *this * y; }
  bigfloat& operator/=(const bigfloat& y) { return *this = *this / y; }

  friend bigfloat operator+(const bigfloat& x, const bigfloat& y) {
    bigfloat a;
    mpf_add(a.x,x.x,y.x);
    return a;
  }

  friend bigfloat operator+(const bigfloat& x, const unsigned long y) {
    bigfloat a;
    mpf_add_ui(a.x,x.x,y);
    return a;
  }

  friend bigfloat operator-(const bigfloat& x, const bigfloat& y) {
    bigfloat a;
    mpf_sub(a.x,x.x,y.x);
    return a;
  }
  
  friend bigfloat operator-(const unsigned long x, const bigfloat& y) {
    bigfloat a;
    mpf_ui_sub(a.x,x,y.x);
    return a;
  }
  
  friend bigfloat operator-(const bigfloat& x, const unsigned long y) {
    bigfloat a;
    mpf_sub_ui(a.x,x.x,y);
    return a;
  }

  friend bigfloat operator-(const bigfloat& x) {
    bigfloat a;
    mpf_neg(a.x,x.x);
    return a;
  }

  friend bigfloat operator*(const bigfloat& x, const bigfloat& y) {
    bigfloat a;
    mpf_mul(a.x,x.x,y.x);
    return a;
  }

  friend bigfloat operator*(const bigfloat& x, const unsigned long y) {
    bigfloat a;
    mpf_mul_ui(a.x,x.x,y);
    return a;
  }

  friend bigfloat operator/(const bigfloat& x, const bigfloat& y){
    bigfloat a;
    mpf_div(a.x,x.x,y.x);
    return a;
  }

  friend bigfloat operator/(const unsigned long x, const bigfloat& y){
    bigfloat a;
    mpf_ui_div(a.x,x,y.x);
    return a;
  }

  friend bigfloat operator/(const bigfloat& x, const unsigned long y){
    bigfloat a;
    mpf_div_ui(a.x,x.x,y);
    return a;
  }

  friend bigfloat sqrt_bf(const bigfloat& x){
    bigfloat a;
    mpf_sqrt(a.x,x.x);
    return a;
  }

  friend bigfloat sqrt_bf(const unsigned long x){
    bigfloat a;
    mpf_sqrt_ui(a.x,x);
    return a;
  }

  friend bigfloat abs_bf(const bigfloat& x){
    bigfloat a;
    mpf_abs(a.x,x.x);
    return a;
  }

  friend bigfloat pow_bf(const bigfloat& a, long power) {
    bigfloat b;
    mpf_pow_ui(b.x,a.x,power);
    return b;
  }

#ifdef USE_MPFR
  friend bigfloat pow_bf(const bigfloat& a, const bigfloat &power) {
    bigfloat b;
    mpfr_pow(b.x,a.x,power.x,GMP_RNDN);
    return b;
  }
#endif

  /* Comparison Functions */

  friend int operator>(const bigfloat& x, const bigfloat& y) {
    int test;
    test = mpf_cmp(x.x,y.x);
    if (test > 0) return 1;
    else return 0;
  }

  friend int operator<(const bigfloat& x, const bigfloat& y) {
    int test;
    test = mpf_cmp(x.x,y.x);
    if (test < 0) return 1;
    else return 0;
  }

  friend int sgn(const bigfloat&);

  /* Miscellaneous Functions */

//  friend bigfloat& random(void);
};

#endif

CPS_END_NAMESPACE
