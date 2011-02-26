#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of complex floating point data type.

  $Id: rcomplex.h,v 1.5 2011-02-26 00:19:27 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2011-02-26 00:19:27 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/rcomplex.h,v 1.5 2011-02-26 00:19:27 chulwoo Exp $
//  $Id: rcomplex.h,v 1.5 2011-02-26 00:19:27 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: rcomplex.h,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/rcomplex.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// rcomplex.h

#ifndef INCLUDED_RCOMPLEX_H
#define INCLUDED_RCOMPLEX_H

class Rcomplex;
CPS_END_NAMESPACE
#include <util/data_types.h>
CPS_START_NAMESPACE

//! A class defining a complex floating point data type.
/*! The usual operations, arithmetic and otherwise, on complex numbers
  are defined.

  \e N.B. the friend functions do not actually reference private data or methods.
*/

class Rcomplex {
  IFloat re;
  IFloat im;

  public:
  /*!
    The complex number is initialised with zero real and imaginary
    parts by default.
  */
  Rcomplex(IFloat a = 0, IFloat b = 0);
  Rcomplex(const Rcomplex& a);
  ~Rcomplex();

  Rcomplex& operator=(const Rcomplex& a);

  //! Get the real part
  /*! \return the real part */
  IFloat real() const { return re; }
  //! Get the imaginary part
  /*! \return the imaginary part */
  IFloat imag() const { return im; }

  IFloat& realx() { return re; }
  IFloat& imagx() { return im; }

  //! Get the imaginary part
  /*! \return the imaginary part */
  
  //! Assign to the real part
  /*! \param r the real part */
  void real(IFloat r) { re = r; }
  //! Assign to the imaginary part
  /*! \param i the imaginary part */
  void imag(IFloat i) { im = i; }

  //! Returns the square norm of the complex number.
  IFloat norm() const;
  //! Returns the norm (modulus) of the complex number.  
  IFloat abs() const;

  //! Overloaded complex sum 
  Rcomplex& operator+=(const Rcomplex& a);
  //! Overloaded real sum
  Rcomplex& operator+=(IFloat a);

  //! Overloaded complex subtraction
  Rcomplex& operator-=(const Rcomplex& a);
  //! Overloaded real subtraction
  Rcomplex& operator-=(IFloat a);

  //! Overloaded complex multiplication
  Rcomplex& operator*=(const Rcomplex& a);
  //! Overloaded real multiplication
  Rcomplex& operator*=(IFloat a);
    
  //! Overloaded complex division
  Rcomplex& operator/=(const Rcomplex& a);
  //! Overloaded real division
  Rcomplex& operator/=(IFloat a);

// These functions are specified as friends, but
// they do not actually reference private members.

  //! Complex conjugation
friend  Rcomplex conj(const Rcomplex& c);

  //! Complex addition
friend  Rcomplex operator + (const Rcomplex& c1, const Rcomplex& c2)
  {  return Rcomplex(c1) += c2; }

  //! Real addition 
friend  Rcomplex operator + (const Rcomplex& c, IFloat f)
  { return Rcomplex(c) += f; }
   
  //! Real addition 
friend  Rcomplex operator + (IFloat f, const Rcomplex& c)
  { return Rcomplex(c) += f; }
  
  //! Complex subtraction 
friend  Rcomplex operator - (const Rcomplex& c1, const Rcomplex& c2)
  {  return Rcomplex(c1) -= c2; }

  //! Real subtraction  
friend  Rcomplex operator - (const Rcomplex& c, IFloat f)
  {  return Rcomplex(c) -= f; }
   
  //! Real subtraction  
friend  Rcomplex operator - (IFloat f, const Rcomplex& c)
  {  return Rcomplex(f) -= c; }

  //! Complex multiplication
friend  Rcomplex operator * (const Rcomplex& c1, const Rcomplex& c2)
  {  return Rcomplex(c1) *= c2; }

  //! Real multiplication
friend  Rcomplex operator * (const Rcomplex& c, IFloat f)
  {  return Rcomplex(c) *= f; }

  //! Real multiplication 
friend  Rcomplex operator * (IFloat f, const Rcomplex& c)
  {  return Rcomplex(c) *= f; }

  //! Complex division
friend  Rcomplex operator / (const Rcomplex& c1, const Rcomplex& c2)
  {  return Rcomplex(c1) /= c2; }

  //! Real division 
friend  Rcomplex operator / (const Rcomplex& c, IFloat f)
  {  return Rcomplex(c) /= f; }
   
  //! Real division 
friend  Rcomplex operator / (IFloat f, const Rcomplex& c)
  {  return Rcomplex(f) /= c; }

  //! Unary negation 
friend  Rcomplex operator - (const Rcomplex& c);

// More friends for the purpose of injections. e.g.
// real(c1+c2) instead of tmp = c1+c2; tmp.real();

//! Real part
friend IFloat real(const Rcomplex& c) { return c.real(); }
//! Imaginary part
friend IFloat imag(const Rcomplex& c) { return c.imag(); }
//! Square norm
 friend IFloat norm(const Rcomplex& c) { return c.norm(); }
 //! Complex norm (modulus)
friend IFloat abs(const Rcomplex& c) { return c.abs(); }

};


#endif

CPS_END_NAMESPACE
