#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of Rcomplex methods,

  $Id: rcomplex_sun.C,v 1.4 2004-08-18 11:58:08 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rcomplex/noarch/rcomplex_sun.C,v 1.4 2004-08-18 11:58:08 zs Exp $
//  $Id: rcomplex_sun.C,v 1.4 2004-08-18 11:58:08 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rcomplex/noarch/rcomplex_sun.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/rcomplex.h>
CPS_START_NAMESPACE

IFloat Rcomplex::norm() const
{  return re*re + im*im; }
 
Rcomplex& Rcomplex::operator+=(const Rcomplex& a)
{  re += a.re;
   im += a.im;
   return *this;
}
 
Rcomplex& Rcomplex::operator+=(IFloat a)
{  re += a;
   return *this;
}

Rcomplex& Rcomplex::operator-=(const Rcomplex& a)
{  re -= a.re;
   im -= a.im;
   return *this;
}
 
Rcomplex& Rcomplex::operator-=(IFloat a)
{  re -= a;
   return *this;
}

Rcomplex& Rcomplex::operator*=(const Rcomplex& a)
{ 
   IFloat sre = re * a.re - im * a.im;
   im = re * a.im + im * a.re;
   re = sre;
   return *this;
}
 
Rcomplex& Rcomplex::operator*=(IFloat a)
{  re *= a;
   im *= a;
   return *this;
}
 
Rcomplex& Rcomplex::operator/=(const Rcomplex& a)
{
  IFloat norm2_a_1 = 1.0/(a.re * a.re + a.im * a.im);
  IFloat sre = (re * a.re + im * a.im) * norm2_a_1;
  im = (im * a.re - re * a.im) * norm2_a_1;
  re = sre;
  return *this;
}

Rcomplex& Rcomplex::operator/=(IFloat a)
{  re /= a;
   im /= a;
   return *this;
}

Rcomplex conj(const Rcomplex& c)
{  return Rcomplex(c.re, -c.im); }

Rcomplex operator - (const Rcomplex& c)
{  return Rcomplex(-c.re, -c.im); }


CPS_END_NAMESPACE
