#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rcomplex/noarch/rcomplex_sun.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: rcomplex_sun.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:39  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:36  anj
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
//  Revision 1.2  2001/05/25 06:16:10  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: rcomplex_sun.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rcomplex/noarch/rcomplex_sun.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<util/rcomplex.h>
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
