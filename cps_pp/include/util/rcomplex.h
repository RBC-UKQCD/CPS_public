#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/rcomplex.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: rcomplex.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2001/08/16 10:50:30  anj
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
//  $RCSfile: rcomplex.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/rcomplex.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// rcomplex.h

#ifndef INCLUDED_RCOMPLEX_H
#define INCLUDED_RCOMPLEX_H

class Rcomplex;
CPS_END_NAMESPACE
#include<util/data_types.h>
CPS_START_NAMESPACE

class Rcomplex {
  IFloat re;
  IFloat im;

public:
  Rcomplex(IFloat a = 0, IFloat b = 0);
  Rcomplex(const Rcomplex& a);
  ~Rcomplex();

  Rcomplex& operator=(const Rcomplex& a);

  IFloat real() const { return re; }
  IFloat imag() const { return im; }

  void real(IFloat r) { re = r; }
  void imag(IFloat i) { im = i; }

  IFloat norm() const;	
  IFloat abs() const;

  Rcomplex& operator+=(const Rcomplex& a);
  Rcomplex& operator+=(IFloat a);

  Rcomplex& operator-=(const Rcomplex& a);
  Rcomplex& operator-=(IFloat a);

  Rcomplex& operator*=(const Rcomplex& a);
  Rcomplex& operator*=(IFloat a);
    
  Rcomplex& operator/=(const Rcomplex& a);
  Rcomplex& operator/=(IFloat a);

// These functions are specified as friends, but
// they do not actually reference private members.

friend  Rcomplex conj(const Rcomplex& c);

friend  Rcomplex operator + (const Rcomplex& c1, const Rcomplex& c2)
  {  return Rcomplex(c1) += c2; }

friend  Rcomplex operator + (const Rcomplex& c, IFloat f)
  { return Rcomplex(c) += f; }
   
friend  Rcomplex operator + (IFloat f, const Rcomplex& c)
  { return Rcomplex(c) += f; }
  
friend  Rcomplex operator - (const Rcomplex& c1, const Rcomplex& c2)
  {  return Rcomplex(c1) -= c2; }

friend  Rcomplex operator - (const Rcomplex& c, IFloat f)
  {  return Rcomplex(c) -= f; }
   
friend  Rcomplex operator - (IFloat f, const Rcomplex& c)
  {  return Rcomplex(f) -= c; }

friend  Rcomplex operator * (const Rcomplex& c1, const Rcomplex& c2)
  {  return Rcomplex(c1) *= c2; }

friend  Rcomplex operator * (const Rcomplex& c, IFloat f)
  {  return Rcomplex(c) *= f; }
   
friend  Rcomplex operator * (IFloat f, const Rcomplex& c)
  {  return Rcomplex(c) *= f; }

friend  Rcomplex operator / (const Rcomplex& c1, const Rcomplex& c2)
  {  return Rcomplex(c1) /= c2; }

friend  Rcomplex operator / (const Rcomplex& c, IFloat f)
  {  return Rcomplex(c) /= f; }
   
friend  Rcomplex operator / (IFloat f, const Rcomplex& c)
  {  return Rcomplex(f) /= c; }

friend  Rcomplex operator - (const Rcomplex& c);

// More friends for the purpose of injections. e.g.
// real(c1+c2) instead of tmp = c1+c2; tmp.real();
//
friend IFloat real(const Rcomplex& c) { return c.real(); }
friend IFloat imag(const Rcomplex& c) { return c.imag(); }
friend IFloat norm(const Rcomplex& c) { return c.norm(); }
friend IFloat abs(const Rcomplex& c) { return c.abs(); }

};


#endif
CPS_END_NAMESPACE
