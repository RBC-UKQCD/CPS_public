#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of Rcomplex methods,

  $Id: rcomplex.C,v 1.4 2004-08-18 11:58:07 zs Exp $ 
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:07 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rcomplex/comsrc/rcomplex.C,v 1.4 2004-08-18 11:58:07 zs Exp $
//  $Id: rcomplex.C,v 1.4 2004-08-18 11:58:07 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rcomplex/comsrc/rcomplex.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/rcomplex.h>
CPS_START_NAMESPACE

#ifdef _TARTAN
CPS_END_NAMESPACE
#include <math64.h>
CPS_START_NAMESPACE
#else
CPS_END_NAMESPACE
#include <math.h>
CPS_START_NAMESPACE
#endif

Rcomplex::Rcomplex(IFloat a, IFloat b)
: re(a), im(b){}

Rcomplex::Rcomplex(const Rcomplex& a)
: re(a.re), im(a.im){}

Rcomplex::~Rcomplex()
{}

Rcomplex& Rcomplex::operator=(const Rcomplex& a)
{ re = a.re; im = a.im; return *this; }

IFloat Rcomplex::abs() const { 
#ifdef _TARTAN
  return double( sqrt( norm() ) ); 
#else
  return sqrt(norm()); 
#endif

}

CPS_END_NAMESPACE
