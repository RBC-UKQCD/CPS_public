#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of rfloat methods.

  $Id: rfloat.C,v 1.4 2004-08-18 11:58:08 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/comsrc/rfloat.C,v 1.4 2004-08-18 11:58:08 zs Exp $
//  $Id: rfloat.C,v 1.4 2004-08-18 11:58:08 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/comsrc/rfloat.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//  rfloat.C


CPS_END_NAMESPACE
#include <util/rfloat.h>
CPS_START_NAMESPACE


rfloat::rfloat(IFloat a) : x(a) {}

rfloat::rfloat(const rfloat& a) : x(a.x) {}

rfloat::~rfloat() {}


//-----------------------------------------------------------------
//  free friend functions
//-----------------------------------------------------------------

rfloat operator+(const rfloat& a, const rfloat& b)
	{ return rfloat(a.x) += b.x; }
rfloat operator+(double a, const rfloat& b)
	{ return rfloat(a) += b.x; }
rfloat operator+(const rfloat& a, double b)
	{ return rfloat(a.x) += b; }

rfloat operator-(const rfloat& a, const rfloat& b)
	{ return rfloat(a.x) -= b.x; }
rfloat operator-(double a, const rfloat& b)
	{ return rfloat(a) -= b.x; }
rfloat operator-(const rfloat& a, double b)
	{ return rfloat(a.x) -= b; }

rfloat operator*(const rfloat& a, const rfloat& b)
	{ return rfloat(a.x) *= b.x; }
rfloat operator*(double a, const rfloat& b)
	{ return rfloat(a) *= b.x; }
rfloat operator*(const rfloat& a, double b)
	{ return rfloat(a.x) *= b; }

rfloat operator/(const rfloat& a, const rfloat& b)
	{ return rfloat(a.x) /= b.x; }
rfloat operator/(double a, const rfloat& b)
	{ return rfloat(a) /= b.x; }
rfloat operator/(const rfloat& a, double b)
	{ return rfloat(a.x) /= b; }


CPS_END_NAMESPACE
