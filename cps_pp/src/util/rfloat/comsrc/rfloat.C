#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of rfloat methods.

  $Id: rfloat.C,v 1.3 2004-06-04 21:14:15 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:15 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/comsrc/rfloat.C,v 1.3 2004-06-04 21:14:15 chulwoo Exp $
//  $Id: rfloat.C,v 1.3 2004-06-04 21:14:15 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
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
