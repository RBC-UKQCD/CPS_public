#ifndef INCLUDED_RCOMPLEX_H
#define INCLUDED_RCOMPLEX_H

#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of complex floating point data type.

  $Id: rcomplex.h,v 1.6 2012-08-10 14:05:33 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-10 14:05:33 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/rcomplex.h,v 1.6 2012-08-10 14:05:33 chulwoo Exp $
//  $Id: rcomplex.h,v 1.6 2012-08-10 14:05:33 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: rcomplex.h,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/rcomplex.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// rcomplex.h

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <math.h>
#include <complex>
CPS_START_NAMESPACE

// We use standard library for Rcomplex
typedef std::complex<IFloat> Rcomplex;

static inline Rcomplex operator/(const Rcomplex &a, IFloat b)
{
    Rcomplex tmp = a;
    tmp /= b;
    return tmp;
}

static inline Rcomplex operator*(const Rcomplex &a, int b) {return a * Float(b);}
static inline Rcomplex operator*(int b, const Rcomplex &a) {return a * Float(b);}

CPS_END_NAMESPACE

#endif
