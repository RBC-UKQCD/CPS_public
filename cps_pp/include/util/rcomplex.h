#ifndef INCLUDED_RCOMPLEX_H
#define INCLUDED_RCOMPLEX_H

#if 1
#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of complex floating point data type.

*/
//--------------------------------------------------------------------
// rcomplex.h

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <math.h>
#include <complex>
CPS_START_NAMESPACE

// We use standard library for Rcomplex
typedef std::complex<IFloat> Rcomplex;

//#if __cplusplus <= 201103L
#if 0
static inline Rcomplex operator/(const Rcomplex &a, IFloat b)
{
    Rcomplex tmp = a;
    tmp /= b;
    return tmp;
}

static inline Rcomplex operator*(const Rcomplex &a, int b) {return a * Float(b);}
static inline Rcomplex operator*(int b, const Rcomplex &a) {return a * Float(b);}
#endif

CPS_END_NAMESPACE
#endif

#endif
