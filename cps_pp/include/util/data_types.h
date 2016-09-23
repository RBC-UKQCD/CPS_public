#ifndef INCLUDED_DATA_TYPES_H
#define INCLUDED_DATA_TYPES_H

#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of basic data types.

*/

//------------------------------------------------------------------
// Global definitions:
//------------------------------------------------------------------

//------------------------------------------------------------------
//! Definition of 'Internal' floating point representation.
//------------------------------------------------------------------
//typedef INTERNAL_LOCALCALC_TYPE IFloat; 

//------------------------------------------------------------------
// Definition of rfloat and rcomplex classes:
//------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/enum.h>
#include <util/rcomplex.h>
#include <complex>
CPS_START_NAMESPACE

//------------------------------------------------------------------
//! Definition of general local floating point type.
//------------------------------------------------------------------
//typedef LOCALCALC_TYPE Float;

//------------------------------------------------------------------
//! Definition of Complex type.
//------------------------------------------------------------------
//typedef Rcomplex Complex;
typedef std::complex<IFloat> Complex;

CPS_END_NAMESPACE
#endif
