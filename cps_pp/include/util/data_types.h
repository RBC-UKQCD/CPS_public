#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of basic data types.

  $Id: data_types.h,v 1.2 2003-07-24 16:53:53 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/data_types.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: data_types.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: data_types.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/data_types.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#ifndef INCLUDED_DATA_TYPES_H
#define INCLUDED_DATA_TYPES_H                //!< Prevent multiple inclusion

//------------------------------------------------------------------
// Global definitions:
//------------------------------------------------------------------
CPS_END_NAMESPACE
CPS_START_NAMESPACE

//------------------------------------------------------------------
//! Definition of 'Internal' floating point representation.
//------------------------------------------------------------------
typedef INTERNAL_LOCALCALC_TYPE IFloat; 

//------------------------------------------------------------------
// Definition of rfloat and rcomplex classes:
//------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/rfloat.h>
#include <util/rcomplex.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
//! Definition of general local floating point type.
//------------------------------------------------------------------
typedef LOCALCALC_TYPE Float;

//------------------------------------------------------------------
//! Definition of Complex type.
//------------------------------------------------------------------
typedef Rcomplex Complex;


#endif


CPS_END_NAMESPACE
