#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of basic data types.

  $Id: data_types.h,v 1.5 2004-12-15 07:32:06 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-12-15 07:32:06 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/data_types.h,v 1.5 2004-12-15 07:32:06 chulwoo Exp $
//  $Id: data_types.h,v 1.5 2004-12-15 07:32:06 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: data_types.h,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/data_types.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#ifndef INCLUDED_DATA_TYPES_H
#define INCLUDED_DATA_TYPES_H                //!< Prevent multiple inclusion

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
#include <util/rfloat.h>
#include <util/rcomplex.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
//! Definition of general local floating point type.
//------------------------------------------------------------------
//typedef LOCALCALC_TYPE Float;

//------------------------------------------------------------------
//! Definition of Complex type.
//------------------------------------------------------------------
typedef Rcomplex Complex;


#endif


CPS_END_NAMESPACE
