#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of Double64 class/type.

  $Id: double64.h,v 1.5 2007-06-05 15:44:19 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2007-06-05 15:44:19 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/double64.h,v 1.5 2007-06-05 15:44:19 chulwoo Exp $
//  $Id: double64.h,v 1.5 2007-06-05 15:44:19 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/double64.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//
// WARNING WARNING WARNING
//-------------------------
//
// Do not use this file for other purposes except for the glb_sum
// If you want to use it for other applications you must 
// allocate memory and set the pointer  D64CRAM below.

/****************************************************************
 * Double64, written by Roy 					*
 ****************************************************************/

/* header file to support 64-bit (Double64) mathematical functions */

#ifndef INCLUDED_DOUBLE64_H_
#define INCLUDED_DOUBLE64_H_

//Get the global options:

//! A type used (solely?) to accumulate global sums in double precision.
/*!
  On platforms other than QCDSP, this should default to \c double.
  On QCDSP this is a custom class implementing 64-bit arithmetic.
*/


typedef GLOBALSUM_TYPE Double64;

#endif 

CPS_END_NAMESPACE
