#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/data_types.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: data_types.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: data_types.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/data_types.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// data_types.h
//
// This files contains the typedef of the various data types.
//
//------------------------------------------------------------------


#ifndef INCLUDED_DATA_TYPES_H
#define INCLUDED_DATA_TYPES_H

//------------------------------------------------------------------
// Global definitions:
//------------------------------------------------------------------
CPS_END_NAMESPACE
#include<config.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Definition of 'Internal' floats (IFloat).
//------------------------------------------------------------------
typedef INTERNAL_LOCALCALC_TYPE IFloat;

//------------------------------------------------------------------
// Definition of rfloat and rcomplex classes:
//------------------------------------------------------------------
CPS_END_NAMESPACE
#include<util/rfloat.h>
#include<util/rcomplex.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Definition of general local floats (Float).
//------------------------------------------------------------------
typedef LOCALCALC_TYPE Float;

//------------------------------------------------------------------
// Definition of Complex.
//------------------------------------------------------------------
typedef Rcomplex Complex;


#endif

CPS_END_NAMESPACE
