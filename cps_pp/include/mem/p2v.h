#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:37 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/mem/p2v.h,v 1.4 2004-08-18 11:57:37 zs Exp $
//  $Id: p2v.h,v 1.4 2004-08-18 11:57:37 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/mem/p2v.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------//
//
// p2v.h
//
// Physical to virtual memory copy routines.
//
//------------------------------------------------------------------//


#ifndef INCLUDED_P2V_H
#define INCLUDED_P2V_H


//------------------------------------------------------------------//
// Copies into CRAM the Vector utilities assembly code
//------------------------------------------------------------------//
void p2vVector();


//------------------------------------------------------------------//
// Copies into CRAM the Wilson library 
//------------------------------------------------------------------//
void p2vWilsonLib();


//------------------------------------------------------------------//
// Copies into CRAM the Staggered dirac_serial assembly code
//------------------------------------------------------------------//
void p2vStagDs();


//------------------------------------------------------------------//
// Copies into CRAM the Gauge heat bath (ghb) assembly code
//------------------------------------------------------------------//
void p2vGhb();


//------------------------------------------------------------------//
// Copies into CRAM the Matrix Multiplication assembly code
//------------------------------------------------------------------//
void p2vCloverLib();

#endif


CPS_END_NAMESPACE
