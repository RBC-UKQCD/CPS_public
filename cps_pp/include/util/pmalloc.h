#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

  $Id: pmalloc.h,v 1.4 2004-08-18 11:57:37 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:37 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pmalloc.h,v 1.4 2004-08-18 11:57:37 zs Exp $
// $Id: pmalloc.h,v 1.4 2004-08-18 11:57:37 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pmalloc.h,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pmalloc.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#ifndef _pmalloc_h
#define _pmalloc_h             //!< Prevent multiple inclusion


void* pmalloc(int request);

void pfree(void* p);

void pclear(void);


#endif /* !_pmalloc_h */

CPS_END_NAMESPACE
