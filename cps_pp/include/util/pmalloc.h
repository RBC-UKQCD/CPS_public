#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

  $Id: pmalloc.h,v 1.3 2004-07-01 17:43:40 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-07-01 17:43:40 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pmalloc.h,v 1.3 2004-07-01 17:43:40 chulwoo Exp $
// $Id: pmalloc.h,v 1.3 2004-07-01 17:43:40 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pmalloc.h,v $
//  $Revision: 1.3 $
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
