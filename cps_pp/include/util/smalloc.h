#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

  $Id: smalloc.h,v 1.3 2004-07-01 17:43:40 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-07-01 17:43:40 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/smalloc.h,v 1.3 2004-07-01 17:43:40 chulwoo Exp $
//  $Id: smalloc.h,v 1.3 2004-07-01 17:43:40 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: smalloc.h,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/smalloc.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#ifndef _smalloc_h
#define _smalloc_h                //!< Prevent multiple inclusion


void* smalloc(int request);

void sfree(void* p);

void sclear(void);
#if TARGET == QCDOC
void* fmalloc(int request);

void ffree(void* p);
#endif


#endif /* !_smalloc_h */

CPS_END_NAMESPACE
