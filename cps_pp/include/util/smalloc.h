#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

  $Id: smalloc.h,v 1.5 2004-08-18 11:57:37 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:37 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/smalloc.h,v 1.5 2004-08-18 11:57:37 zs Exp $
//  $Id: smalloc.h,v 1.5 2004-08-18 11:57:37 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: smalloc.h,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/smalloc.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#ifndef _smalloc_h
#define _smalloc_h                //!< Prevent multiple inclusion


void* smalloc(int request);

void sfree(void* p);

void sclear(void);

void* fmalloc(int request);

void ffree(void* p);



#endif /* !_smalloc_h */

CPS_END_NAMESPACE
