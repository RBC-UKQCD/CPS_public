#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

  $Id: pmalloc.h,v 1.5 2004-09-02 16:57:04 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-09-02 16:57:04 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pmalloc.h,v 1.5 2004-09-02 16:57:04 zs Exp $
// $Id: pmalloc.h,v 1.5 2004-09-02 16:57:04 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pmalloc.h,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pmalloc.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#ifndef _pmalloc_h
#define _pmalloc_h             //!< Prevent multiple inclusion

/*!\defgroup mem_alloc Memory allocation */
/*! \addtogroup mem_alloc 
  @{
*/

//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* pmalloc(int request);

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void pfree(void* p);

//! Doesn't appear to do anything.
void pclear(void);

/*! @} */

#endif /* !_pmalloc_h */

CPS_END_NAMESPACE
