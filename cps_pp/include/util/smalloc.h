#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

  $Id: smalloc.h,v 1.6 2004-09-02 16:57:05 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-09-02 16:57:05 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/smalloc.h,v 1.6 2004-09-02 16:57:05 zs Exp $
//  $Id: smalloc.h,v 1.6 2004-09-02 16:57:05 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: smalloc.h,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/smalloc.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#ifndef _smalloc_h
#define _smalloc_h                //!< Prevent multiple inclusion

/*! \addtogroup mem_alloc Memory allocation
  @{
*/


//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* smalloc(int request);

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void sfree(void* p);

//! Doesn't appear to do anything.
void sclear(void);


//! Allocate memory
/*!
  Allocates fast memory (EDRAM) on the QCDOC.
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* fmalloc(int request);

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void ffree(void* p);

/*! @} */

#endif /* !_smalloc_h */

CPS_END_NAMESPACE
