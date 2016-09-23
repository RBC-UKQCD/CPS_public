#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

  $Id: pmalloc.h,v 1.6 2004/10/27 14:30:25 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/10/27 14:30:25 $
//  $Header: /space/cvs/cps/cps++/include/util/pmalloc.h,v 1.6 2004/10/27 14:30:25 zs Exp $
// $Id: pmalloc.h,v 1.6 2004/10/27 14:30:25 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: pmalloc.h,v $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/include/util/pmalloc.h,v $
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
  Exits with the appropriate Error code if allocation fails.
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* pmalloc(size_t request);

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
