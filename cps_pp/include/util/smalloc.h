#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

  $Id: smalloc.h,v 1.8 2004-09-07 05:20:52 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-09-07 05:20:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/smalloc.h,v 1.8 2004-09-07 05:20:52 chulwoo Exp $
//  $Id: smalloc.h,v 1.8 2004-09-07 05:20:52 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: smalloc.h,v $
//  $Revision: 1.8 $
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
void* smalloc(char *cname, char *fname, char *vname, int request);

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void sfree(void* p);
void sfree(char *cname, char *fname, char *vname, void *p);

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
