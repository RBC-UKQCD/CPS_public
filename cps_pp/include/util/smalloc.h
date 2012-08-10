#ifndef _smalloc_h
#define _smalloc_h                //!< Prevent multiple inclusion

#include<config.h>
#include <stdlib.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

  $Id: smalloc.h,v 1.14 2012-08-10 14:05:33 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-10 14:05:33 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/smalloc.h,v 1.14 2012-08-10 14:05:33 chulwoo Exp $
//  $Id: smalloc.h,v 1.14 2012-08-10 14:05:33 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: smalloc.h,v $
//  $Revision: 1.14 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/smalloc.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


/*! \addtogroup mem_alloc Memory allocation
  @{
*/


//! Allocate memory
/*!
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.
  \param request The amount of memory (in bytes) to allocate
  \param vname The name of the variable pointing to the allocated memory.
  \param fname The name of the calling function or method
  \param cname The name of the calling class
  \return A pointer to the allocated memory
  \post Program exits with the appropriate Error code if allocation fails.  
*/
void* smalloc(size_t request,
	      const char vname[], const char fname[]="smalloc", const char cname[]="");
void* smalloc(size_t request);

//! Allocate memory
/*!
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.
  \param cname The name of the calling class
  \param fname The name of the calling function or method
  \param vname The name of the variable pointing to the allocated memory.
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
  \post  Exits with the appropriate Error code if allocation fails.
*/

void* smalloc(const char cname[], const char fname[], const char vname[], size_t request);

//! Free allocate memory
/*!
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.
  \param p Pointer to the memory to be freed.
  \param vname The name of the variable pointing to the allocated memory.
  \param fname The name of the calling function or method
  \param cname The name of the calling class
  \post  Exits with the appropriate Error code if allocation fails.  
*/
void sfree(void* p,
	   const char vname[], const char fname[]="sfree", const char cname[]="");
void sfree(void* p);

//! Free memory
/*!
  With verbose reporting of details, if the appropriate Verbose
  level is enabled.
  \param cname The name of the calling class
  \param fname The name of the calling function or method
  \param vname The name of the variable pointing to the allocated memory.
  \param p Pointer to the memory to be freed.
*/
void sfree(const char cname[], const char fname[], const char vname[], void *p);

//! Doesn't appear to do anything.
void sclear();


//! Allocate memory
/*!
  Allocates in fast memory (EDRAM) on the QCDOC.
  If this fails, then allocation of transient DDR memory is attempted.
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.
  \param request The amount of memory (in bytes) to allocate
  \param vname The name of the variable pointing to the allocated memory.
  \param fname The name of the calling function or method
  \param cname The name of the calling class
  \return A pointer to the allocated memory
  \post  Exits with the appropriate Error code if allocation fails.  
*/
void* fmalloc(size_t request,
	      const char vname[], const char fname[]="fmalloc", const char cname[]="");
void* fmalloc(size_t request);

//! Free allocate memory
/*!
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.  
  \param p Pointer to the memory to be freed.
  \param vname The name of the variable pointing to the allocated memory.
  \param fname The name of the calling function or method
  \param cname The name of the calling class
*/
void ffree(void* p,
	   const char vname[], const char fname[]="ffree", const char cname[]="");
void ffree(const char cname[], const char fname[], const char vname[], void *p);
void ffree(void* p);

//! Allocate memory
/*!
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.
  \param cname The name of the calling class
  \param fname The name of the calling function or method
  \param vname The name of the variable pointing to the allocated memory.
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
  \post  Exits with the appropriate Error code if allocation fails.
*/
void* fmalloc(const char cname[], const char fname[], const char vname[], size_t request);

//! Free memory
/*!
  With verbose reporting of details, if the appropriate Verbose
  level is enabled.
  \param cname The name of the calling class
  \param fname The name of the calling function or method
  \param vname The name of the variable pointing to the allocated memory.
  \param p Pointer to the memory to be freed.
*/

/*! @} */


CPS_END_NAMESPACE

#endif /* !_smalloc_h */
