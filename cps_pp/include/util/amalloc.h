#include<config.h>
/*!\file
  \brief  Declaration of dynamic memory allocation routine for arrays.

  $Id: amalloc.h,v 1.3 2004-08-18 11:57:37 zs Exp $
*/
//--------------------------------------------------------------------
#ifndef AMALLOC_H
#define AMALLOC_H                //!< Prevent multiple inclusion

#include <stdarg.h>
#include <stddef.h>

CPS_START_NAMESPACE

void* amalloc(size_t, int, ...);

CPS_END_NAMESPACE

#endif 
