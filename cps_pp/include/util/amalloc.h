#include<config.h>
/*!\file
  \brief  Declaration of dynamic memory allocation routine for arrays.

  $Id: amalloc.h,v 1.2 2004-05-10 15:26:54 zs Exp $
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
