#include<config.h>
/*!\file
  \brief  Declaration of dynamic memory allocation routine for arrays.

  $Id: amalloc.h,v 1.4 2004-10-27 14:30:25 zs Exp $
*/
//--------------------------------------------------------------------
#ifndef AMALLOC_H
#define AMALLOC_H                //!< Prevent multiple inclusion


CPS_START_NAMESPACE

void* amalloc(void*  (*allocator)(size_t, const char *vname="",
			  const char *fname="smalloc", const char *cname=""),
	      size_t, int, ...);

CPS_END_NAMESPACE

#endif 
