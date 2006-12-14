#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: fmalloc.C,v 1.2 2006-12-14 17:53:51 chulwoo Exp $
*/

#include <stdlib.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE

void* fmalloc(size_t request,
	      const char *vname, const char *fname, const char *cname){

    return smalloc(request, vname, fname, cname);
    
}

void* fmalloc(size_t request){
    return smalloc(request);
}

void ffree(void* p, const char *vname, const char *fname, const char *cname){
    sfree(p, cname, fname, vname);
}


CPS_END_NAMESPACE
