#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc_common.C,v 1.5 2004-10-27 14:24:31 zs Exp $
*/

#include <util/smalloc.h>

CPS_START_NAMESPACE

void* smalloc(const char *cname, const char *fname, const char *vname, size_t request){
    return smalloc(request, vname, fname, cname);
}

void sfree(const char *cname, const char *fname, const char *vname, void *p){
    return sfree(cname, fname, vname, p);
}

void* fmalloc(const char *cname, const char *fname, const char *vname, size_t request){
  return fmalloc(cname, fname, vname, request);
}

void ffree(const char *cname, const char *fname, const char *vname, void *p){
    return ffree(cname, fname, vname, p);
}

CPS_END_NAMESPACE
