#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc_common.C,v 1.6 2004-11-01 10:26:28 zs Exp $
*/

#include <util/smalloc.h>

CPS_START_NAMESPACE

void* smalloc(const char *cname, const char *fname, const char *vname, size_t request){
    return smalloc(request, vname, fname, cname);
}

void sfree(const char *cname, const char *fname, const char *vname, void *p){
    sfree(p, cname, fname, vname);
}

void* fmalloc(const char *cname, const char *fname, const char *vname, size_t request){
    return fmalloc(request, cname, fname, vname);
}

void ffree(const char *cname, const char *fname, const char *vname, void *p){
    ffree(p, cname, fname, vname);
}

CPS_END_NAMESPACE
