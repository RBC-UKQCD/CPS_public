#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc_common.C,v 1.9 2012-03-26 13:50:11 chulwoo Exp $
*/

#include <util/smalloc.h>

CPS_START_NAMESPACE

void* smalloc(const char cname[], const char fname[], const char vname[], size_t request){
    return smalloc(request, vname, fname, cname);
}

//void* smalloc(size_t request){
//    return smalloc(request, "", "smalloc", "");
//}

void sfree(const char cname[], const char fname[], const char vname[], void *p){
    sfree(p, cname, fname, vname);
}

void* fmalloc(const char *cname, const char *fname, const char *vname, size_t request){
    return fmalloc(request, vname, fname, cname);
}

void ffree(const char *cname, const char *fname, const char *vname, void *p){
    ffree(p, cname, fname, vname);
}

CPS_END_NAMESPACE
