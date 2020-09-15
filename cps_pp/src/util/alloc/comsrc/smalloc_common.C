#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc_common.C,v 1.9 2012/03/26 13:50:11 chulwoo Exp $
*/

#include <util/smalloc.h>

CPS_START_NAMESPACE

#if 0
void* smalloc(size_t request, const char vname[], const char fname[], const char cname[]){
    return smalloc(cname, fname, vname, request);
}
void sfree(void *p, const char vname[], const char fname[], const char cname[]){
    sfree(cname, fname, vname, p);
}

void* fmalloc(size_t request, const char vname[], const char fname[], const char cname[]){
    return fmalloc(cname, fname, vname, request);
}
void ffree(void *p, const char vname[], const char fname[], const char cname[]){
    ffree(cname, fname, vname, p);
}
#endif

CPS_END_NAMESPACE
