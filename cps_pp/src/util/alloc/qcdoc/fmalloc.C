#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: fmalloc.C,v 1.10 2005-06-16 07:23:12 chulwoo Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <qcdoc_align.h>
#include <qalloc.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE

void* fmalloc(size_t request,
	      const char *vname, const char *fname, const char *cname){
    void *p = qalloc(QFAST,request);
    if (!p) return smalloc(request, cname, fname, vname);
    VRB.Smalloc(cname, fname, vname, p, request);
    return p;
}

void* fmalloc(size_t request){
    return fmalloc(request, "", "fmalloc", "");
}


void ffree(void* p,
	   const char *vname, const char *fname, const char *cname){
    sfree(p, cname, fname, vname);
}



CPS_END_NAMESPACE
