#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: fmalloc.C,v 1.7 2004-10-27 14:24:31 zs Exp $
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
    return p;
}

void ffree(void* p,
	   const char *vname, const char *fname, const char *cname){
    VRB.Sfree(cname, fname, vname, p);
    qfree(p);
}



CPS_END_NAMESPACE
